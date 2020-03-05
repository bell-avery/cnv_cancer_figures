import cptac
import pandas as pd
import statsmodels.api as sm
import time
import concurrent.futures
import multiprocessing
import argparse
import scipy.stats

parser = argparse.ArgumentParser()
parser.add_argument("CNVfile", help = "The CNVfile")
parser.add_argument("ProteinFile", help = "The ProteinFile")
parser.add_argument("-N","--NumFromCNV", help = "The number of columns from CNV file", type=int, default = 1)
parser.add_argument("-S","--StartingIndex", help = "The starting index)", type=int, default = 0)

args = parser.parse_args()

CNVfile = args.CNVfile
proteinFile = args.ProteinFile
numCol = args.NumFromCNV
startIndex = args.StartingIndex


CNVsub = pd.read_csv(CNVfile, sep = '\t', header=[0], index_col = 0)
proSub = pd.read_csv(proteinFile, sep = '\t', header = [0], index_col = 0)
CNVsub = CNVsub.iloc[:,startIndex:startIndex+numCol]


def linRegress(cnvName, cnvCol):
    regressedResults = []
    for proName, proCol in proSub.iteritems():
        tableMerged = pd.merge(cnvCol, proCol, how = 'inner', right_index = True, left_index = True)
        independent = tableMerged.iloc[:,0]
        dependent = tableMerged.iloc[:,1]
        independent = independent.astype('float64')
        dependent = dependent.astype('float64')
        print(type(independent[0]))
        print(independent[0])
        print(type(dependent[0]))
        print(dependent[0])
        slope, intercept, rValue, pValue, stdErr = scipy.stats.linregress(dependent, independent)
        if pValue < .05:
            regressedResults.append([cnvName, proName, slope, intercept, rValue, pValue, stdErr])
        #model = sm.OLS(dependent, independent).fit()
        #if model.pvalues[0] < .05 and model.rsquared > .3:
            #returnable = [cnvName, proName, model.pvalues[0], model.rsquared]
            #regressedResults.append(returnable)
    return regressedResults


startTime = time.time()

with concurrent.futures.ProcessPoolExecutor() as executor:
        output = []
        results = [executor.submit(linRegress, cnvName, cnvCol) for cnvName, cnvCol in CNVsub.iteritems()]

        for f in concurrent.futures.as_completed(results):
                output = output + f.result()
        tableToWrite = pd.DataFrame.from_records(output, columns = ['CNV', 'Protien', 'slope', 'intercept', 'rValue', 'pValue', 'stdErr'])
        outName = ''.join([CNVfile,'linRegressOutput.tsv'])