#http://bioconductor.org/packages/devel/bioc/vignettes/ggbio/inst/doc/ggbio.pdf

#you may have to install bunch of this stuff
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggbio")
BiocManager::install("Homo.sapiens")
BiocManager::install("karyoploteR")


library(karyoploteR)
library(ggbio)
library(Homo.sapiens)
class(Homo.sapiens)
library(tidyverse)
library(readxl)

#these are all the significant genes
sigGenes = read_excel("sig_genes.xlsx")

#the bonforoni corrected data
brca = read_tsv('brca_CNV_up.tsv')
ccrcc = read_tsv('ccrcc_CNV_up.tsv')
endo = read_tsv('en_CNV_up.tsv')
gbm = read_tsv('gbm_CNV_up.tsv')
hnscc = read_tsv('hnscc_CNV_up.tsv')
luad = read_tsv('luad_CNV_up.tsv')
ovarian = read_tsv('ovarian_CNV_up.tsv')

#the bonforoni corrected data with the CNV/Protien pairs
brca2 = mutate(brca, name = paste(CNV, Protien, sep = '/'))
ccrcc2 = mutate(ccrcc,name = paste(CNV, Protien, sep = '/'))
endo2 = mutate(endo, name = paste(CNV, Protien, sep = '/'))
gbm2 = mutate(gbm, name = paste(CNV, Protien, sep = '/'))
hnscc2 = mutate(hnscc, name = paste(CNV, Protien, sep = '/'))
luad2 = mutate(luad, name = paste(CNV, Protien, sep = '/'))
ovarian2 = mutate(ovarian, name = paste(CNV, Protien, sep = '/'))

#list of CNV/Protien pairs that are common between the cancers (excluding gbm)
common = Reduce(intersect, list(brca2$name,ccrcc2$name, endo2$name, hnscc2$name, luad2$name, ovarian2$name))

nCommon = Reduce(intersect,list(luad$CNV, ovarian$CNV))

#myList is a list of CNV from common
com = str_split(common, '/')
myList = unlist(com)
myList = myList[c(TRUE,FALSE)]


#takes in list of CNVs, return ones that work with "genesymbol":
getList = function(inList){
  returnable = c()
  numFilteredOut = 0
  for(i in c(1:length(inList))){
    worked = try(genesymbol[inList[i]])
    if(isTRUE(class(worked)=="try-error")) {
      print(i)
      #print(inList[i])
      numFilteredOut = numFilteredOut+1
      next } 
    else {
      returnable = c(returnable, inList[i])
    } 
  }
  print(paste0("number of CNVs filtered out: ", numFilteredOut, " of ", length(inList), " (",round((numFilteredOut/length(inList))*100, 3),"%)"))
  returnable
}

#plot list of genes onto chromosomes
#usage: plotList(brca$CNV), plotList(sigGenes$BRCA)
plotList = function(inList){
 
  data(genesymbol, package = "biovizBase")
  tryList = getList(inList)
  wh <- genesymbol[tryList] 
  wh <- range(wh, ignore.strand = TRUE)
  kp <- plotKaryotype(genome="hg19",plot.type=2, 
                      chromosomes=c("chr1", "chr3", "chr7", "chr8", "chr20"))
  kpPlotRegions(kp, wh, col="green")
}









