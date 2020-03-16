require(ggplot2)
require(gridExtra)
library(tidyverse)
library(grid)
all_data <- read_csv('all_data.csv')
up_or_down <- read_csv('up_or_down.csv')
cutoffs <- read_csv('cutoffs_table.csv')
brca_counts <- read_csv('brca_counts.csv')

# Plot 1 
all_data = all_data %>%
  as_tibble() %>%
  gather(key="Y", value="Z", -1) %>%
  group_by(X1, Z) %>%
  summarise(
    counts = n()
  )

plt1 = ggplot(total, aes(fill=factor(Z), y=counts, x=X1)) +
  geom_bar(position="dodge", stat="identity") +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  labs(fill="Copy Number", title = 'B') +
  scale_fill_discrete(breaks=c("high", "low", "normal"), labels=c("Significant Increase", "Significant Decrease", "Average"))


# Plot 2
plt2 = ggplot(brca_counts, aes(high)) +
  geom_histogram() +
  geom_vline(xintercept=54.5, color='red') +
  labs(title="A")


# Plot 3
plt3 = up_or_down %>%
  as_tibble() %>%
  gather(key="Y", value="Z", -1, -9) %>%
  group_by(chromosome) %>%
  ggplot(aes(X1, factor(Y), fill=factor(Z))) +
  geom_tile(mapping=aes(group=factor(chromosome))) +
  theme(axis.text.x=element_blank()) +
  xlab("Gene") +
  ylab("Cancer Type") +
  labs(fill="Copy Number", title="C") +
  scale_fill_discrete(breaks=c("high", "low", "normal"), labels=c("Significant Increase", "Significant Decrease", "Average"))


g = grid.arrange(
  plt2, plt1, plt3,
  widths = c(1,1),
  layout_matrix = rbind(c(1, 2),
                        c(3, 3)),
  bottom=textGrob("Figure 1: Similarities in Copy Number Variation among cancer types. A) Histogram showing a count of genes 
above the given threshold of 0.2 in Breast Cancer. In red is the calculated threshold for breast cancer (n=54.5). 
Any gene with >54.5 tumor samples above the threshold is considered to have a significant increase in copy 
number variation. B) Comparison of number of significant copy number increase and decreases in each 
cancer type. C) A Heatmap showing the patterns of increased and decreased copy number across cancer types
in 394 genes. Genes along the x-axis are ordered alphabetically.", 
                  x=0.05, 
                  y=0.5, 
                  just="left")
)

ggsave("figure1_final.pdf", g)
