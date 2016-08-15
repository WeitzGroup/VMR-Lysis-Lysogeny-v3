# Part of the code used in:
# Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
# 
# From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny-v2
# MIT License


##  R file for analysis of Knowles et al. 2016 virome hallmark genes in figures 4a,b and S4a


### Dependencies ###


library(resample)  # required for bootstrap sampling
# Tim Hesterberg (2015). resample: Resampling Functions. R package
#  version 0.4. https://CRAN.R-project.org/package=resample

library(MASS)  # requried for robust regression
# Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with
#  S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0

library(robustbase) # required for robust regression
#Peter Rousseeuw, Christophe Croux, Valentin Todorov, Andreas
#  Ruckstuhl, Matias Salibian-Barrera, Tobias Verbeke, Manuel Koller,
#  Martin Maechler (2016). robustbase: Basic Robust Statistics. R
#  package version 0.92-6. URL
#  http://CRAN.R-project.org/package=robustbase



### Import data ###


# This is downloaded from Knowles et al. 2016 (see "Figure 4: Temperate features in viromes increase with host density." from http://dx.doi.org/10.1038/nature17193)
	
	data4 <- read.csv("nature17193-f4.csv", header=TRUE)


# This is downloaded from the SI in Knowles et al. 2016 (see "Extended Data Figure 4: Temperateness of viral communities increases with host density and viral functional composition change." from http://dx.doi.org/10.1038/nature17193)

	dataSI4a <- read.csv("nature17193-sf4.csv", header=TRUE)
	
### Prepare data for analysis ###


	DATA = as.data.frame(cbind(data4$Microbes.per.ml,data4$X..integrase.reads,data4$X..excisionase.reads,data4$Functional.Diversity..H..,data4$X..pathogenicity.reads,dataSI4a$X..prophage.reads))
	COLNAMES = c("Cells","Integrase","Excisionase","FuncDiv","Pathogenicity","PProphage")
	colnames(DATA) = COLNAMES


### Reanalyse confidence intervals for data ###

# Samples
	SAMPLES <- 10000 #want to use 10,000 bootstrap samples

# load the function to perform analysis
	source("CItesting.R")

# set seed for reproducibility
	SEED = 1

#Check log10(Microbial Cells) ~ Provirus like reads -- output is saved to PProphage.Rdata and presented in PProphage.csv
print("Checking log10(Microbial Cells) ~ Provirus like reads")
DATA$Interest = DATA$PProphage
CItesting(DATA,SEED,SAMPLES,"PProphage")

#Check log10(Microbial Cells) ~ Integrase -- output is saved to Integrase.Rdata and presented in Integrase.csv
print("Checking log10(Microbial Cells) ~ Integrase")
DATA$Interest = DATA$Integrase
CItesting(DATA,SEED,SAMPLES,"Integrase")

#Check log10(Microbial Cells) ~ Excisionase -- output is saved to Excisionase.Rdata and presented in Excisionase.csv
print("Checking log10(Microbial Cells) ~ Excisionase")
DATA$Interest = DATA$Excisionase
CItesting(DATA,SEED,SAMPLES,"Excisionase")

