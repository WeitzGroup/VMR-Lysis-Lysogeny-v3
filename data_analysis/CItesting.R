# Part of the code used in:
# Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
# 
# From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny-v3
# MIT License


CItesting <- function(DATA,SEED,SAMPLES,TextInterest) {

### Create bootstrap ensembles for each metric that is considered ###
## Standard approaches

#1) Pearson
	print("...performing Pearson bootstrap...")
	
	set.seed(SEED) # for reproductibility
	PearsonBS = boot(DATA, FPearson, SAMPLES) #perform a 10,000 sample bootstrap


#2) Kendall
	print("...performing Kendall bootstrap...")

	set.seed(SEED) # for reproductibility
	KendallBS = boot(DATA, FKendall, SAMPLES) #perform a 10,000 sample bootstrap
		
#3) Spearman
	print("...performing Spearman bootstrap...")

	set.seed(SEED) # for reproductibility
	SpearmanBS = boot(DATA, FSpearman, SAMPLES) #perform a 10,000 sample bootstrap


## Robust regression techniques using the MASS package

#1) Tukey bisquare
	print("...performing Bisquare M (MASS) bootstrap...")
	
	set.seed(SEED) # for reproductibility
	TbisBS = boot(DATA, FTbisquare, SAMPLES) #perform a 10,000 sample bootstrap
	

#2) Tukey bisquare -MM
	print("...performing Bisquare MM (MASS) bootstrap...")

	set.seed(SEED) # for reproductibility
	TbisBSmm = boot(DATA, FTbisquare_mm, SAMPLES) #perform a 10,000 sample bootstrap


#3) Huber
	print("...performing Huber M (MASS) bootstrap...")

	set.seed(SEED) # for reproductibility
	HuberBS = boot(DATA, FHuber, SAMPLES) #perform a 10,000 sample bootstrap

#4) Hampel
	print("...performing Hampel M (MASS) bootstrap...")

	set.seed(SEED) # for reproductibility
	HampelBS = boot(DATA, FHampel, SAMPLES) #perform a 10,000 sample bootstrap
	

## Robust regression techniques using the robustbase package

#1) Bisquare	
	print("...performing Bisquare MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility
	n_BisquareBS = boot(DATA, newBisquare, SAMPLES)


#2) GGW
	print("...performing GWW MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility	
	n_GGWBS = boot(DATA, newGGW, SAMPLES)


#3) LQQ
	print("...performing LQQ MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility
	n_lqqBS = boot(DATA, newlqq, SAMPLES)


#4) Hampel
	print("...performing Hampel MM (robustbase) bootstrap...")
	
	set.seed(SEED) # for reproductibility
	n_HampelBS = boot(DATA, newHampel, SAMPLES)


#5) Optimal
	print("...performing Optimal MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility
	n_OptBS = boot(DATA, newOptimal, SAMPLES)


#6) Welsh
	print("...performing Welsh MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility
	n_welshBS = boot(DATA, newWelsh, SAMPLES)	

		

### Analyse the bootstraps for each metric ###


	names = c("Pearson","Kendall","Spearman","M Bisquare","MM Bisquare (rlm)","M Hampel","M Huber","MM Bisquare(lmrob)","MM GGW","MM LQQ","MM Hampel","MM Optimal","MM Welsh")


## Find:
# a) The correlation coefficient / slope coefficient associated with each analysis
# b) The p-value associated with each correlation coefficient/ slope coefficient
# c) for the robust regression analyses, the number of datapoints which are assigned a weighting of zero.	

#correlation:
	X = cor.test(log10(DATA$Cells), DATA$Interest, method="pearson")
	OC1 = X$estimate	
	OP1 = X$p.value
	OZ1 = 0 # not a robust regression analysis. Data points are not assigned weights here.
	

	X = cor.test(log10(DATA$Cells), DATA$Interest, method="kendall")
	OC2 = X$estimate
	OP2 = X$p.value
	OZ2 = 0 # not a robust regression analysis. Data points are not assigned weights here.
	

	X = cor.test(log10(DATA$Cells), DATA$Interest, method="spearman")
	OC3 = X$estimate
	OP3 = X$p.value
	OZ3 = 0 # not a robust regression analysis. Data points are not assigned weights here.
	
#rlm
	X = rlm(Interest ~ log10(Cells), DATA, psi="psi.bisquare",method="M", maxit=1000)
	AA = data.frame(summary(X)$coefficients)$t.value # t-values
	BC1 = X$coefficients[2] #slope coefficient
	BP1 = 2*pt(abs(AA[2]), summary(X)$df[2], lower.tail=FALSE)#p-value from t-value for slope coefficient
	BZ1 = sum(X$w==0)# number of datapoints with weight 0
	
	
	X = rlm(Interest ~ log10(Cells), DATA, psi="psi.bisquare",method="MM", maxit=1000)
	AA = data.frame(summary(X)$coefficients)$t.value # t-values
	BC2 = X$coefficients[2] #slope coefficient
	BP2 = 2*pt(abs(AA[2]), summary(X)$df[2], lower.tail=FALSE) #p-value from t-value for slope coefficient
	BZ2 = sum(X$w==0) # number of datapoints with weight 0
	
	
	
	X = rlm(Interest ~ log10(Cells), DATA, psi="psi.hampel",method="M", maxit=1000)
	AA = data.frame(summary(X)$coefficients)$t.value # t-values
	BC3 = X$coefficients[2] #slope coefficient
	BP3 = 2*pt(abs(AA[2]), summary(X)$df[2], lower.tail=FALSE)#p-value from t-value for slope coefficient
	BZ3 = sum(X$w==0)# number of datapoints with weight 0
	
		

	X = rlm(Interest ~ log10(Cells), DATA, psi="psi.huber",method="M", maxit=1000)
	AA = data.frame(summary(X)$coefficients)$t.value # t-values
	BC4 = X$coefficients[2] #slope coefficient
	BP4 = 2*pt(abs(AA[2]), summary(X)$df[2], lower.tail=FALSE)#p-value from t-value for slope coefficient
	BZ4 = sum(X$w==0)# number of datapoints with weight 0
	
	
	

#lmrob
	X = lmrob(Interest ~ log10(Cells), DATA, psi="bisquare",method="MM")
	AC1 = X$coefficients[2] #slope coefficient
	AP1 = data.frame(summary(X)$coefficients)$Pr...t..[2] # p-value of slope coefficient
	AZ1 = sum(X$rweights==0) # number of datapoints with weight 0
	
	
	X = lmrob(Interest ~ log10(Cells), DATA, psi="GGW",method="MM")
	AC2 = X$coefficients[2] #slope coefficient
	AP2 = data.frame(summary(X)$coefficients)$Pr...t..[2] # p-value of slope coefficient
	AZ2 = sum(X$rweights==0) # number of datapoints with weight 0
		

	X = lmrob(Interest ~ log10(Cells), DATA, psi="LQQ",method="MM")
	AC3 = X$coefficients[2] #slope coefficient
	AP3 = data.frame(summary(X)$coefficients)$Pr...t..[2] # p-value of slope coefficient
	AZ3 = sum(X$rweights==0) # number of datapoints with weight 0
	

	X = lmrob(Interest ~ log10(Cells), DATA, psi="Hampel",method="MM")
	AC4 = X$coefficients[2] #slope coefficient
	AP4 = data.frame(summary(X)$coefficients)$Pr...t..[2] # p-value of slope coefficient
	AZ4 = sum(X$rweights==0) # number of datapoints with weight 0
	
	
	X = lmrob(Interest ~ log10(Cells), DATA, psi="Optimal",method="MM")
	AC5 = X$coefficients[2] #slope coefficient
	AP5 = data.frame(summary(X)$coefficients)$Pr...t..[2] # p-value of slope coefficient
	AZ5 = sum(X$rweights==0) # number of datapoints with weight 0
		

	X = lmrob(Interest ~ log10(Cells), DATA, psi="Welsh",method="MM")
	AC6 = X$coefficients[2] #slope coefficient
	AP6 = data.frame(summary(X)$coefficients)$Pr...t..[2] # p-value of slope coefficient
	AZ6 = sum(X$rweights==0) # number of datapoints with weight 0
	
	
## Find the 95% confidence intervals ##
	#BOOTway
	GetCIs <- function(x) {
		NN = x
		if(dim(x$t)[2]>1){ #if robust regression need to watch out for non-convergence.
		NN$t0 = x$t0[2] # just the slope as observation
		NN$t = matrix(x$t[which(x$t[,3]==1),2]) # remove non-converged replicates. (only look at slope replicates)
		NN$R = length(NN$t) #change number of replicates; as non-converged replicates have been removed.
		} 

		A = boot.ci(NN,conf=0.95,type="bca")$bca[4:5] #finds 95% CI using BCa confidence intervals.
	return(list(CI=A, SampleSize = NN$R))
	}


	#Standard
	O1 = GetCIs(PearsonBS) 
	O2 = GetCIs(KendallBS)
	O3 = GetCIs(SpearmanBS)

	#RLM
	B1 = GetCIs(TbisBS)
	B2 = GetCIs(TbisBSmm)
	B3 = GetCIs(HampelBS)
	B4 = GetCIs(HuberBS)

	#lmrob
	A1 = GetCIs(n_BisquareBS)
	A2 = GetCIs(n_GGWBS)
	A3 = GetCIs(n_HampelBS)
	A4 = GetCIs(n_lqqBS)
	A5 = GetCIs(n_OptBS)
	A6 = GetCIs(n_welshBS)


	ZEROS = rbind(OZ1,OZ2,OZ3,BZ1,BZ2,BZ3,BZ4,AZ1,AZ2,AZ3,AZ4,AZ5,AZ6)
	COEFFS = rbind(OC1,OC2,OC3,BC1,BC2,BC3,BC4,AC1,AC2,AC3,AC4,AC5,AC6)
	PVALS = rbind(OP1,OP2,OP3,BP1,BP2,BP3,BP4,AP1,AP2,AP3,AP4,AP5,AP6)
	CIcombined=rbind(O1[[1]],O2[[1]],O3[[1]],B1[[1]],B2[[1]],B3[[1]],B4[[1]],A1[[1]],A2[[1]],A3[[1]],A4[[1]],A5[[1]],A6[[1]])
	SAMSIZE = rbind(O1$SampleSize,O2$SampleSize,O3$SampleSize,B1$SampleSize,B2$SampleSize,B3$SampleSize,B4$SampleSize,A1$SampleSize,A2$SampleSize,A3$SampleSize,A4$SampleSize,A5$SampleSize,A6$SampleSize)
	

	### save the data ###
	
	OUT = cbind(COEFFS, PVALS, ZEROS, SAMSIZE, CIcombined)
	colnames(OUT)=c("coefficient","pval","Zeros","bootstrap sample size","bootstrap CI 2.5%","bootstrap CI 97.5%")
	rownames(OUT)=names
	SaveOutput = OUT
	
	
	write.csv(SaveOutput, file = paste(TextInterest,".csv",sep="")) #save as csv.

}



## Functions used within the different bootstrap ensembles
FPearson <- function(data,ind){
		data=data[ind,]
		cor.test(log10(data$Cells),data$Interest, method="pearson")$estimate
	}

FKendall <- function(data,ind){
		data=data[ind,]
		cor.test(log10(data$Cells),data$Interest, method="kendall")$estimate
	}

FSpearman <-function(data,ind){
		data=data[ind,]
		cor.test(log10(data$Cells),data$Interest, method="spearman")$estimate
	}

FTbisquare <- function(data,ind){
		data=data[ind,]
		A<-rlm(Interest ~ log10(Cells), data, method="M", maxit=1000,psi="psi.bisquare")
		return( c(A$coefficients, A$converged) )
	}

FTbisquare_mm <- function(data,ind){
		data=data[ind,]
		A<-rlm(Interest ~ log10(Cells), data, method="MM", maxit=1000,psi="psi.bisquare")
		return( c(A$coefficients, A$converged) )
	}

FHuber <- function(data,ind){
		data=data[ind,]
		A<-rlm(Interest ~ log10(Cells), data, method="M", maxit=1000,psi="psi.huber")
		return( c(A$coefficients, A$converged) )
	}

FHampel <- function(data,ind){
		data=data[ind,]
		A <- rlm(Interest ~ log10(Cells), data, method="M", maxit=1000,psi="psi.hampel")
		return( c(A$coefficients, A$converged) )
	}
newBisquare <- function(data,ind){
		data=data[ind,]
		A=lmrob(Interest ~ log10(Cells), data, psi="bisquare",method="MM")
		return( c(A$coefficients, A$converged) )
	}
newGGW <- function(data,ind){
		data=data[ind,]
		A=lmrob(Interest ~ log10(Cells), data, psi="ggw",method="MM")
		return( c(A$coefficients, A$converged) )
	}

newlqq <- function(data,ind){
		data=data[ind,]
		A=lmrob(Interest ~ log10(Cells), data, psi="lqq",method="MM")
		return( c(A$coefficients, A$converged) )
	}

newHampel <- function(data,ind){
		data=data[ind,]
		A=lmrob(Interest ~ log10(Cells), data, psi="hampel",method="MM")
		return( c(A$coefficients, A$converged) )
	}
newOptimal <- function(data,ind){
		data=data[ind,]
		A=lmrob(Interest ~ log10(Cells), data, psi="optimal",method="MM")
		return( c(A$coefficients, A$converged) )
	}
newWelsh <- function(data,ind){	
		data=data[ind,]
		A=lmrob(Interest ~ log10(Cells), data, psi="welsh",method="MM")
		return( c(A$coefficients, A$converged) )
	}
