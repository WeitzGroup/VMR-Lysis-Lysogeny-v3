# Part of the code used in:
# Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
# 
# From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny-v2
# MIT License


CItesting <- function(DATA,SEED,SAMPLES,TextInterest) {

### Create bootstrap ensembles for each metric that is considered ###

## Standard approaches

#1) Pearson
	print("...performing Pearson bootstrap...")
	
	set.seed(SEED) # for reproductibility
	PearsonBS = bootstrap(DATA, FPearson, SAMPLES) #perform a 10,000 sample bootstrap


#2) Kendall
	print("...performing Kendall bootstrap...")

	set.seed(SEED) # for reproductibility
	KendallBS = bootstrap(DATA, FKendall, SAMPLES) #perform a 10,000 sample bootstrap
		
#3) Spearman
	print("...performing Spearman bootstrap...")

	set.seed(SEED) # for reproductibility
	SpearmanBS = bootstrap(DATA, FSpearman, SAMPLES) #perform a 10,000 sample bootstrap


## Robust regression techniques using the MASS package

#1) Tukey bisquare
	print("...performing Bisquare M (MASS) bootstrap...")
	
	set.seed(SEED) # for reproductibility
	TbisBS = bootstrap(DATA, FTbisquare, SAMPLES) #perform a 10,000 sample bootstrap
	

#2) Tukey bisquare -MM
	print("...performing Bisquare MM (MASS) bootstrap...")

	set.seed(SEED) # for reproductibility
	TbisBSmm = bootstrap(DATA, FTbisquare_mm, SAMPLES) #perform a 10,000 sample bootstrap


#3) Huber
	print("...performing Huber M (MASS) bootstrap...")

	set.seed(SEED) # for reproductibility
	HuberBS = bootstrap(DATA, FHuber, SAMPLES) #perform a 10,000 sample bootstrap

#4) Hampel
	print("...performing Hampel M (MASS) bootstrap...")

	set.seed(SEED) # for reproductibility
	HampelBS = bootstrap(DATA, FHampel, SAMPLES) #perform a 10,000 sample bootstrap
	

## Robust regression techniques using the robustbase package

#1) Bisquare	
	print("...performing Bisquare MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility
	n_BisquareBS = bootstrap(DATA, newBisquare, SAMPLES)


#2) GGW
	print("...performing GWW MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility	
	n_GGWBS = bootstrap(DATA, newGGW, SAMPLES)


#3) LQQ
	print("...performing LQQ MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility
	n_lqqBS = bootstrap(DATA, newlqq, SAMPLES)


#4) Hampel
	print("...performing Hampel MM (robustbase) bootstrap...")
	
	set.seed(SEED) # for reproductibility
	n_HampelBS = bootstrap(DATA, newHampel, SAMPLES)


#5) Optimal
	print("...performing Optimal MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility
	n_OptBS = bootstrap(DATA, newOptimal, SAMPLES)


#6) Welsh
	print("...performing Welsh MM (robustbase) bootstrap...")

	set.seed(SEED) # for reproductibility
	n_welshBS = bootstrap(DATA, newWelsh, SAMPLES)	

		

### Analyse the bootstraps for each metric ###


	names = c("Pearson","Kendall","Spearman","M Bisquare","MM Bisquare (rlm)","M Hampel","M Huber","MM Bisquare(lmrob)","MM GGW","MM LQQ","MM Hampel","MM Optimal","MM Welsh")

## Find the 95% confidence intervals ##
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

## Find the number of datapoints which are assigned a weighting of zero by each method	
	A = rlm(Interest ~ log10(Cells), DATA, psi="psi.bisquare",method="M", maxit=1000)
	OZ1 = sum(A$w==0)
	
	A = rlm(Interest ~ log10(Cells), DATA, psi="psi.bisquare",method="MM", maxit=1000)
	OZ2 = sum(A$w==0)
	
	A = rlm(Interest ~ log10(Cells), DATA, psi="psi.hampel",method="M", maxit=1000)
	OZ3 = sum(A$w==0)
	
	A = rlm(Interest ~ log10(Cells), DATA, psi="psi.huber",method="M", maxit=1000)
	OZ4 = sum(A$w==0)
	

	A = lmrob(Interest ~ log10(Cells), DATA, psi="bisquare",method="MM")
	Z1 = sum(A$rweights==0)
	
	A = lmrob(Interest ~ log10(Cells), DATA, psi="GGW",method="MM")
	Z2 = sum(A$rweights==0)
	
	A = lmrob(Interest ~ log10(Cells), DATA, psi="LQQ",method="MM")
	Z3 = sum(A$rweights==0)
	
	A = lmrob(Interest ~ log10(Cells), DATA, psi="Hampel",method="MM")
	Z4 = sum(A$rweights==0)
	
	A = lmrob(Interest ~ log10(Cells), DATA, psi="Optimal",method="MM")
	Z5 = sum(A$rweights==0)
	
	A = lmrob(Interest ~ log10(Cells), DATA, psi="Welsh",method="MM")
	Z6 = sum(A$rweights==0)
	

	ZEROS = rbind(0,0,0,OZ1,OZ2,OZ3,OZ4,Z1,Z2,Z3,Z4,Z5,Z6)
	CIcombined=rbind(O1[[1]][1,],O2[[1]][1,],O3[[1]][1,],B1[[1]][2,],B2[[1]][2,],B3[[1]][2,],B4[[1]][2,],A1[[1]][2,],A2[[1]][2,],A3[[1]][2,],A4[[1]][2,],A5[[1]][2,],A6[[1]][2,])
	SAMSIZE = rbind(O1$SampleSize,O2$SampleSize,O3$SampleSize,B1$SampleSize,B2$SampleSize,B3$SampleSize,B4$SampleSize,A1$SampleSize,A2$SampleSize,A3$SampleSize,A4$SampleSize,A5$SampleSize,A6$SampleSize)

	### save the data ###
	
	OUT = cbind(ZEROS,SAMSIZE,CIcombined)
	colnames(OUT)=c("Zeros","sample size","CI 2.5%","CI 97.5%")
	rownames(OUT)=names
	SaveOutput = OUT
	
	
	write.csv(SaveOutput, file = paste(TextInterest,".csv",sep="")) #save as csv.

}

#Function to find the CIs of a bootstrap only using the members when the measure is valid
GetCIs <- function(x) {

	#x is a bootstrap object

	#1. Find the replicates within the bootstrap which converged
	if(dim(x$replicates)[2]==3) {
	USE <- which(x$replicates[,3]==1)  #row indices of converged regressions
	L_USE <- length(USE) # number converged

	#2. Make new bootstrap object containing only converged values
	x2 <- x
	x2$replicates <- x$replicates[USE,]
	x2$R <- L_USE
	}
	else {
	x2 <- x
	L_USE <- x2$R
	}


	#3. CI's
	A <- CI.percentile(x2,expand=FALSE,probs=c(0.025,0.975)) 

return(list(Exact=A, SampleSize = L_USE))
}

## Functions used within the different bootstrap ensembles
FPearson <- function(data){
		cor.test(log10(data$Cells),data$Interest, method="pearson")$estimate
	}

FKendall <- function(data){
		cor.test(log10(data$Cells),data$Interest, method="kendall")$estimate
	}

FSpearman <-function(data){
		cor.test(log10(data$Cells),data$Interest, method="spearman")$estimate
	}

FTbisquare <- function(data){
		A<-rlm(Interest ~ log10(Cells), data, method="M", maxit=1000,psi="psi.bisquare")
		return( c(A$coefficients, A$converged) )
	}

FTbisquare_mm <- function(data){
		A<-rlm(Interest ~ log10(Cells), data, method="MM", maxit=1000,psi="psi.bisquare")
		return( c(A$coefficients, A$converged) )
	}

FHuber <- function(data){
		A<-rlm(Interest ~ log10(Cells), data, method="M", maxit=1000,psi="psi.huber")
		return( c(A$coefficients, A$converged) )
	}

FHampel <- function(data){
		A <- rlm(Interest ~ log10(Cells), data, method="M", maxit=1000,psi="psi.hampel")
		return( c(A$coefficients, A$converged) )
	}
newBisquare <- function(data){
		A=lmrob(Interest ~ log10(Cells), data, psi="bisquare",method="MM")
		return( c(A$coefficients, A$converged) )
	}
newGGW <- function(data){
		A=lmrob(Interest ~ log10(Cells), data, psi="ggw",method="MM")
		return( c(A$coefficients, A$converged) )
	}

newlqq <- function(data){
		A=lmrob(Interest ~ log10(Cells), data, psi="lqq",method="MM")
		return( c(A$coefficients, A$converged) )
	}

newHampel <- function(data){
		A=lmrob(Interest ~ log10(Cells), data, psi="hampel",method="MM")
		return( c(A$coefficients, A$converged) )
	}
newOptimal <- function(data){
		A=lmrob(Interest ~ log10(Cells), data, psi="optimal",method="MM")
		return( c(A$coefficients, A$converged) )
	}
newWelsh <- function(data){
		A=lmrob(Interest ~ log10(Cells), data, psi="welsh",method="MM")
		return( c(A$coefficients, A$converged) )
	}
