# Part of the code used in:
# Weitz et al. Lysis, Lysogeny, and Virus-Microbe Ratios
# 
# From https://github.com/WeitzGroup/VMR-Lysis-Lysogeny-v3
# MIT License



##  R file for plotting the lysogeny hallmark gene relationships with log10(Microbial density) from Knowles et al. 2016.

### Import data ###


# This is downloaded from Knowles et al. 2016 (see "Figure 4: Temperate features in viromes increase with host density." from http://dx.doi.org/10.1038/nature17193)
	
	data4 <- read.csv("nature17193-f4.csv", header=TRUE)


# This is downloaded from the SI in Knowles et al. 2016 (see "Extended Data Figure 4: Temperateness of viral communities increases with host density and viral functional composition change." from http://dx.doi.org/10.1038/nature17193)

	dataSI4a <- read.csv("nature17193-sf4.csv", header=TRUE)
	
### Prepare data for analysis ###


	DATA = as.data.frame(cbind(data4$Microbes.per.ml,data4$X..integrase.reads,data4$X..excisionase.reads,data4$Functional.Diversity..H..,data4$X..pathogenicity.reads,dataSI4a$X..prophage.reads))
	COLNAMES = c("Cells","Integrase","Excisionase","FuncDiv","Pathogenicity","PProphage")
	colnames(DATA) = COLNAMES


#### PERFORM REGRESSIONS and find p-values

#provirus-like reads vs. log10(Cells)
PP_MODEL=lm(PProphage ~ log10(Cells), data=DATA)
PP_SUM = summary(PP_MODEL)
PP_Rsq = PP_SUM$r.squared
PP_pval = pf(PP_SUM$fstatistic[1],PP_SUM$fstatistic[2],PP_SUM$fstatistic[3],lower.tail=FALSE) #gives same p-value as shown in summary(PP_MODEL)
PP_FIT = coef(PP_MODEL)[1] + coef(PP_MODEL)[2]*log10(DATA$Cells)

PP_cor = cor.test(log10(DATA$Cells),DATA$PProphage,method="pearson")
PP_cor_pval = PP_cor$p.value  #equivalent to PP_pval
PP_cor_r2 = PP_cor$estimate^2 # equivalent to PP_Rsq

#integrase reads vs. log10(Cells)
INT_MODEL=lm(Integrase ~ log10(Cells), data=DATA)
INT_SUM = summary(INT_MODEL)
INT_Rsq = INT_SUM$r.squared
INT_pval = pf(INT_SUM$fstatistic[1],INT_SUM$fstatistic[2],INT_SUM$fstatistic[3],lower.tail=FALSE)  #gives same p-value as shown in summary(PP_MODEL)
INT_FIT = coef(INT_MODEL)[1] + coef(INT_MODEL)[2]*log10(DATA$Cells)

INT_cor = cor.test(log10(DATA$Cells),DATA$Integrase,method="pearson")
INT_cor_pval = INT_cor$p.value  #equivalent to INT_pval
INT_cor_r2 = INT_cor$estimate^2 # equivalent to INT_Rsq

#excisionase reads vs. log10(Cells)
EXC_MODEL=lm(Excisionase ~ log10(Cells), data=DATA)
EXC_SUM = summary(EXC_MODEL)
EXC_Rsq = EXC_SUM$r.squared
EXC_pval = pf(EXC_SUM$fstatistic[1],EXC_SUM$fstatistic[2],EXC_SUM$fstatistic[3],lower.tail=FALSE) #gives same p-value as shown in summary(PP_MODEL)
EXC_FIT = coef(EXC_MODEL)[1] + coef(EXC_MODEL)[2]*log10(DATA$Cells)

EXC_cor = cor.test(log10(DATA$Cells),DATA$Excisionase,method="pearson")
EXC_cor_pval = EXC_cor$p.value  #equivalent to EXC_pval
EXC_cor_r2 = EXC_cor$estimate^2 # equivalent to EXC_Rsq


#### PLOT THE DATA
COL="steelblue4" #line colour
ROUND=3 # number of significant digits
PCH=19 # point shape
LWD=2 #line width
#plot height to write text values
 
h_rsq = 0.95
h_p = 0.85

dev.new()
par(mfrow=c(3,1))
par(mar=c(5,4,1,1))
XX = 5.4 #text x-position in fig

#Panel 1: provirus-like reads
MAX=8.5 #Ylim

plot(log10(DATA$Cells), DATA$PProphage, ylim=c(0,MAX), ylab="% provirus-like reads in virome", xlab = expression(paste("log10(Microbial cells ml" ^-1, ")")),pch=PCH)
#lines(log10(DATA$Cells), PP_FIT, col=COL,lwd=LWD)
text(XX,h_rsq*MAX,bquote("R"^2 == .(round(PP_Rsq,ROUND))),pos=4)
text(XX,h_p*MAX,bquote(paste(italic(p)," = ", .(round(PP_pval,ROUND)))),pos=4)


#Panel 2: Integrase reads
MAX=3.5 #Ylim

plot(log10(DATA$Cells), DATA$Integrase, ylim=c(0,MAX), ylab="% integrase reads in virome", xlab = expression(paste("log10(Microbial cells ml" ^-1, ")")),pch=PCH)
#lines(log10(DATA$Cells), INT_FIT, col=COL,lwd=LWD)
text(XX,h_rsq*MAX,bquote("R"^2 == .(round(INT_Rsq,ROUND))),pos=4)
text(XX,h_p*MAX,bquote(paste(italic(p)," = ", .(round(INT_pval,ROUND)))),pos=4)


#Panel 3: Excisionase reads
MAX=0.2 #Ylim

plot(log10(DATA$Cells), DATA$Excisionase, ylim=c(0,MAX), ylab="% excisionase reads in virome", xlab = expression(paste("log10(Microbial cells ml" ^-1, ")")),pch=PCH)
#lines(log10(DATA$Cells), EXC_FIT, col=COL,lwd=LWD)
text(XX,h_rsq*MAX,bquote("R"^2 == .(round(EXC_Rsq,ROUND))),pos=4)
text(XX,h_p*MAX,bquote(paste(italic(p)," = ", .(round(EXC_pval,ROUND)))),pos=4)


dev.copy2pdf(file="ViromePanel.pdf",width=89*(1/25.4),height=190*(1/25.4)) #converted mm to inches
dev.copy2eps(file="ViromePanel.eps",,width=89*(1/25.4),height=190*(1/25.4))
dev.copy(jpeg, file="ViromePanel.jpg",units="mm",width=89,height=190, res=300)
dev.off()









