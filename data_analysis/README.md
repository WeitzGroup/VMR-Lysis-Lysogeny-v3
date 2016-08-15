
Virome evaluation
=============================


Virome data presented in Knowles et al. 2016 figure 4a, 4b and Extended figure 4a is evaluated using R (version 3.3.1) and packages *resample* (version 0.4), *MASS* (version 7.3-45) and *robustbase* (version 0.92-6).

## Data

To perform the analysis first you will need to aquire the data aquired by Knowles et al. 2016 from [doi:10.1038/nature17193](http://dx.doi.org/10.1038/nature17193).
The data used by **Figure 4** should be placed within this folder as:
 * *nature17193-f4.csv*  

Similarly, the data  used by **Extended Figure 4a** should be placed within this folder as:
 * *nature17193-sf4.csv*  

## Analysis

 * *DataAnalysis.R*   &nbsp; - &nbsp;   script to load and prepare data for analysis using bootstrap methods
 * *CItesting.R*   &nbsp; - &nbsp;   function called by DataAnalysis.R to perform bootstrap analysis and return confidence intervals

To perform the analysis (once the data has been aquired and placed in this folder) start R in this folder and run:

```r
	source("DataAnalysis.R")
```


## Outputs
The analysis for each dataset shown below returns a table displaying the method used, Zeroes -- the number of datapoints (out of 24) assigned a weighting of zero by each method, the sample size -- the number of members of a bootstrap ensemble (out of 10,000) which successfully converged and the 95% confidence intervals calculated from these members of the bootstrap ensemble.

 * *PProphage.csv*   &nbsp; - &nbsp;   Table returning the results of the analysis for the relationship between log10(microbial abundance per ml) and percentage pro-virus like reads in viromes
 * *Integrase.csv*   &nbsp; - &nbsp;   Table returning the results of the analysis for the relationship between log10(microbial abundance per ml) and percentage integrase reads in viromes
 * *Excisionase.csv*   &nbsp; - &nbsp;   Table returning the results of the analysis for the relationship between log10(microbial abundance per ml) and percentage excisionase reads in viromes


###References


Knowles B. et al. 2016. Lytic to temperate switching of viral communities. *Nature*. 531: 466–470. [doi:10.1038/nature17193](http://dx.doi.org/10.1038/nature17193)

R Core Team. 2016. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria.  URL https://www.R-project.org/.

Hesterberg T. 2015. resample: Resampling Functions. R package version 0.4. https://CRAN.R-project.org/package=resample

Venables W. N., Ripley B. D. 2002. Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0

Rousseeuw P. et al. 2016. robustbase: Basic Robust Statistics. R  package version 0.92-6. URL  http://CRAN.R-project.org/package=robustbase



