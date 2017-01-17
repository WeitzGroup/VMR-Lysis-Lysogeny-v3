

Model evaluation
=============================


Three models are considered in the paper. Latin hypercube sampling was performed to find 10,000 sets of parameters (within reasonable biological ranges) for each model. The corresponding steady state densities found using each of these parameter sets for microbial hosts and their viruses are then plotted. These are shown in Figure S1 of Weitz et al. (In submission). To recreate each of the subplots showing the relationship between the virus-to-microbe ratio (VMR) and microbial density at steady state for a particular model run each of the corresponding MATLAB files:

<br><br>

LV — Lotka-Volterra type model

 * *fig_LV_ratio.m*

PtW — Piggyback-the-Winner model  (Knowles et al. 2016)

 * *fig_PtW_ratio.m*

WD — Weitz & Dushoff model  (Weitz & Dushoff, 2008)

 * *fig_WD_ratio.m*
<br><br>

*lv_stats.mat*, *ptw_stats.mat* and *wd_stats.mat* contain the data used, for each respective model, in the manuscript.

Additional functions: 
 * *LHSmid.m*  &nbsp; - &nbsp;    function to perform Latin hypercube sampling
 * *array2vstruct.m*  &nbsp;  - &nbsp;   function to convert matrix array object into struct object     
 * *hist2d_mat.m*   &nbsp; -  &nbsp;   function used for creating contours using data density
 * *psprintc.m* &nbsp; -  &nbsp;   function used when saving figures

 
<br><br>

#### References

Weitz J.S., Beckett S.J., Brum J.R., Cael B.B., Dushoff J. Lysis, Lysogeny, and Virus-Microbe Ratios. *In submission.*

Knowles B. et al. 2016. Lytic to temperate switching of viral communities. *Nature*. 531: 466–470. [doi:10.1038/nature17193](http://dx.doi.org/10.1038/nature17193)

Weitz J.S., Dushoff J. 2008. Alternative stable states in host-phage dynamics. *Theoretical Ecology*. 1: 13-19. [doi:10.1007/s12080-007-0001-1](http://dx.doi.org/10.1007/s12080-007-0001-1)
