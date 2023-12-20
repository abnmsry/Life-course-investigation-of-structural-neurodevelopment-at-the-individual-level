# Life course investigation of structural neurodevelopment at the individual level

*demo.txt*: Due to data usage regulations, we are unable to upload raw data. However, researchers can obtain access to all the datasets by submitting applications. Here, we provided a simulated dataset to demo the clustering analysis. Group labels were assigned to 2,138 pariticipants randomly with corresponding probability identified in the manuscript (0.461 for Group 1,0.496 for Group 2 and 0.043 for Group 3). For different groups, we summarized the average and standard deviation for every brain regions in different age times using raw data. Then, we generated simulated regional volumes using summary-level regional statistics.

*Group_Trajectory.R*: 
> 0. Trajectory calculation
> 1. Group clustering
> 2.1. Estimation of age and region-specific GMV development among groups (Fig. 2b)
> 2.2. Estimation of group-specific developmental curve of total GMV in IMAGEN (Fig. 2c and Supplementary Fig. 4)
> 3. Estimation of peak total GMV in IMAGEN (Supplementary Fig. 5-6)

*Group_Comparison.R*:
> Comparison neurocognition measurements among groups (Fig. 2a and Supplementary Fig. 7)

*GWAS.R*:
> Group-reweighted GMV calculation (Fig. 3)

*EWAS.R*:
> 1. DMP analysis: epigenome-wide association study (Fig. 4a)
> 2. Mediation analysis (Fig. 4d/e)

total running time should be less than 1 day.

Operating system: 11th Gen Intel(R) Core(TM) i7-11800H (32.0 GB RAM); Windows 10
Software: R-4.2.2
Dependencies: dplyr-1.1.0; lme4-1.1-34; nlme-3.1-163; splines-4.2.2; mgcv-1.9.0; ggplot2-3.4.3; lavaan-0.6-16; (installtime should be less than 2 mins)

Welcome for any questions about analyses!
