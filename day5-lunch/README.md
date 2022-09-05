I did the entire assignment in python despite instructions because Rajiv said if we felt more comfortable in Python we could take that liberty.

3. 
Yes, from the sample space, there appears to be a weak positive linear relationship shared between the mother's age and number of mutations with a slope of +15.43. However, this is assuming the father's age shares a linear relationship with n_muts; the R^2 is relatively weak.
However, the pval is much below .05 so the relationship is strongly positively significant.

Intercept     15.436470 (n_muts/Mother_age)
Mother_age     1.885151
R-squared:     0.547

Pval:  Intercept     5.012000e-10
Mother_age    1.152020e-69

4.
Yes, from the sample space, there appears to be a positie linear relationship shared between the father's age and number of mutations with a slope of +16.85. However, the pval is much below .05 so the relationship is strongly positively significant.
Intercept     16.845773 (n_muts/Father_age)
Father_age     1.620475
R-squared:     0.607

Pval:  Intercept     8.653418e-15
Father_age    5.542092e-82

6.
T-TEST
Ttest_indResult(statistic=-53.40356528726923, pvalue=2.1986031


7.
OLS MODEL PREDICTS 
 Father_age=50.5 --> n_muts=98.67976753830877