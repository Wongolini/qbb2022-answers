Really great work! The plots you have look great. There are just a few minor issues.

1. I don't think `freqx` is giving you quite what you want. It looks like this is a lot of allele counts, but you need to do a little more work to get the allele frequencies themselves if you do it this way. I would recommend the `--freq` flag to get the AFs more easily. Also, I have your AF code, but not the plot itself (-1 point)
2. I'm not sure the `--independent-pairwise` flag is quite what we want for this regression. It looks like this does an LD-based analysis, and not necessarily the linear regression you want. Although it also looks like you might have used the `--assoc` flag, because your results look like `.qassoc`. This though I think is also not the right test as this does something called a Wald Test which is a bit different. I would try the `--linear` flag. Make sure you also include the top PCs as covariates! This messes up your manhattan plot a bit, but the manhattan plot code looks great otherwise! (no points deducted)
3. We still need your boxplot and boxplot code (-1.5 points)

Awesome job so far.

(7.5/10)
