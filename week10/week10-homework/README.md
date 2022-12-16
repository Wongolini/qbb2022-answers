python xpresso_predict.py pretrained_models/humanMedian_trainepoch.11-0.426.h5 input_fasta/human_promoters.fa.gz xpresso.human.out

Are the predictions different, or the same? Are there any trends that you notice looking at the data? Looking at the xpresso_predict.py script you downloaded, why do you think the predictions may differ?

Predictions differ. There are to puncta at (-1,-1) and at (1,1) but the paper predicted model and my predictions disagree at a puncta (-1,1). Inbetween (-1,1) and (1,1) is a blob of points. Each puncta has a noisy gaussian distribution of points.

Predictions could differ due to a number of reasons:
1. Pre-processing of the data differed between xpressor_predict.py and the original paper's
2. The CNN uses introduces sort of stochastic noise to the data internally. Despite being pre-trained, the model will not always make the same exact decision for a given set of data to reduce over-classifying something.
