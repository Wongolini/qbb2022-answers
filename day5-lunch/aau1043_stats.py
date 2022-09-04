#%%
import numpy as np 
import scipy as sp 
import matplotlib.pyplot as plt
from day5_lunch import get_header
'''
Use numpy genfromtxt to load the “joined” data from step 3 into a numpy array. 
Use the Names = option to give your fields informative names.

Use matplotlib to plot:
the count of maternal de novo mutations vs. maternal age (upload as ex2_a.png)
the count of paternal de novo mutations vs. paternal age (upload as ex2_b.png)
    >Use ordinary least squares smf.ols() to test for an association between maternal age 
    and maternally inherited de novo mutations.

    >Is this relationship significant?
    >What is the size of this relationship?
    >Use ordinary least squares smf.ols() to test for an association between paternal age and 
     paternally inherited de novo mutations.
    >Is this relationship significant?
    >What is the size of this relationship?
    >Plot a histogram of the number of maternal de novo mutations and paternal de novo mutations 
     per proband on a single plot with semi-transparency (and upload as ex2_c.png).

    > Test whether the number of maternally inherited de novo mutations per proband is significantly 
     different than the number of paternally inherited de novo mutations per proband.

    > Predict the number of paternal de novo mutations for a proband with a father who was 50.5 years 
     old at the proband’s time of birth.

'''
file='aau1043_merged.tsv'
data = np.genfromtxt(file,delimiter='\t',dtype=None,encoding='utf-8',names=get_header(file))
# %%
