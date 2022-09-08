# day 4 homework feedback

I really like how you've labeled the heatmap squares with the actual observed powers! Did you set the `vmin` argument for both heatmaps?

(Completely optional feedback: Do you like having the values displayed such that the y-axis goes from the top down, increasing number of coin tosses as you move down the axis? and relatedly, such that the x-axis goes from the left to the right, decreasing the probability of heads as you move across the axis? If you had a 2D matrix where the rows corresponded to the probability of heads (first axis) and the columns corresponded to the number of tosses (second axis) in your numpy 2d storage array, will the axes be displayed differently and more intuitively?)

great code overall! The assignment does ask you to incorporate the nested for loop and storing the power within the `run_experiment()` function rather than calling that function repeatedly. What you've done works and leads to the same conclusions, but consider how you might edit the function to incorporate the new code. While your results will be the same every time your script is run, will your results be the same compared to someone who only calls the `run_experiment()` function once therefore only sets the random seed once?


You've got great thoughts within the README, but you don't compare the homework simulation to the S13 figure like requested. The S13 figure doesn't use the Hidden Markov Model at all. Happy to discuss rhapsodi in more detail if you're interested, but for this assignment, please compare to figure S13.
