Solid work. I think the way you're doing traceback is probably a little bit redundant. It looks like you're essentially re-doing the work you did to fill out the F-matrix, as you move backwards. And then you're going forwards again, through the path that you just calculated. It should be faster to just keep track of traceback as you move forwards (although, to be fair, more memory intensive), and then you just have to go backwards through the traceback a single time.

Your alignment looks awesome though. Nice work, really. Unfortunately, I have to point out that you might be missing the last bit of each alignment (the remainder, after dividing by 50).

The only thing we're missing is the number of gaps in each sequence (-0.5)

9.5/10
