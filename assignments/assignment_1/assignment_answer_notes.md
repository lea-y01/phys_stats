notes for assignment answers: 

* Which pair of dimensions most strongly separates high- and low-Prob stars? **Choose exactly one: RA–Dec, Parallax–pmRA, pmRA–pmDE**
* Which dimension contributes least to separation? **Choose exactly one: RA, Dec, Parallax, pmRA, pmDE**
* Which statement matches best the evidence in your plots?
    1. The algorithm primarily exploits spatial clustering (i.e. in position)
    2. The algorithm primarily exploits kinematic coherence (i.e. in stellar motion)
    3. The separation is weak and likely noise-dominated.

### Cluster NGC_6402: 
i. I would say that parallax-pmDE most strongly separates high- and low-prob stars. I know that the effect is also visible in Parallax-RAdeg and Parallax-DEdeg, but I found it more interesting that the distribution of points with my choice (Plx-PM) appears more asymmetric and spread out (there appears some spatial uniformity in the Plx-RAdeg&DEdeg plots).  
ii. either RAdeg or DEdeg contributes the least to separation
iii. The algorithm primarily exploits kinematic coherence (i.e. stellar motion), as the separation and spreading is more visible in the proper motion plots.

### Cluster NGC_1039: 
i. Parallax-DEdeg seems to separate the two groups the most. 
ii. pmRA seems to contribute the least to separation, which I notice the most in hte pmDE-pmRA plot, where the scatter distribution looks quite spread out and shared. 
iii. The algorithm primarily exploits spatial clustering. 

### Cluster NGC_362: 
i. Parallax-pmDE seems to most strongly separate the groups; there is a noticeable tilt in the scatter plot. I notice that the relative difference in max probability density between the groups for DEdeg's histogram is the largest (~9/2), but the spatial PDFs are so closely overlapped compared to the other three. When looking to see which of the proper motion seemed to be more influential, I compare the proper motion vs. RA/DEC plots (bottom left four) to try to distinguish influence. pmRA-RAdeg & pmRA-DEdeg seem almost reflections of each other, while in the two pmDE plots there is clear contrast in symmetry.  
ii.  RAdeg seem to contribute the least. 
iii. Algorithm primarily exploits kinematic coherence (based on previous answers and big influence of pmDE).

### Trumpler_20: 
i. The separation here feels harder to distinguish. But again, Parallax-pmDE seems to separate the groups the most. However, I also notice that between RAdeg and DEdeg, there is a clear spreading with the declination dependence.
ii. RAdeg seems to contribute the least, as its relevant scatter plots seem consistent regardless of right ascension. 
iii. Both spatial and kinematic aspects seem to be at play here, so I will say the separation is too weak and likely noise-dominated. Especially also because the PDFs between the groups overlap so much compared to other plots. 

* Self-check: Briefly, answer the following two questions: Could this explanation have been written without looking at any Gaia data? If yes, go back and make sure your answers are backed up with observations from your figures. 

I don't think these explanations could be made without looking at Gaia data, since there is some nuance in my answers.

* Which observation changes the most between different clusters? 

This feels also hard to answer definitively, as the plots seem to all change a lot between different clusters. One thing I did notice changing was the range of possible values (the x-axis ) in the histograms, specifically with the proper motion plots. At first I though this might just be a question of sample size (because the cluster with the largest number of candidates, NGC_362, had noticeably wider distriutions especially in the proper motion parameters). But then comparing the second highest count cluster Trumpler_20 (count = 1692) and the lowest, NGC_1039 (1095), I cannot infer a relationship between candidate count and the possibility range. For example: for NGC_1039 (c = 1095) RAdeg hist is between 36-44. NGC_6402 (c = 1335)'s range is about 264.2 - 264.8 for the same plot, while NGC_632 (c = 8941)'s is 10-25. This difference is observable also in the proper motion histograms as well, and also does not appear to correspond to sample size (pmDE histogram range for lowest count ~ 3, highest count ~ 5, second highest ~ 0.5) 

I now realize this might have to do with the observations, which still does not feel completely intuitive given that these are absolute measurements being derived from observations made with the same instrument. 

