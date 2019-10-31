# EALOPS
Analysis codes for Zhao et al. (2019) Pupillometry as an Objective Measure of Sustained Attention in Young and Older Listeners. Trends in Hearing.

Before running the scripts, check the line 2 of all analysis scripts. Please make sure to change the variable 'path_in' to the directory of the processed data.

1. Plot Behaviour: run 'plot_main_behaviour.m' to plot behavioural bar charts (hit, false alarm and number of bad trials). 
* Uncomment line 3 ('Expt1 - Young') to plot Figure 2
* Uncomment line 4 ('Expt2 - Old') to plot Figure 6
2. Plot PDR: run 'plot_main_PDR_behaviourTimeBinned_corr.m' to plot PDR results and time-binned behavioural results
* Uncomment line 3 ('Expt1 - Young') to plot Figure 3
* Uncomment line 4 ('Expt2 - Old') to plot Figure 7
3. Plot pupilmetrics' pupil data: run 'plot_pupilmetrics.m' to plot Figure 4.
