#### 
## 1. load_data.R loads the data and makes the data frame oxy and the
###     reduced data frame oxy_reduced.  The reduced data frame
###     contains the most abundant microbes as well as a known oxalate 
##      degrader.

### 2.  some_data_plots.R has some code to make a few plots of the data
##      just to get a handle on it.

##  3. fit.R estimates the interaction parameters, growth rates, etc 
##      using regularized least squares and optimx.  There are way better
##      ways to do this in the future..  This takes forever (maybe 15 hours)
##		to run.

##  4.  plot_fits.R has some code to make some plots of the fits.

##  5.  heatmap_beta_matrix.R makes some heatmaps of the interaction matrix.

#an edit made in ben_working branch
