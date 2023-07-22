# timeseries-analysis
* For this project 2 different datasets have been used, the first one is the weather analysis (https://www.kaggle.com/datasets/mastmustu/weather-analysis?resource=download) and the second one is sunspots (https://www.kaggle.com/datasets/tavoglc/sunspots). Firstly choose the years of interest (for this project 2014 and 2019) and then we remove the redundant features from the first dataset. We retime the timetable to daily values for both datasets and then we merge them. We plot the values of every feature on the final dataset. Then we explore the autocorrelation of the features and apply white noise process to remove it.
* We create and plot correlation network with weighted edges (that are important according to their p-values) connecting the nodes, then we plot directed weighted graphs to visualize Granger causality and conditional Granger causality values between features.
* We then use regression (without feature selection, then with feature selection), training a model with the data from one year and testing the model on the data from the other.

fitAR.m, GCI.m, CGCI.m are the white noise process, Granger causality and conditional Granger causality implemetations respectively used in the main_project.m.