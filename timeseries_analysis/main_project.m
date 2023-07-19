% setting parameters
p=2;
alpha=0.05;
zalpha=1.96;
alpha2=0.08;

% inserting the dataset
sunspots=readtable("sunspots.csv");
weather=readtable("climate_data.csv");
weather=removevars(weather,{'Date1'}); % removing the redundant column (Date1=Date)
weather=removevars(weather,{'Month'}); % removing the redundant column (month is included in the date)
weather=removevars(weather,{'MaximumRainPerMinute'}); % removing the redundant column (value = 0 for every day)
weather.Date.Format='uuuu-MM-dd';% setting the format so the data of the date are accessible
sunspots.date.Format='uuuu-MM-dd';

% choosing the years of interest and filling missing days
weather_2014=weather(year(weather.Date)==2014,:);
sunspots_2014=sunspots(year(sunspots.date)==2014,:);
DailySunspots=array2table(sunspots_2014.dailysunspots);
weather_2014=table2timetable(weather_2014);
weather_2014=retime(weather_2014,"daily","spline");
weather_2014=[weather_2014,DailySunspots];
weather_2014=renamevars(weather_2014,["Var1"],["DailySunspots"]);
weather_2019=weather(year(weather.Date)==2019,:);
sunspots_2019=sunspots(year(sunspots.date)==2019,:);
DailySunspots=array2table(sunspots_2019.dailysunspots);
weather_2019=table2timetable(weather_2019);
weather_2019=retime(weather_2019,"daily","spline");
weather_2019=[weather_2019,DailySunspots];
weather_2019=renamevars(weather_2019,["Var1"],["DailySunspots"]);
clear weather sunspots_2014 sunspots_2019 DailySunspots sunspots
n14=length(weather_2014.Date);
n19=length(weather_2019.Date);

% Plotting Average direction (deg)
fig_1=figure(1);
t = tiledlayout(2,1,'TileSpacing','Compact');
subplot(2,1,1)
plot(weather_2014.Date,weather_2014.AverageDirection__deg_,'Color',[0.5,0.1,0.5])
ylim([0 360])
ylabel('Wind direction (deg)');
title('Average direction of the wind for 2014');
subplot(2,1,2)
plot(weather_2019.Date,weather_2019.AverageDirection__deg_,'Color',[0.5,0.1,0.5])
ylim([0 360])
ylabel('Wind direction (deg)');
title('Average direction of the wind for 2019');

% Plotting Temperature related variables
fig_2=figure(2);
t = tiledlayout(2,1,'TileSpacing','Compact');
subplot(2,1,1)
plot(weather_2014.Date,weather_2014.AverageTemperature__F_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2014.Date,weather_2014.MaximumTemperature__F_,'Color',[0.7,0.2,0.2])
plot(weather_2014.Date,weather_2014.MinimumTemperature__F_,'Color',[0.2,0.2,0.7])
legend('Average Temperature','Maximum Temperature','Minimum Temperature')
ylabel('Temperature (F)');
title('Average, Maximum and Minimum Temperature for 2014');
ylim([-40 100])
hold off
subplot(2,1,2)
plot(weather_2019.Date,weather_2019.AverageTemperature__F_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2019.Date,weather_2019.MaximumTemperature__F_,'Color',[0.7,0.2,0.2])
plot(weather_2019.Date,weather_2019.MinimumTemperature__F_,'Color',[0.2,0.2,0.7])
legend('Average Temperature','Maximum Temperature','Minimum Temperature')
ylabel('Temperature (F)');
title('Average, Maximum and Minimum Temperature for 2019');
ylim([-40 100])
hold off

% Plotting Average Dwe Point and Maximum heat index variables
fig_3=figure(3);
t = tiledlayout(2,1,'TileSpacing','Compact');
subplot(2,1,1)
plot(weather_2014.Date,weather_2014.AverageDewpoint__F_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2014.Date,weather_2014.MaximumHeatIndex__F_,'Color',[0.7,0.2,0.2])
legend('Average Dewpoint','Maximum Heat Index')
ylabel('Temperature (F)');
title('Average Dewpoint and Maximum Heat Index for 2014');
hold off
subplot(2,1,2)
plot(weather_2019.Date,weather_2019.AverageDewpoint__F_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2019.Date,weather_2019.MaximumHeatIndex__F_,'Color',[0.7,0.2,0.2])
legend('Average Dewpoint','Maximum Heat Index')
ylabel('Temperature (F)');
title('Average Dewpoint and Maximum Heat Index for 2019');
hold off

% Plotting Humidity related variables
fig_4=figure(4);
t = tiledlayout(2,1,'TileSpacing','Compact');
subplot(2,1,1)
plot(weather_2014.Date,weather_2014.AverageHumidity___,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2014.Date,weather_2014.MaximumHumidity___,'Color',[0.7,0.2,0.2])
plot(weather_2014.Date,weather_2014.MinimumHumidity___,'Color',[0.2,0.2,0.7])
legend('Average Humidity','Maximum Humidity','Minimum Humidity')
ylabel('% Humidity');
title('Average, Maximum and Minimum Humidity for 2014');
ylim([0 100])
hold off
subplot(2,1,2)
plot(weather_2019.Date,weather_2019.AverageHumidity___,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2019.Date,weather_2019.MaximumHumidity___,'Color',[0.7,0.2,0.2])
plot(weather_2019.Date,weather_2019.MinimumHumidity___,'Color',[0.2,0.2,0.7])
legend('Average Humidity','Maximum Humidity','Minimum Humidity')
ylabel('% Humidity');
title('Average, Maximum and Minimum Humidity for 2019');
ylim([0 100])
hold off

% Plotting Pressure related variables
fig_5=figure(5);
t = tiledlayout(4,1,'TileSpacing','Compact');
subplot(4,1,1)
plot(weather_2014.Date,weather_2014.AverageBarometer_in_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2014.Date,weather_2014.MaximumPressure,'Color',[0.7,0.2,0.2])
plot(weather_2014.Date,weather_2014.MinimumPressure,'Color',[0.2,0.2,0.7])
legend('Average Barometer','Maximum Pressure','Minimum Pressure')
ylabel('Pressure (in)');
title('Average, Maximum and Minimum Pressure for 2014');
ylim([27 31.5])
hold off
subplot(4,1,2)
plot(weather_2014.Date,weather_2014.diff_pressure,'Color',[0.5,0.1,0.5])
title('Pressure Difference for 2014')
ylabel('Pressure (in)');
ylim([0 2.5])
subplot(4,1,3)
plot(weather_2019.Date,weather_2019.AverageBarometer_in_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2019.Date,weather_2019.MaximumPressure,'Color',[0.7,0.2,0.2])
plot(weather_2019.Date,weather_2019.MinimumPressure,'Color',[0.2,0.2,0.7])
legend('Average Barometer','Maximum Pressure','Minimum Pressure')
ylabel('Pressure (in)');
title('Average, Maximum and Minimum Pressure for 2019');
ylim([27 31.5])
hold off
subplot(4,1,4)
plot(weather_2019.Date,weather_2019.diff_pressure,'Color',[0.5,0.1,0.5])
title('Pressure Difference for 2019')
ylabel('Pressure (in)');
ylim([0 2.5])

% Plotting Wind and Gust related variables
fig_6=figure(6);
t = tiledlayout(4,1,'TileSpacing','Compact');
subplot(4,1,1)
plot(weather_2014.Date,weather_2014.AverageWindspeed_mph_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2014.Date,weather_2014.MaximumWindspeed_mph_,'Color',[0.7,0.2,0.2])
legend('Average Windspeed','Maximum Windspeed')
ylabel('Windspeed (mph)');
title('Average and Maximum Windspeed for 2014');
ylim([0 60])
hold off
subplot(4,1,2)
plot(weather_2014.Date,weather_2014.AverageGustspeed_mph_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2014.Date,weather_2014.MaximumGustSpeed_mph_,'Color',[0.7,0.2,0.2])
legend('Average Gustspeed','Maximum Gustspeed')
ylabel('Gustspeed (mph)');
title('Average and Maximum Gustspeed for 2014');
ylim([0 60])
hold off
subplot(4,1,3)
plot(weather_2019.Date,weather_2019.AverageWindspeed_mph_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2019.Date,weather_2019.MaximumWindspeed_mph_,'Color',[0.7,0.2,0.2])
legend('Average Windspeed','Maximum Windspeed')
ylabel('Windspeed (mph)');
title('Average and Maximum Windspeed for 2019');
ylim([0 60])
subplot(4,1,4)
plot(weather_2019.Date,weather_2019.AverageGustspeed_mph_,'Color',[0.2,0.7,0.2])
hold on
plot(weather_2019.Date,weather_2019.MaximumGustSpeed_mph_,'Color',[0.7,0.2,0.2])
legend('Average Gustspeed','Maximum Gustspeed')
ylabel('Gustspeed (mph)');
title('Average and Maximum Gustspeed for 2019');
ylim([0 60])
hold off

% Plotting Rainfall related variables
fig_7=figure(7);
t = tiledlayout(4,1,'TileSpacing','Compact');
subplot(4,1,1)
plot(weather_2014.Date,weather_2014.RainfallForMonth_in_,'Color',[0.5,0.1,0.5])
ylabel('Rainfall (in)');
title('Rainfall per month for 2014');
ylim([0 4])
subplot(4,1,2)
plot(weather_2014.Date,weather_2014.RainfallForYear_in_,'Color',[0.5,0.1,0.5])
ylabel('Rainfall (in)');
title('Rainfall per year for 2014');
ylim([0 15])
subplot(4,1,3)
plot(weather_2019.Date,weather_2019.RainfallForMonth_in_,'Color',[0.5,0.1,0.5])
ylabel('Rainfall (in)');
title('Rainfall per month for 2019');
ylim([0 4])
subplot(4,1,4)
plot(weather_2019.Date,weather_2019.RainfallForYear_in_,'Color',[0.5,0.1,0.5])
ylabel('Rainfall (in)');
title('Rainfall per year for 2019');
ylim([0 15])

% Plotting Sunspots
fig_8=figure(8);
t = tiledlayout(2,1,'TileSpacing','Compact');
subplot(2,1,1)
plot(weather_2014.Date,weather_2014.DailySunspots,'Color',[0.7,0.2,0.2])
ylabel('Number of sunspots');
title('Daily number of sunspots for 2014');
ylim([0 (max(weather_2014.DailySunspots)+1)])
hold off
subplot(2,1,2)
plot(weather_2019.Date,weather_2019.DailySunspots,'Color',[0.7,0.2,0.2])
ylabel('Number of sunspots');
title('Daily number of sunspots for 2019');
ylim([0 (max(weather_2014.DailySunspots)+1)])
hold off


% converting timetables to tables
weather_2014=timetable2table(weather_2014);
weather_2019=timetable2table(weather_2019);

%calculating and plotting autocorrelation for 2014 values
fig_9=figure(9);
for i=2:21
    [acf,lags]=autocorr(weather_2014(:,i));
    subplot(5,4,(i-1))
    plot(acf.Lags,acf.ACF,'Color',[0.5,0.1,0.5])
    hold on
    plot(acf.Lags,((zalpha/sqrt(n14))*ones(21,1)),'c--','Color',[0.9,0.5,0.9])
    plot(acf.Lags,(-(zalpha/sqrt(n14))*ones(21,1)),'c--','Color',[0.9,0.5,0.9])
    xlabel('lag')
    ylabel('Autocorrelation')
    title(weather_2014.Properties.VariableNames{i})
    ylim([-1 1])
    xlim([0 20])
    hold off
end


% calculating and plotting autocorrelation for 2014 values after AR fitting
fig_10=figure(10);
for i=2:21
    arr=table2array(weather_2014(:,i));
    v=fitAR(arr,p);
    [acf,lags]=autocorr(v);
    v=array2table(v);
    weather_2014(:,i)=v;
    subplot(5,4,(i-1))
    plot(lags,acf,'Color',[0.5,0.1,0.5])
    hold on
    plot(lags,((zalpha/sqrt(n14))*ones(21,1)),'c--','Color',[0.9,0.5,0.9])
    plot(lags,(-(zalpha/sqrt(n14))*ones(21,1)),'c--','Color',[0.9,0.5,0.9])
    xlabel('lag')
    ylabel('Autocorrelation')
    title(weather_2014.Properties.VariableNames{i})
    ylim([-1 1])
    xlim([0 20])
    hold off
    j=j+1;
end

%calculating and plotting autocorrelation for 2014 values

fig_11=figure(11);
for i=2:21
    [acf,lags]=autocorr(weather_2019(:,i));
    subplot(5,4,(i-1))
    plot(acf.Lags,acf.ACF,'Color',[0.5,0.1,0.5])
    hold on
    plot(acf.Lags,((zalpha/sqrt(n19))*ones(21,1)),'c--','Color',[0.9,0.5,0.9])
    plot(acf.Lags,(-(zalpha/sqrt(n19))*ones(21,1)),'c--','Color',[0.9,0.5,0.9])
    xlabel('lag')
    ylabel('Autocorrelation')
    title(weather_2019.Properties.VariableNames{i})
    ylim([-1 1])
    xlim([0 20])
    hold off
end


% calculating and plotting autocorrelation for 2019 values after AR fitting
fig_12=figure(12);
for i=2:21
    arr=table2array(weather_2019(:,i));
    v=fitAR(arr,p);
    [acf,lags]=autocorr(v);
    v=array2table(v);
    weather_2019(:,i)=v;
    subplot(5,4,(i-1))
    plot(lags,acf,'Color',[0.5,0.1,0.5])
    hold on
    plot(lags,((zalpha/sqrt(n14))*ones(21,1)),'c--','Color',[0.9,0.5,0.9])
    plot(lags,(-(zalpha/sqrt(n14))*ones(21,1)),'c--','Color',[0.9,0.5,0.9])
    xlabel('lag')
    ylabel('Autocorrelation')
    title(weather_2019.Properties.VariableNames{i})
    ylim([-1 1])
    xlim([0 20])
    hold off
end


%plotting 2014 network 
arr=table2array(weather_2014(:,2:21));
[R,P]=corrcoef(arr);
fig_13=figure(13);
adj=(P<alpha);
graph1=graph(adj);
wr=R.';
w=abs(wr(triu(adj,1).'));
graph1.Edges.Weight=w;
graph1.Nodes.Name={weather_2014.Properties.VariableNames{2:21}}.';
LWidths = 5*graph1.Edges.Weight/max(graph1.Edges.Weight);
plot(graph1,'Layout','circle','LineWidth',LWidths,'NodeColor',[0.5,0.1,0.5],'EdgeColor',[0.7,0.3,0.9])
title('2014 weather network')

%plotting 2019 network
arr=table2array(weather_2019(:,2:21));
[R,P]=corrcoef(arr);
fig_14=figure(14);
adj=(P<alpha);
graph1=graph(adj);
wr=R.';
w=abs(wr(triu(adj,1).'));
graph1.Edges.Weight=w;
graph1.Nodes.Name={weather_2019.Properties.VariableNames{2:21}}.';
LWidths = 5*graph1.Edges.Weight/max(graph1.Edges.Weight);
plot(graph1,'Layout','circle','LineWidth',LWidths,'NodeColor',[0.5,0.1,0.5],'EdgeColor',[0.7,0.3,0.9])
title('2019 weather network')

% calculating and plotting Granger causality 2014
arr=table2array(weather_2014(:,2:21));
[GCIM,pGCIM] = GCI(arr,5,1);
fig_15=figure(15);
adj=(pGCIM<alpha);
graph1=digraph(adj);
wr=GCIM.';
w=abs(wr(adj.'));
graph1.Edges.Weight=w;
graph1.Nodes.Name={weather_2014.Properties.VariableNames{2:21}}.';
LWidths = 5*graph1.Edges.Weight/max(graph1.Edges.Weight);
plot(graph1,'Layout','circle','LineWidth',LWidths,'ArrowSize',10,'NodeColor',[0.5,0.1,0.5],'EdgeColor',[0.7,0.3,0.9])
title('2014 weather Granger causality network')

% calculating and plotting conditional Granger causality 2014
[CGCIM,pCGCIM] = CGCI(arr,5,1);
fig_16=figure(16);
adj=(pCGCIM<alpha);
graph1=digraph(adj);
wr=CGCIM.';
w=abs(wr(adj.'));
graph1.Edges.Weight=w;
graph1.Nodes.Name={weather_2014.Properties.VariableNames{2:21}}.';
LWidths = 5*graph1.Edges.Weight/max(graph1.Edges.Weight);
plot(graph1,'Layout','circle','LineWidth',LWidths,'ArrowSize',10,'NodeColor',[0.5,0.1,0.5],'EdgeColor',[0.7,0.3,0.9])
title('2014 weather conditional Granger causality network')

[CGCIM,pCGCIM] = CGCI(arr,5,1);
fig_17=figure(17);
adj=(pCGCIM<alpha2);
graph1=digraph(adj);
wr=CGCIM.';
w=abs(wr(adj.'));
graph1.Edges.Weight=w;
graph1.Nodes.Name={weather_2014.Properties.VariableNames{2:21}}.';
LWidths = 5*graph1.Edges.Weight/max(graph1.Edges.Weight);
plot(graph1,'Layout','circle','LineWidth',LWidths,'ArrowSize',10,'NodeColor',[0.5,0.1,0.5],'EdgeColor',[0.7,0.3,0.9])
title('2014 weather conditional Granger causality network (alpha=0.08)')

% calculating and plotting Granger causality 2019
arr=table2array(weather_2019(:,2:21));
[GCIM,pGCIM] = GCI(arr,5,1);
fig_18=figure(18);
adj=(pGCIM<alpha);
graph1=digraph(adj);
wr=GCIM.';
w=abs(wr(adj.'));
graph1.Edges.Weight=w;
graph1.Nodes.Name={weather_2019.Properties.VariableNames{2:21}}.';
LWidths = 5*graph1.Edges.Weight/max(graph1.Edges.Weight);
plot(graph1,'Layout','circle','LineWidth',LWidths,'ArrowSize',10,'NodeColor',[0.5,0.1,0.5],'EdgeColor',[0.7,0.3,0.9])
title('2019 weather Granger causality network')

% calculating and plotting conditional Granger causality 2019
[CGCIM,pCGCIM] = CGCI(arr,5,1);
fig_19=figure(19);
adj=(pCGCIM<alpha);
graph1=digraph(adj);
wr=CGCIM.';
w=abs(wr(adj.'));
graph1.Edges.Weight=w;
graph1.Nodes.Name={weather_2019.Properties.VariableNames{2:21}}.';
LWidths = 5*graph1.Edges.Weight/max(graph1.Edges.Weight);
plot(graph1,'Layout','circle','LineWidth',LWidths,'ArrowSize',10,'NodeColor',[0.5,0.1,0.5],'EdgeColor',[0.7,0.3,0.9])
title('2019 weather conditional Granger causality network')

[CGCIM,pCGCIM] = CGCI(arr,5,1);
fig_20=figure(20);
adj=(pCGCIM<alpha2);
graph1=digraph(adj);
wr=CGCIM.';
w=abs(wr(adj.'));
graph1.Edges.Weight=w;
graph1.Nodes.Name={weather_2019.Properties.VariableNames{2:21}}.';
LWidths = 5*graph1.Edges.Weight/max(graph1.Edges.Weight);
plot(graph1,'Layout','circle','LineWidth',LWidths,'ArrowSize',10,'NodeColor',[0.5,0.1,0.5],'EdgeColor',[0.7,0.3,0.9])
title('2019 weather conditional Granger causality network (alpha=0.08)')

% rescaling parameters
arr=table2array([weather_2014(:,2:15),weather_2014(:,17:21);weather_2019(:,2:15),weather_2019(:,17:21)]);
colmin=min(arr);
colmax=max(arr);
arr=rescale(arr,'InputMin',colmin,'InputMax',colmax);
weather_2014(:,2:15)=array2table(arr(1:365,1:14));
weather_2014(:,17:21)=array2table(arr(1:365,15:19));
weather_2019(:,2:15)=array2table(arr(366:end,1:14));
weather_2019(:,17:21)=array2table(arr(366:end,15:19));

% Regression
% Dependant variable minimum pressure col 16
% Ridge
model_2014=fitrlinear(weather_2014(:,2:21),'MinimumPressure','Learner','leastsquares','Regularization','ridge');
model_2019=fitrlinear(weather_2019(:,2:21),'MinimumPressure','Learner','leastsquares','Regularization','ridge');
X_2014=[weather_2014(:,2:15),weather_2014(:,16:21)];
X_2019=[weather_2019(:,2:15),weather_2019(:,16:21)];
Y_2014=weather_2014.MinimumPressure;
Y_2019=weather_2019.MinimumPressure;
Y_1414=predict(model_2014,X_2014);
Y_1419=predict(model_2014,X_2019);
Y_1914=predict(model_2019,X_2014);
Y_1919=predict(model_2019,X_2019);
Er_1414=abs((Y_1414-Y_2014));
Er_1419=abs((Y_1419-Y_2019));
Er_1914=abs((Y_1914-Y_2014));
Er_1919=abs((Y_1919-Y_2019));
MAE_1414=mean(Er_1414);
MAE_1419=mean(Er_1419);
MAE_1914=mean(Er_1914);
MAE_1919=mean(Er_1919);

fig_21=figure(21);
subplot(2,1,1)
plot(weather_2014.Date,Y_2014,'Color',[0.3,0.8,0.3])
hold on
plot(weather_2014.Date,Y_1414,'Color',[0.3,0.3,0.8])
plot(weather_2014.Date,Y_2019,'Color',[1,1,0])
plot(weather_2014.Date,Y_1419,'Color',[0.8,0.3,0.3])
legend('2014 real','2014 predicted','2019 real','2019 predicted')
ylabel('IID Minimum Pressure')
title('Predicted and real values of minimum pressure (model 2014, Ridge)')
hold off
subplot(2,1,2)
plot(weather_2014.Date,Er_1414,'Color',[0.3,0.3,0.8])
hold on
plot(weather_2014.Date,Er_1419,'Color',[0.8,0.3,0.3])
title('Error of the predicted minimum pressure values')
legend(sprintf('2014 MAE=%1.4f', MAE_1414),sprintf('2019 MAE=%1.4f', MAE_1419))
ylabel('Absolute error')
hold off

fig_22=figure(22);
subplot(2,1,1)
plot(weather_2014.Date,Y_2019,'Color',[0.3,0.8,0.3])
hold on
plot(weather_2014.Date,Y_1919,'Color',[0.3,0.3,0.8])
plot(weather_2014.Date,Y_2014,'Color',[1,1,0])
plot(weather_2014.Date,Y_1914,'Color',[0.8,0.3,0.3])
ylabel('IID Minimum Pressure')
title('Predicted and real values of minimum pressure (model 2019, Ridge)')
legend('2019 real','2019 predicted','2014 real','2014 predicted')
hold off
subplot(2,1,2)
plot(weather_2014.Date,Er_1919,'Color',[0.3,0.3,0.8])
hold on
plot(weather_2014.Date,Er_1914,'Color',[0.8,0.3,0.3])
title('Error of predicted minimum pressure values')
legend(sprintf('2019 MAE=%1.4f', MAE_1919),sprintf('2014 MAE=%1.4f', MAE_1914))
ylabel('Absolute error')
hold off

% lasso
model_2014=fitrlinear(weather_2014(:,2:21),'MinimumPressure','Learner','leastsquares','Regularization','lasso');
model_2019=fitrlinear(weather_2019(:,2:21),'MinimumPressure','Learner','leastsquares','Regularization','lasso');
X_2014=[weather_2014(:,2:15),weather_2014(:,16:21)];
X_2019=[weather_2019(:,2:15),weather_2019(:,16:21)];
Y_2014=weather_2014.MinimumPressure;
Y_2019=weather_2019.MinimumPressure;
Y_1414=predict(model_2014,X_2014);
Y_1419=predict(model_2014,X_2019);
Y_1914=predict(model_2019,X_2014);
Y_1919=predict(model_2019,X_2019);
Er_1414=abs((Y_1414-Y_2014));
Er_1419=abs((Y_1419-Y_2019));
Er_1914=abs((Y_1914-Y_2014));
Er_1919=abs((Y_1919-Y_2019));
MAE_1414=mean(Er_1414);
MAE_1419=mean(Er_1419);
MAE_1914=mean(Er_1914);
MAE_1919=mean(Er_1919);

fig_23=figure(23);
subplot(2,1,1)
plot(weather_2014.Date,Y_2014,'Color',[0.3,0.8,0.3])
hold on
plot(weather_2014.Date,Y_1414,'Color',[0.3,0.3,0.8])
plot(weather_2014.Date,Y_2019,'Color',[1,1,0])
plot(weather_2014.Date,Y_1419,'Color',[0.8,0.3,0.3])
legend('2014 real','2014 predicted','2019 real','2019 predicted')
ylabel('IID Minimum Pressure')
title('Predicted and real values of minimum pressure (model 2014, LASSO)')
hold off
subplot(2,1,2)
plot(weather_2014.Date,Er_1414,'Color',[0.3,0.3,0.8])
hold on
plot(weather_2014.Date,Er_1419,'Color',[0.8,0.3,0.3])
title('Error of the predicted minimum pressure values')
legend(sprintf('2014 MAE=%1.4f', MAE_1414),sprintf('2019 MAE=%1.4f', MAE_1419))
ylabel('Absolute error')
hold off

fig_24=figure(24);
subplot(2,1,1)
plot(weather_2014.Date,Y_2019,'Color',[0.3,0.8,0.3])
hold on
plot(weather_2014.Date,Y_1919,'Color',[0.3,0.3,0.8])
plot(weather_2014.Date,Y_2014,'Color',[1,1,0])
plot(weather_2014.Date,Y_1914,'Color',[0.8,0.3,0.3])
ylabel('IID Minimum Pressure')
title('Predicted and real values of minimum pressure (model 2019, LASSO)')
legend('2019 real','2019 predicted','2014 real','2014 predicted')
hold off
subplot(2,1,2)
plot(weather_2014.Date,Er_1919,'Color',[0.3,0.3,0.8])
hold on
plot(weather_2014.Date,Er_1914,'Color',[0.8,0.3,0.3])
title('Error of predicted minimum pressure values')
legend(sprintf('2019 MAE=%1.4f', MAE_1919),sprintf('2014 MAE=%1.4f', MAE_1914))
ylabel('Absolute error')
hold off

% Ridge with feature selection
X_2014=[weather_2014.MinimumTemperature__F_,weather_2014.MaximumHumidity___,weather_2014.AverageTemperature__F_];
X_2019=[weather_2019.AverageBarometer_in_,weather_2019.AverageHumidity___,weather_2019.AverageTemperature__F_,weather_2019.DailySunspots];
Y_2014=weather_2014.MinimumPressure;
Y_2019=weather_2019.MinimumPressure;
model_2014=fitrlinear(X_2014,Y_2014,'Learner','leastsquares','Regularization','ridge');
model_2019=fitrlinear(X_2019,Y_2019,'Learner','leastsquares','Regularization','ridge');
X_1419=[weather_2019.MinimumTemperature__F_,weather_2019.MaximumHumidity___,weather_2019.AverageTemperature__F_];
X_1914=[weather_2014.AverageBarometer_in_,weather_2014.AverageHumidity___,weather_2014.AverageTemperature__F_,weather_2014.DailySunspots];
Y_1414=predict(model_2014,X_2014);
Y_1419=predict(model_2014,X_1419);
Y_1914=predict(model_2019,X_1914);
Y_1919=predict(model_2019,X_2019);
Er_1414=abs((Y_1414-Y_2014));
Er_1419=abs((Y_1419-Y_2019));
Er_1914=abs((Y_1914-Y_2014));
Er_1919=abs((Y_1919-Y_2019));
MAE_1414=mean(Er_1414);
MAE_1419=mean(Er_1419);
MAE_1914=mean(Er_1914);
MAE_1919=mean(Er_1919);

fig_25=figure(25);
subplot(2,1,1)
plot(weather_2014.Date,Y_2014,'Color',[0.3,0.8,0.3])
hold on
plot(weather_2014.Date,Y_1414,'Color',[0.3,0.3,0.8])
plot(weather_2014.Date,Y_2019,'Color',[1,1,0])
plot(weather_2014.Date,Y_1419,'Color',[0.8,0.3,0.3])
legend('2014 real','2014 predicted','2019 real','2019 predicted')
ylabel('IID Minimum Pressure')
title('Predicted and real values of minimum pressure (model 2014, Ridge with feature selection)')
hold off
subplot(2,1,2)
plot(weather_2014.Date,Er_1414,'Color',[0.3,0.3,0.8])
hold on
plot(weather_2014.Date,Er_1419,'Color',[0.8,0.3,0.3])
title('Error of the predicted minimum pressure values')
legend(sprintf('2014 MAE=%1.4f', MAE_1414),sprintf('2019 MAE=%1.4f', MAE_1419))
ylabel('Absolute error')
hold off

fig_26=figure(26);
subplot(2,1,1)
plot(weather_2014.Date,Y_2019,'Color',[0.3,0.8,0.3])
hold on
plot(weather_2014.Date,Y_1919,'Color',[0.3,0.3,0.8])
plot(weather_2014.Date,Y_2014,'Color',[1,1,0])
plot(weather_2014.Date,Y_1914,'Color',[0.8,0.3,0.3])
ylabel('IID Minimum Pressure')
title('Predicted and real values of minimum pressure (model 2019, Ridge with feature selection)')
legend('2019 real','2019 predicted','2014 real','2014 predicted')
hold off
subplot(2,1,2)
plot(weather_2014.Date,Er_1919,'Color',[0.3,0.3,0.8])
hold on
plot(weather_2014.Date,Er_1914,'Color',[0.8,0.3,0.3])
title('Error of predicted minimum pressure values')
legend(sprintf('2019 MAE=%1.4f', MAE_1919),sprintf('2014 MAE=%1.4f', MAE_1914))
ylabel('Absolute error')
hold off
