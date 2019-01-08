# -*- coding: utf-8 -*-
"""
Autoregressive Model of Reservoir inflows for Lucky Peak Complex
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.tsa.api as sm
import seaborn as sns
from pandas.core import datetools
from scipy.stats import truncnorm
import os

#os.chdir("C:/Users/kendrakaiser/Dropbox/BSU/Python/Data")
os.chdir("/Users/kendrakaiser/Documents/GitRepos/ReservoirModeling/Data")
doy=pd.read_csv("DOY.csv")
df =pd.read_csv("BRB_reservoir_data_1997-2018_noleap.csv")
df['Date'] =pd.to_datetime(df['Date'])
df['Y']=df['Date'].dt.year
df['M']=df['Date'].dt.month
df['D']=df['Date'].dt.day


# Create a Water Year column for our pandas data frame. This is a pretty 
# simple thing to do, but may not be necessary if you're not dealing with
# discharge data. Here's how it goes:
# 1. Create an empty array that is full of zeros and equal in length to 
#    the number of days in the record
WY = np.zeros(len(df['Y'].index)) 
# 2. For those records where the month is less than 10, their associated 
#    year is the correct water year
WY[df['M'].values < 10] = df['Y'].values[df['M'].values < 10] 
# 3. For those records where the month is greater than or equal to 10, 
#    the correct water year is one more than the current calendar year
WY[df['M'].values >= 10] = df['Y'].values[df['M'].values >= 10] + 1
# 4. Save the water year array as a column in the pandas data frame, as an
#    integer
df['WY'] = WY.astype(int)
df['doy']= doy.astype(int)

# Print fist and last 10 records to make sure we did it right
qrows = np.concatenate([np.arange(0,10,1),np.arange(-11,-1,1)])
df.iloc[qrows,:]

#Stats on every DOY for a variable distribution to draw from 


#doy_ss = np.zeros(365)
doy_std= np.zeros(365)
doy_var= np.zeros(365)
doy_mean= np.zeros(365)

for i in range(0,365):
    doy_inflw = np.zeros(21) 
    doy_inflw = df["in_unreg"].values[doy["doy"].values == i+1]
    #for j in range(0, len(doy_inflw)):
        #doy_ss[i] = sum([doy_inflw[j]**2])
    doy_std[i]=np.std(doy_inflw)
    doy_var[i]=np.var(doy_inflw)
    doy_mean[i]=np.mean(doy_inflw)

#plt.plot(range(1,366), doy_ss)
plt.subplot(221)
plt.plot(range(1,366), doy_std)
plt.subplot(222)
plt.plot(range(1,366), doy_var)
plt.subplot(223)
plt.plot(range(1,366), doy_mean)
plt.subplot(224)
# Density Plot and Histogram of mean flows
sns.distplot(doy_inflw, hist=True, kde=True, 
             bins=30, color = 'darkblue', 
             hist_kws={'edgecolor':'black'},
             kde_kws={'linewidth': 3})

np.min(df["in_unreg"])

stats=pd.DataFrame()
stats['doy']=range(1,366)
stats['std']=doy_std.astype(int)
stats['mean']=doy_mean.astype(int)

#place std of each doy in the main df
df2=df.merge(stats, left_on='doy', right_on='doy', how='inner')
df2=df2.sort_values(by=['Date'])

#draw random samples from normal distribution 
smp= np.zeros([365,100])
for j in range(0,100):
    for i in range(0,365):
        smp[i,j]=np.random.normal(0, doy_std[i])
 
plt.plot(range(0,365), smp)



#segment the whole dataframe into a training dataset (WY 1999- 2008)
# and a test dataset (everything that after Water Year 2008).
df_train = df2[df2.WY < 2008]
df_test  = df2[df2.WY >= 2008]

#examine strucutre of partial autocorrelation function
Qt = pd.Series(df_train['in_unreg'].values,df_train['Date'].values)

Qpacf = sm.pacf(Qt)

plt.figure(figsize=(14,10))
plt.stem(Qpacf)
plt.ylabel('Partial Autocorrelation Function (PACF) [-]',fontsize=16)
plt.xlabel('Lag Distance [day]',fontsize=16)
ax = plt.gca()
ax.tick_params('both',labelsize=16)
plt.show()

from statsmodels.graphics.tsaplots import plot_pacf
f, axarr = plt.subplots(1, 1, figsize=(14,10))
_ = plot_pacf(Qt,method='ols',lags=40,ax=axarr.axes)
axarr.set_title('Statsmodel plot of PACF',fontsize=16)
axarr.set_ylabel('Partial Autocorrelation Function (PACF) [-]',fontsize=16)
axarr.set_xlabel('Lag Distance [day]',fontsize=16)
axarr = plt.gca()
axarr.tick_params('both',labelsize=16)


#Fit AR(1) Model to Data

Q_AR1_model = sm.AR(Qt).fit(1)
print(Q_AR1_model.params)


#Use developed model to make predictions of the test data
Qtrain = df_train['in_unreg'].values
Qtest  = df_test['in_unreg'].values

DatesTest = df_test['Date'].values

Qttm1 = np.concatenate([Qtrain[-2:-1],Qtest[0:-1]]) 
#slice notation: -1 (last item in array) [-2:] last two items in array - so this is taking the 2nd to last value of the training set and the adding it to the beginning of the test set

AR1_mu   = Q_AR1_model.params[0] 
AR1_phi1 = Q_AR1_model.params[1]

e= np.random.normal(0, df_test['std'].values)

def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    e=truncnorm(
            (low-mean)/sd, (upp-mean)/sd, loc=mean, scale=sd)
    return e.rvs()
    
e2= get_truncated_normal(mean=0, sd=df_test['std'].values, low = -df_test['in_unreg'].values, upp=60000)
err2=e2.rvs()


reps= 100
Qpred = np.zeros((len(df_test['std']), reps), dtype=float)
err= np.zeros((len(df_test['std'])), dtype=float)      
            

#is there a way to doe this so its a function instead of a for loop? 
sns.distplot(err, hist=True, kde=True, 
    bins=30, color = 'darkblue', 
    hist_kws={'edgecolor':'black'},
    kde_kws={'linewidth': 3})
#Predictive MODEL 
QAR1 = AR1_mu + AR1_phi1*Qttm1
QhatAR1= QAR1*err

for n in range(0, len(df_test['std'])):
    err[n]= get_truncated_normal(mean=1, sd=.1, low = 0.7, upp=1.3)
    for j in range(0, reps): 
        Qpred[:,j]=(AR1_mu + AR1_phi1*Qttm1)*err

plt.figure(figsize=(14,10))

plt.plot(DatesTest,Qpred)
plt.plot(DatesTest,QhatAR1,'b-')

plt.plot(DatesTest,Qtest,'r-')
plt.ylabel('Discharge [m${}^3$/s]',fontsize=16)
plt.xlabel('Date',fontsize=16)
plt.legend(('AR(1)','Obs.'),fontsize=16)
ax = plt.gca()
ax.tick_params('both',labelsize=16)
plt.savefig('Inflow_AR1.pdf')
plt.show()

plt.figure(figsize=(14,10))
plt.plot(DatesTest,QhatAR1,'b-')

plt.plot(DatesTest,Qtest,'r-')
plt.xlim([DatesTest[250],DatesTest[300]])
plt.ylabel('Discharge [m${}^3$/s]',fontsize=16)
plt.xlabel('Date',fontsize=16)
plt.legend(('AR(1)','Obs.'),fontsize=16)
ax = plt.gca()
ax.tick_params('both',labelsize=16)
plt.show()

#compare predictions and observations
# Compute the R^2 values for each prediction 
R2AR1  = np.corrcoef(QhatAR1,Qtest)**2
muAR1  = np.mean(QhatAR1-Qtest)
stdAR1 = np.std(QhatAR1-Qtest)


# Plot the AR(1) and AR(2) model results 
plt.figure(figsize=(14,10))

plt.plot(QhatAR1,Qtest,'bo')
plt.plot([20, 80], [20, 80], 'k-')
plt.title('AR(1) Model',fontsize=16)
plt.ylabel('Observed Discharge [m${}^3$/s]',fontsize=16)
plt.xlabel('Observed Discharge [m${}^3$/s]',fontsize=16)
ax = plt.gca()
ax.tick_params('both',labelsize=16)
ax.annotate('R${}^2$ = %.3f'%R2AR1[0,1], xy=(20,80), fontsize=16)
ax.annotate('Avg. error = %.3f'%muAR1, xy=(20,78), fontsize=16)
ax.annotate('Std. error = %.3f'%stdAR1, xy=(20,76), fontsize=16)

plt.show()

