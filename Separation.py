# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 19:32:03 2019

@author: Hangjian Zhao
"""

#This script is used for continous hydrograph separation over a long term period using digital filter approach
#1)Single parameter digital filter (Nathan and McMahon, 1990)
#2)Two parameters digital filter (Eckhardt,2005))

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


class Components():
    def __init__(self, Quickflow, Baseflow,Totalflow):
        self.Quickflow = Quickflow
        self.Baseflow = Baseflow
        self.Totalflow = Totalflow

def Singlefilter(flow,a,ratio):  #0.9-0.95
    #Read in the continous flow data
    flow=pd.read_csv(flow).values
    #Pre-allocate the matrices
    #Runoff component
    R=np.zeros((flow.shape[0],flow.shape[1]))
    #Baseflow component
    B=np.zeros((flow.shape[0],flow.shape[1]))
    #Setting the starting condtion (Assuming at the first time step, runoff is zero and baseflow is quantile from the flow data)
    R[0,:]=0
    B[0,:]=np.quantile(flow,q=ratio,axis=0)
    
    for k in range(flow.shape[1]):    
        for t in np.arange(0,(flow.shape[0]-1),1):
            R[t+1,k]=a*R[t,k]+(1+a)/2*(flow[t+1,k]-flow[t,k])
            if R[t+1,k]<0:
                R[t+1,k]=0
            elif R[t+1,k]>flow[t+1,k]:
                R[t+1,k]=flow[t+1,k]     
            B[t+1,k]=flow[t+1,k]-R[t+1,k]
    
    return Components(R,B,flow)

def Twofilter(flow,a,BFIMAX,ratio):
    #Read in the continous flow data
    flow=pd.read_csv(flow).values
    #Pre-allocate the matrices
    #Runoff component
    R=np.zeros((flow.shape[0],flow.shape[1]))
    #Baseflow component
    B=np.zeros((flow.shape[0],flow.shape[1]))
    #Setting the starting condtion (Assuming at the first time step, runoff is zero and baseflow is quantile from the flow data)
    R[0,:]=0
    B[0,:]=np.quantile(flow,q=ratio,axis=0)
    
    for k in range(flow.shape[1]):    
        for t in np.arange(0,(flow.shape[0]-1),1): 
            B[t+1,k]=((1-BFIMAX)*a*B[t,k]+(1-a)*BFIMAX*flow[t+1,k])/(1-a*BFIMAX)
            
            if B[t+1,k]>flow[t+1,k]:
                B[t+1,k]=flow[t+1,k]
            else:
                B[t+1,k]=B[t+1,k]
            R[t+1,k]=flow[t+1,k]-B[t+1,k]
    
    return Components(R,B,flow)

def plotting(Results,method):
    f,(ax1,ax2)=plt.subplots(2,figsize=(6.5,8.5))
    ax1.plot(Results.Totalflow,lw=1.0)
    ax1.plot(Results.Baseflow,lw=1.5,c='red')
    ax1.set_title('Total and Baseflow-%s Parameter Digital Filter'%method)
    ax2.plot(Results.Quickflow)
    ax2.set_title('Quickflow-%s Parameter Digital Filter'%method)


#%%Funcation call
#Specify the filepath for the flow data
flow=r"\TotalFlow.csv"
#Baseflow filter patameter (default a=0.98)
a=0.95
#Quntiles of the flow to be used as the first day baseflow (setting up the initial condition)
ratio=0.01
#Specify maximum value of long term ratio of baseflow to total flow -addtional parameter in the Two parameters digital filter approach
BFIMAX=0.8  #0.80 for perennial streams with porous aquifers #0.50 for ephemeral streams with porous aquifers #0.25 for perennial streams with hard rock aquifers

#Calling single digital filter separation method
FlowSingle=Singlefilter(flow,a,ratio)

#Calling two parameters filter separation method
FlowDouble=Twofilter(flow,a,BFIMAX,ratio)

#Extracting the baseflow and quickflow component
QF1=FlowSingle.Quickflow
BF1=FlowSingle.Baseflow

QF2=FlowDouble.Quickflow
BF2=FlowDouble.Baseflow

#Plotting
plotting(FlowSingle,'Single')
plotting(FlowDouble,'Two')


#References:
#Nathan, R.J. and T.A. McMahon, 1990. Evaluation of Automated Techniques for Baseflow and Recession Analysis. Water Resources Research, 26(7):1465-1473.
#Eckhardt, K., 2005. How to Construct Recursive Digital Filters for Baseflow Separation. Hydrological Processes, 19(2):507-515.









































    
    
    
    
    
    