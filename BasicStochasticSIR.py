#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import random
import matplotlib.pyplot as plt
import statistics as st

np.set_printoptions(precision=2)


# In[32]:


##Quarentine Measures
#Quarantine function step function

def ExpQuarentine(d,lockdownday=20,effective=0.8, speed=0.05):
    dd=d
    val=0
    if d< lockdownday:
        val= 1 
    else:
        val= ((1-effective)+ effective*np.exp(- (d-lockdownday)*speed))
    
    return val


# In[53]:


##Plotting Tools
#Set 

#Plot with variance
def var_plot(traj,ax1,last_d=None,color=[0, 102, 20]):
    plt.rcParams.update({'font.size': 16})
    
    rc=color[0]/255
    gc=color[1]/255
    bc=color[2]/255
    meantraj= curve_mean(traj)
    stdtraj = curve_std(traj)
    if last_d==None:
        last_d = min([len(meantraj),len(stdtraj)])
    x = np.arange(1, last_d, 1)
    y1 =np.asarray(meantraj[1:last_d]) + np.asarray(stdtraj[1:last_d])
    y2 =np.asarray(meantraj[1:last_d]) - np.asarray(stdtraj[1:last_d])
    y3= np.asarray(meantraj[1:last_d])

    ax1.fill_between(x, y1, y2, facecolor=(rc,gc,bc,0.3))
    ax1.plot(x,y3, color=(rc,gc,bc,1),linewidth=1)
    ax1.set_ylabel('Count')
    ax1.set_xlabel('Days')

    ratio = 0.7
    xleft, xright = ax1.get_xlim()
    ybottom, ytop = ax1.get_ylim()
    ax1.set_aspect(abs((xright-xleft)/(ybottom-ytop))*ratio)
    
#Plot many
def plot_many_trajectories(trajs, last_d=None):
    if last_d == None:
        last_d=len(trajs[0])
    plt.figure(figsize = (20,5));
    for t in trajs:
        plt.plot(range(1,last_d),t[1:last_d])


# In[4]:


##Statistics Functions
def curve_mean(traj):
    return  [sum(col) / float(len(col)) for col in zip(*traj)]
def curve_std(traj):
    traja=np.asarray(traj)
    cmean=curve_mean(traj)
    dlen=len(cmean)
    
    std=[np.std(traja[:,i]) for i in range(0,dlen)]
    return std
    


# In[63]:


#Simulation: 
#Takes Initial Susceptible, Initial Infected, Spread per day per person, 
#Recovery per day per person, max sim days, quarentine function
def SIR_sim(S0, I0, SPDPP, RPDPP, d_end, quarentine):
    #Initialize
    S=S0
    I=I0
    R=0
    eventList=[1,2,3]

    S_d=[0 for i in range(0,d_end)]
    I_d=[0 for i in range(0,d_end)]
    R_d=[0 for i in range(0,d_end)]

    #SIR Loop

    #print("Running: ", end =" ")
    for d in range(1,d_end):
        stepsperday=100*int((SPDPP*S*I +RPDPP*I))+1
        spreadrate= (SPDPP/stepsperday)*quarentine[d]
        recrate=RPDPP/stepsperday
        #print(".", end =" ",sep=' ')

        for step in range(0,stepsperday):
            spreadprob=spreadrate*S*I
            recprob=recrate*I
            probevent=[1-spreadprob-recprob,spreadprob,recprob]
            #Error Catching
            if probevent[0]<0:
                print("negative prob! on day",d ," step ", step, ": ", probevent[0])
                print(S, I, R)
                break;

            event=random.choices(eventList,probevent)[0]
            if event==1:
                S=S
            elif event==2:
                S=S-1
                I=I+1
            elif event==3:
                I=I-1
                R=R+1

        #Record Values at end of day        
        S_d[d]=S
        I_d[d]=I
        R_d[d]=R

        #If all recovered, end simulation
        if R + S == S0+I0:
            d_end=d
            #print("\rPandemic Ends on day",d)
            break
            
    #print("\r----Final :",S,I,R)
    return S_d,I_d,R_d,d_end     

    
    


# In[64]:


##Sim Loop For Agregate Data

def MultiSIRloop(S0=10000,I0=1, SR=0.0001, RR=1/14, sims=10, maxdays=200, quar=None):
    if quar==None:
        quar=[1 for i in range(0, maxdays)]

    S_traj=[[0 for i in range(0,maxdays)] for j in range(0, sims)] 
    I_traj=[[0 for i in range(0,maxdays)] for j in range(0, sims)] 
    R_traj=[[0 for i in range(0,maxdays)] for j in range(0, sims)]
    daycount=[0 for j in range(0, sims)]

    for s in range(0,sims):
        S_traj[s], I_traj[s], R_traj[s],daycount[s]=SIR_sim(S0,I0,SR, RR, maxdays, quar)

    return S_traj, I_traj, R_traj, daycount, quar

def quick_means(S0=10000,I0=1, SR=0.0001, RR=1/14, sims=1, maxdays=200):
    S, I, R, days, quar=MultiSIRloop(S0,I0, SR, RR, sims, maxdays)
    return S, I, R

def quick_infected(S0=10000,I0=1, SR=0.0001, RR=1/14, sims=1, maxdays=200):
    S, I, R, days, quar=MultiSIRloop(S0,I0, SR, RR, sims, maxdays)
    return S, I, R


# In[74]:


md=200
my_quar=[ExpQuarentine(d) for d in range(0,md)]
S, I, R, days, quar= MultiSIRloop(sims=2, maxdays=md)


# In[75]:



last_days=80
fig_mine, (ax_mine) = plt.subplots(1, 1, sharex=True)
fig_mine.set_size_inches(18.5, 10.5)
    
var_plot(S,ax_mine, last_d=last_days,color=[61, 61, 142]) #blue
var_plot(I,ax_mine,last_d=last_days,color=[178, 76, 76]) #red
var_plot(R,ax_mine, last_d=last_days,color=[71, 132, 71]) #green

plt.plot(range(1,last_days),10000*np.asarray(quar.copy())[1:last_days], color='gray', linewidth=2);


# In[ ]:




