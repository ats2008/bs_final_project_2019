
# coding: utf-8

# In[6]:

from random import random,randint
from numpy import exp,var,average,sqrt
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from statsmodels.tsa.stattools import acf 

def __fu(x,a,b):
    return a*exp(-x/b)
def update_state(state,N,J,beta):
    sitex=randint(0,N-1)
    sitey=randint(0,N-1)
    newsp=state[sitex][sitey]*-1
    osp=state[sitex][sitey]

    if sitex==0:
        xm=N-1
    else:
        xm=sitex-1
    if sitey==0:
        ym=N-1
    else:
        ym=sitey-1
        
    if sitex==N-1:
        Xm=0
    else:
        Xm=sitex+1
    if sitey==N-1:
        Ym=0
    else:
        Ym=sitey+1

    de=(state[xm][sitey]+state[Xm][sitey]+
        state[sitex][ym]+state[sitex][Ym])*(newsp-state[sitex][sitey])*(-J)
    if de<=0:
        state[sitex][sitey]*=(-1)
    else:
        ex=exp(-de*beta)
        t=random()
        if t<ex:
            state[sitex][sitey]*=(-1)
        else:
            de=0
    dsp=state[sitex][sitey]-osp
    return state,de,dsp


# Initializes the Ising model with N x N lattice sites @T=inf state and returns [state,energy,magnetiztion]

# In[7]:

def initialize_HT(N,J):
    state=[[randint(0,1) for i in range(N)] for j in range(N)]
    for i in range(N):
        for j in range(N):
            if state[i][j]==0:
                state[i][j]=-1
    e=0
    for i in range(N):
        x1,x2=i-1,i+1
        x1=(x1+N)%N
        x2=(x2+N)%N
        for j in range(N):
            y1,y2=i-1,i+1
            y1=(y1+N)%N
            y2=(y2+N)%N
            e+=(state[x1][j]+state[x2][j]+state[i][y1]+state[i][y2])*state[i][j]*(-J)
    m=sum(sum(state[i][j] for i in range(N)) for j in range(N))
    return state,e,m   


# Initializes the Ising model with N x N lattice sites @T=0 state and returns [state,energy,magnetiztion]

# In[9]:

def initialize_ZT(N,J):
    state=[[-1 for i in range(N)] for j in range(N)]
    e=0
    for i in range(N):
        x1,x2=i-1,i+1
        x1=(x1+N)%N
        x2=(x2+N)%N
        for j in range(N):
            y1,y2=i-1,i+1
            y1=(y1+N)%N
            y2=(y2+N)%N
            e+=(state[x1][j]+state[x2][j]+state[i][y1]+state[i][y2])*state[i][j]*(-J)
    m=sum(sum(state[i][j] for i in range(N)) for j in range(N))
    return state,e,m


# Given fname returns the dictionary : {'params':params,'istate':istate,'fstate':fstate,'time':sw,'energy':en,'magentization':mag}

# In[10]:


# In[1]:

def get_data_OnModel(fname):
    f=open(fname,'r')
    l=f.readline()
    params=dict()
    sw,en,mag,magn=[],[],[],[]
    while "#"==l[0]:
        l=l[:-1].replace("#","")
#         print(l)
        it=l.split(",")
        if "sweep," in l:
        	for i in range(len(it)-3):
	        	magn.append([])
        	break
        try:
            params.update({it[0]:int(it[1])})
        except:
            try:
                params.update({it[0]:float(it[1])})
            except:
                params.update({it[0]:it[1]})
        l=f.readline()
    l=f.readline()   
    while l and "#"!=l[0]:
        it=l[:-1].split(",")
        sw.append(float(it[0]))
        en.append((float(it[1])))
        mag.append((float(it[2])))
        for t in range(len(it)-3):
            magn[t].append((float(it[3+t])))
        l=f.readline()
    istate=[]
    if "#istate" in l:
        l=f.readline()
        while l and "#"!=l[0]:
            it=l[:-1].split(",")
            row=[float(i) for i in it]
            istate.append(row)
            l=f.readline()
    fstate=[]
    if "#fstate" in l:
        l=f.readline()
        while l and "#"!=l[0]:
            it=l[:-1].split(",")
            row=[float(i) for i in it]
            fstate.append(row)
            l=f.readline()
    rslt={'params':params,'istate':istate,'fstate':fstate,'time':sw,'energy':en,'magentization':mag,
                                                         'magentization_n':magn}
    f.close()
    return rslt


# In[ ]:


def get_correlation_values(T,E,M,Max_lags,fit_count_of_acf=None,tol=0.15):
    """takes 'Temp','Energy','Magnetization' and return 'correlation time','its error','tolarence state' 
	Max_lags is maximum laging time to be considered [usually total length of the array given]
	fit_count_of_acf number of acf values to be given to fit function
	tol is fractional tolarence for the error 
    """
    c=acf(M,nlags=Max_lags,fft=True)
    if fit_count_of_acf==None:
        fit_count_of_acf=Max_lags
    K_fit=fit_count_of_acf
    a,b1=curve_fit(__fu,[i for i in range(len(c[:K_fit]))],c[:K_fit])
    t1=int(a[1])+1
    c=acf(E,nlags=Max_lags,fft=True)
#     print(len(c),len([i for i in range(len(c[:K_fit]))]),len(c[:K_fit]),K_fit)
    a,b2=curve_fit(__fu,[i for i in range(len(c[:K_fit]))],c[:K_fit])
    t2=int(a[1])+1
    
    if t1>t2:
        tau= t1
        err=sqrt(b1[1][1])
    else:
        tau=t2
        err=sqrt(b2[1][1])
    
    if err/tau > tol:
        return tau,err,False
    return tau,err,True
    


# In[ ]:

def sample_data(data,Max_lags=None,fit_count_of_acf=None,tol=0.15,return_params=False):
    if Max_lags==None:
        Max_lags=len(data['time'])-1 
    tau,err,status=get_correlation_values(data['time'][-Max_lags:],data['energy'][-Max_lags:],data['magentization'][-Max_lags:],Max_lags,fit_count_of_acf,tol)
    K=Max_lags
#    if not status:
#        return dict(),tau,err
    step=max(2*tau,25)
#     print(tau,err,step)
    L=len(data['time'])-1
    sample=dict()
    for i in data.keys():
        sample.update({i:[]})
    sample['params']=data['params']
    sample['istate']=data['istate']
    sample['fstate']=data['fstate'] 
    for i in range(int(K/step)):
        st=i*step
        sample['time'].append(data['time'][L-st])
        sample['energy'].append(data['energy'][L-st])
        sample['magentization'].append(data['magentization'][L-st])
    if return_params:
        return sample,tau,err
    else:
        return sample


# In[ ]:

def jacknife(A,func):
    c0=func(A)
    sig=0
    for i in range(len(A)):
        temp=A.pop(0)
        c=func(A)
        sig+=(c-c0)**2
        A.append(temp)
    return c0,sqrt(sig)


def write_data(fname,data):
    f=open(fname,'w')
    for i in data['params'].keys():
        f.write("#"+str(i)+" , "+str(data['params'][i])+"\n")
    istate,state=data['istate'],data['fstate']
    f.write("#sweep,E,M\n")
    for i,j,k in zip(data['time'],data['energy'],data['magentization']):
        f.write(str(i)+","+str(j)+","+str(k)+"\n")
    
    f.write("#istate \n")
    for i in range(len(istate)):
        f.write(str(istate[i][0]))
        for j in range(1,len(istate)):
            f.write(","+str(istate[i][j]))
        f.write("\n")
        
    f.write("#fstate \n")
    for i in range(len(state)):
        f.write(str(state[i][0]))
        for j in range(1,len(state)):
            f.write(","+str(state[i][j]))
        f.write("\n")
    f.close()

def stat_analysis(data,Max_lags):
    K_fit=Max_lags
    f,ax=plt.subplots(nrows=3,ncols=2,figsize=(10,15))
    ax[0][0].scatter(data['time'],data['energy'],s=1)
    ax[0][0].scatter(data['time'][-Max_lags:],data['energy'][-Max_lags:],c='g',s=2)
    ax[0][0].set_title('Energy vs Steps [whole]')
    ax[0][1].scatter(data['time'][-Max_lags:],data['energy'][-Max_lags:],s=2)
    ax[0][1].set_title('Energy vs Steps [Analyzed Part]')
    ax[1][0].scatter(data['time'],data['magentization'],s=1)
    ax[1][0].scatter(data['time'][-Max_lags:],data['magentization'][-Max_lags:],c='g',s=2)
    ax[1][0].set_title('Magnetization vs Steps [whole]')
    ax[1][1].scatter(data['time'][-Max_lags:],data['magentization'][-Max_lags:],s=2)
    ax[1][1].set_title('Magnetization vs Steps [Analyzed Part]')
    
    c=acf(data['energy'][-Max_lags:],nlags=Max_lags,fft=True)
    a,b=curve_fit(__fu,[i for i in range(len(c[:K_fit]))],c[:K_fit])
    m=max(c)
    C=[i/m for i in c]
    X=[i for i in range(len(c))]
    ax[2][0].scatter(X,C,s=2)
    ax[2][0].set_title('ACF for Energy')
    print("Energy tau : ",a[1],"\t\ttau err : ",sqrt(b[1][1]))
    
    c=acf(data['magentization'][-Max_lags:],nlags=Max_lags,fft=True)
    a,b=curve_fit(__fu,[i for i in range(len(c[:K_fit]))],c[:K_fit])
    m=max(c)
    C=[i/m for i in c]
    X=[i for i in range(len(c))]
    ax[2][1].scatter(X,C,s=2)
    ax[2][1].set_title('ACF for Magnetization')
    print("Magnetization tau : ",a[1],"\ttau err : ",sqrt(b[1][1]))

	
    

