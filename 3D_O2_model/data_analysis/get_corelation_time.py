from ising_lib import *
from numpy import *

print("\n\n Finding the correlation times of each file now ..\n")

def variance(X):
    return var(X)

def Average(X):
    return average(X)    
    
    
def analyze(temparature,L,E,M):
    u,u_err=jacknife(E,Average)
    mag,mag_err=jacknife(M,Average)
    sp_heat,sp_heat_err=jacknife(E,variance)
    sp_heat/=(temparature*temparature*L)
    sp_heat_err/=(temparature*temparature*L)
    sucep,sucep_err=jacknife(M,variance)
    sucep/=(temparature*L)
    sucep_err/=(temparature*L)
    M2=[i**2 for i in M]
    M4=[i**4 for i in M]
    m2,m2_err=jacknife(M2,Average)
    m4,m4_err=jacknife(M4,Average)
    return [mag,mag_err],[sp_heat,sp_heat_err],[sucep,sucep_err],[m2,m2_err],[m4,m4_err],[u,u_err]

BASE_FOLDER='../multi_histo_data/'

print("reading file list at  : "+BASE_FOLDER+'fnames.txt')
f=open(BASE_FOLDER+'fnames.txt','r')
fnames=[]
l=f.readline()
while(l):
    fnames.append(l[:-1])
    l=f.readline()
f.close()

Max_lags=9000
fit_count_of_acf=Max_lags
tol=0.5

print("Number of files to be analyzed : ",len(fnames))
print("Max_lags = ",Max_lags)
print("fit_count_of_acf = ",fit_count_of_acf)
print("tol = ",tol)
print("\n writing config to : "+BASE_FOLDER+'hist.config')
f=open(BASE_FOLDER+'hist.config','w')
for nme in fnames:
    data=get_data(BASE_FOLDER+nme)
    tau,err,status=get_correlation_values(data['time'][-Max_lags:],data['energy'][-Max_lags:],
                            data['magentization'][-Max_lags:],Max_lags,fit_count_of_acf,tol)
    print(nme," -> ",tau," , ",err," , ",status)
    f.write(nme[:-6]+","+str(tau)+","+str(err)+","+str(status).lower()+'\n')
f.close()   
print(" done with the correlation thing !!")
