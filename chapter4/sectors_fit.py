import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import numpy_indexed as npi
import seaborn as sns
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
from scipy import stats

folder= 'path_to_folder' ##folder with output of simulation



Ispace = np.array([10,20,30,40,50,60,70]) #initial radii
Ninf = [] #empty list for calculated sectors at long time

colors = plt.cm.viridis(np.linspace(0,1,7))


Aspace = np.array([1,1.5,2,2.5]) ###range of values of A examined


dat = []
fitBuff = [450,450,450,450]

A=2
##iniate empty lists for calculated values
data = []
if A % 1 ==0:
    A = int(A)
Ninf = []

I_exps = []

xmeans = []
ymeans = []
popts = []

fig, ax =plt.subplots()## iniatiate plot

for i, I in enumerate(Ispace):
    datax=[]
    datay= []
    sectf = glob.glob(folder+'sect*A*'+str(A)+'_*U'+str(I)+'*') ##sector data
    popf = glob.glob(folder+'pop*A*'+str(A)+'_*U'+str(I)+'*') ##population history
    readcnts=0
    readcntp=0
    I_exp=0
    for f in sectf:
        try:
            datax.append(list(np.loadtxt(f)[:,0])) ##"radius"
            datay.append(list(np.loadtxt(f)[:,1]*2)) ##sectors
            readcnts+=1 
        except:
            None

    vels =0
    pos = np.zeros(400)
    for f in popf:


        posDat = np.sqrt((np.loadtxt(f)[:,1]/50)/np.pi) ## calculate radius
        pos[:len(posDat)] += posDat

        I_exp+=posDat[0]
        readcntp+=1 


    if readcntp>0:
        pos /=readcntp 
        pos = pos[:(np.argwhere(np.diff(pos)<0)[0][0]-20)]
        t_vel = np.arange(0,len(pos)*5,5)

        vels = np.zeros(len(pos)-50)
        for q in range(50,len(pos)):
            #print(np.polyfit(t[i:],pos[i:],1)[0])
            vels[i-50] = np.polyfit(t_vel[i:],pos[i:],1)[0]

        vels = vels[:(np.argwhere(np.diff(vels)<0)[0][0]-1)]
        a_vel=np.mean(vels[-10:])

    else:
        a_vel = 0

    #print(a_vel)
    #plt.plot(t,vs)

    x = np.concatenate(datax[:])
    y = np.concatenate(datay[:])
    x_unique, y_mean = npi.group_by(x).mean(y)



    bump = -1


    I_exps.append(I_exp/readcntp) ### normalized sector count


    mPt = np.argwhere(x_unique>300)[0][0]
    mse=[]

    print(mPt)
    ###perform fit
    popt,popc = curve_fit(sectFit, x_unique[mPt:],y_mean[mPt:],bounds = ([0,-300],[50,0]))


    perr = np.sqrt(np.diag(popc))
    
    Ninf.append(sectFit(10**6,*popt) )
    #plt.scatter(I,sectFit(10**6,*popt))

    ax.plot(np.linspace(200,1000),sectFit(np.linspace(200,1000),*popt),color=colors[i],linestyle= '--')
    ax.plot(x_unique,y_mean,color=colors[i],label=I)
    ## save data
    data.append([I_exp/readcntp] + list(popt)+list(perr))
    popts.append(popt)
    xmeans.append(x_unique)
    ymeans.append(y_mean)

ax.legend()

###plot for  fitted sectors as in thesis

fig,axs= plt.subplots(1,4,figsize=(18,5))##iniiate plot

ax=axs[0]
ax.plot(xmeans[-1],ymeans[-1],'k',label='simulations',lw=5)
ax.plot(np.linspace(100,1000),sectFit(np.linspace(100,1000),*popts[-1]),color='gray',linestyle= '--',
        label=r'N = $\frac{at}{t+b}$ fit',lw=3)

#add only bottom and left spine        
for axis in ['top','right']:
    ax.spines[axis].set_visible(False)
for axis in ['bottom','left']:
    ax.spines[axis].set_linewidth(3)

#adjust tickparameters
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlabel('$r(t)$',fontsize=24)
ax.set_ylabel('$N(r)$',fontsize=24)
ax.legend(fontsize=16,framealpha=0,edgecolor='w')


ax=axs[1]
ax.scatter(np.sqrt(I_exps),Ninf,marker='x',s=200,label='simulations')


for axis in ['top','right']:
    ax.spines[axis].set_visible(False)
for axis in ['bottom','left']:
    ax.spines[axis].set_linewidth(3)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlabel('$\sqrt{r_0}$',fontsize=24)
ax.set_ylabel('$N(\infty)$',fontsize=24)
res = stats.linregress(np.sqrt(I_exps),Ninf)
ax.plot(np.linspace(2.5,9.5,51), res.slope*np.linspace(2.5,9.5,51)+res.intercept,c='gray',linestyle='--'
        ,label='linear fit')
ax.legend(fontsize=16,framealpha=0,edgecolor='w')


ax=axs[2]

ax.scatter([1,1.5,2,2.5],df['velocity'].values*df['Dg/v'].values,marker='x',s=200)


for axis in ['top','right']:
    ax.spines[axis].set_visible(False)
for axis in ['bottom','left']:
    ax.spines[axis].set_linewidth(3)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlabel('$A$',fontsize=24)

ax.set_ylabel('$D_g$',fontsize=24)


ax=axs[3]

ax.scatter([1,1.5,2,2.5],df['velocity'].values*df['Ds/v'].values,marker='x',s=200)


for axis in ['top','right']:
    ax.spines[axis].set_visible(False)
for axis in ['bottom','left']:
    ax.spines[axis].set_linewidth(3)
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xlabel('$A$',fontsize=24)

ax.set_ylabel('$D_s$',fontsize=24)


plt.subplots_adjust(wspace=.4)
plt.tight_layout()
#plt.savefig('Fig_sector_fit.pdf')

