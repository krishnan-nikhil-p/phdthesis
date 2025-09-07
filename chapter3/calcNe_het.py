Nbs = np.array([75,100,125,150])
alphas = np.array([0.005,0.01,0.03])
taus = [[150,300,450,600,750,900,1050,1200] ,[75,150,225,300,375,450,525,600],[25,50,75,100,125,150,175,200]]

dataNe_lin_simpleT =  np.zeros((len(Nbs),len(taus[0]),len(alphas)))
dataNe_exp_simpleT =  np.zeros((len(Nbs),len(taus[0]),len(alphas)))
dataNe_lin_trueT =  np.zeros((len(Nbs),len(taus[0]),len(alphas)))
dataNe_exp_trueT = np.zeros((len(Nbs),len(taus[0]),len(alphas)))
trueTs =  np.zeros((len(Nbs),len(taus[0]),len(alphas)))
simpleTs = np.zeros((len(Nbs),len(taus[0]),len(alphas)))

for N, Nb in enumerate(Nbs):
    for a, alphaB0 in enumerate(alphas):
        for t, tau in enumerate(taus[a]):




                Tau = tau
                ti = np.argmin(np.abs(het_arr[N,t,a,:]/cnts[N,t,a] - .1))
                tf = np.argmin(np.abs(surv[N,t,a,:]-.05))
                #print(ti,tf)
                if (tf - ti <= 400):
                #    #ti= max(0, tf- 400)
                    ti= max(0, tf- 400)

                TgenS = 1/alphaB0 + Tau
                TgenC = np.sum(findGenTime(bprofs[N,t,a],pprofs[N,t,a], Nb,  alphaB0/Nb,Tau))
                m,b = np.polyfit(np.arange(0,tf*dts[N,t,a]-ti*dts[N,t,a],dts[N,t,a])/TgenS,np.log(het_arr[N,t,a,ti:tf]/cnts[N,t,a]),1) 
                Ne = - 1/m
                #print(-1/m)
                if a==2 and Nb ==125:
                    plt.plot(np.arange(0,tf*dts[N,t,a]-ti*dts[N,t,a],dts[N,t,a])/TgenC,  
                                                               het_arr[N,t,a,ti:tf ]/cnts[N,t,a],label=t)
                dataNe_lin_simpleT[N,t,a] = Ne

                Ne = - 1/( np.polyfit(np.arange(0,tf*dts[N,t,a]-ti*dts[N,t,a],dts[N,t,a])/TgenC,np.log(het_arr[N,t,a,ti:tf]/cnts[N,t,a]),1) [0])
                dataNe_lin_trueT[N,t,a] = Ne

                popt, pcov = curve_fit(lambda t,a,b,c: a*np.exp(b*t)+c,np.arange(0,tf*dts[N,t,a]-ti*dts[N,t,a],dts[N,t,a])/TgenS,  
                                                               het_arr[N,t,a,ti:tf ]/cnts[N,t,a],  p0=(.5, -0.01,.01))




                dataNe_exp_simpleT[N,t]=-1/popt[1]

                popt, pcov = curve_fit(lambda t,a,b,c: a*np.exp(b*t)+c,np.arange(0,tf*dts[N,t,a]-ti*dts[N,t,a],dts[N,t,a])/TgenC,  
                                                               het_arr[N,t,a,ti:tf ]/cnts[N,t,a],  p0=(.5, -0.01,.01))


                dataNe_exp_trueT[N,t,a]=-1/popt[1]
                trueTs[N,t,a] = TgenC
                simpleTs[N,t,a] = TgenS
                print(tf-ti,N,t,a,alphaB0)


plt.yscale('log')
plt.legend()