import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.stats import binom
import gc

import scipy.special as sc
from scipy.optimize import minimize
from scipy.optimize import minimize_scalar

import time,datetime
np.random.seed(int(903+time.mktime(datetime.datetime.today().timetuple())))
###############################################

# Model parameters
Vo = 7.5*pow(10,-2)   #TCID_50/ml 
beta = 3.2*pow(10,-5) #((TCID_50/ml)*d)^-1
c = 5.2               #d^-1
p = 4.6*pow(10,-2)    #(TCID_50*ml)*d^-1
k = 4.0               #d^-1
delta = 5.2           #d^-1

beta = beta/24
c = c/24
p = p/24
k = k/24
delta = delta/24

MOI = Vo/(4.0*pow(10,8))

#Initial Conditions 
T0 = 1001365
E0 = 0
I0 = 0
D0 = 0
V0 = MOI*T0*10**2

beta = beta*(4.0*pow(10,8))
c = c
p = p
k = k
delta = delta

print(beta, 'beta')
print(c, 'c')
print(p, 'p')
print(k, 'k')
print(delta, 'delta')
print(V0)

n_simulations = 1
tau =  0.001 #dt
save = int(1/tau)
endtime = 24*10 #maximum days of infection

Time = np.zeros([n_simulations, 1],dtype=tuple)
T = np.zeros([n_simulations, 1],dtype=tuple)
E = np.zeros([n_simulations, 1],dtype=tuple)
I = np.zeros([n_simulations, 1],dtype=tuple)
D = np.zeros([n_simulations, 1],dtype=tuple)
V = np.zeros([n_simulations, 1],dtype=tuple)

for i in range(n_simulations):
    Time[i][0] = [0.0]
    T[i][0] = [T0]
    E[i][0] = [E0]
    I[i][0] = [I0]
    D[i][0] = [D0]
    V[i][0] = [V0]

Rate = np.zeros((5))

for i in range(n_simulations):

    with open('Data/Tau_Leaping_Model/Stochastic/TauLeap_SV_'+str(i)+'.txt', 'w+') as outfile:
        outfile.write(str(Time[i][0][-1])+', '+str(T[i][0][-1])+', '+str(E[i][0][-1])+', '+str(I[i][0][-1])+', '+str(D[i][0][-1])+', '+str(V[i][0][-1])+', '+'\n')

    print('simulation =',i)
    timetime = 0
    while timetime < endtime:

        Rate[0] = beta*V[i][-1][-1]*T[i][-1][-1]/T0
        Rate[1] = k*E[i][-1][-1]
        Rate[2] = delta*I[i][-1][-1]
        Rate[3] = p*I[i][-1][-1]
        Rate[4] = c*V[i][-1][-1]

        leap = tau

#########################################################

        uT  = np.random.random()
        uE  = np.random.random()
        uI  = np.random.random()
        uVp = np.random.random()
        uVc = np.random.random()
        
        NT  = poisson.ppf(uT, Rate[0]*leap)
        NE  = poisson.ppf(uE, Rate[1]*leap)
        NI  = poisson.ppf(uI, Rate[2]*leap)

        def func(q):
            left = uVp
            right = sc.gammaincc(q, Rate[3]*leap)
            return(np.abs(left - right))
        resc = minimize_scalar(func,method='bounded',bounds=(0,Rate[3]*leap*2))
        NVp = resc.x

        def func(q):
            left = uVc
            right = sc.gammaincc(q, Rate[4]*leap)
            return(np.abs(left - right))
        resc = minimize_scalar(func,method='bounded',bounds=(0,Rate[4]*leap*2))
        NVc = resc.x

        NT  = min(NT, T[i][-1][-1])
        NE  = min(NE, E[i][-1][-1])
        NI  = min(NI, I[i][-1][-1])
        NVc = min(NVc, V[i][-1][-1])

        T[i][0] = np.append(T[i][0], T[i][-1][-1] - NT)
        E[i][0] = np.append(E[i][0], E[i][-1][-1] + NT - NE)
        I[i][0] = np.append(I[i][0], I[i][-1][-1] + NE - NI)
        D[i][0] = np.append(D[i][0], D[i][-1][-1] + NI)
        V[i][0] = np.append(V[i][0], V[i][-1][-1] + NVp - NVc)

        Time[i][0] = np.append(Time[i][0], Time[i][-1][-1] + leap)
        timetime = Time[i][-1][-1]

        print('[%d%%]\r'%(timetime/endtime*100), end="")

        if int((timetime)/tau)%save == 0:
            with open('Data/Tau_Leaping_Model/Stochastic/TauLeap_SV_'+str(i)+'.txt', 'a+') as outfile:
                outfile.write(str(Time[i][-1][-1])+', '+str(T[i][-1][-1])+', '+str(E[i][-1][-1])+', '+str(I[i][-1][-1])+', '+str(D[i][-1][-1])+', '+str(V[i][-1][-1])+', '+'\n')
##################################################


Time_ave = Time[0][0]
T_ave = T.mean(axis=0)[0]
E_ave = E.mean(axis=0)[0]
I_ave = I.mean(axis=0)[0]
D_ave = D.mean(axis=0)[0]
V_ave = V.mean(axis=0)[0]

with open('Data/Tau_Leaping_Model/Stochastic/TauLeap_SV_average.txt', 'w') as outfile:
    for i in range(len(Time_ave)):
        if i%save == 0:
            outfile.write(str(Time_ave[i])+', '+str(T_ave[i]/T0)+', '+str(E_ave[i]/T0)+', '+str(I_ave[i]/T0)+', '+str(D_ave[i]/T0)+', '+str(V_ave[i]/T0)+', '+'\n')

gc.collect()

################################################################

