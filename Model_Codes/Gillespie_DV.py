import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
import gc
import time,datetime


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

beta = beta*(4.0*pow(10,8))/T0
c = c
p = p
k = k
delta = delta

n_simulations = 10
endtime = 24*10
save = 1
savetime = np.arange(0,endtime+save,save)

#print(V0)
#print(beta)
#exit(0)

for i in range(n_simulations):
    with open('Data/Gillespie_Model/Deterministic/Gillespie_DV_'+str(i)+'.txt', 'w+') as outfile:
        outfile.write('')

Time = np.zeros([1],dtype=tuple)
T = np.zeros([1],dtype=tuple)
E = np.zeros([1],dtype=tuple)
I = np.zeros([1],dtype=tuple)
D = np.zeros([1],dtype=tuple)
V = np.zeros([1],dtype=tuple)

Rate = np.zeros((5))

num = 5
scale = num/V0

for i in range(n_simulations):
    np.random.seed(int(903+time.mktime(datetime.datetime.today().timetuple())))

    Time[0] = [0.0]
    T[0] = [T0]
    E[0] = [E0]
    I[0] = [I0]
    D[0] = [D0]
    V[0] = [V0*scale]

    print('simulation =',i)
    count = 0
    timetime = 0
    while timetime < endtime:

        Rate[0] = beta*V[-1][-1]/scale*T[-1][-1]
        Rate[1] = k*E[-1][-1]
        Rate[2] = delta*I[-1][-1]
        Rate[3] = p*I[-1][-1]*scale
        Rate[4] = c*V[-1][-1]

        R = Rate[0] + Rate[1] + Rate[2]
        u = np.random.random()
        tau = 1/R * np.log(1/u)

        NT = 0.0
        NE = 0.0
        NI = 0.0
        ND = 0.0
        NV = 0.0

        u = np.random.random()

        RandRateSum = u*R

        if ((RandRateSum >= 0)&(RandRateSum < Rate[0])):
            NT += -1
            NE +=  1
        elif ((RandRateSum >= Rate[0])&(RandRateSum < Rate[0]+Rate[1])):
            NE += -1
            NI +=  1
        elif ((RandRateSum >= Rate[0]+Rate[1])&(RandRateSum < Rate[0]+Rate[1]+Rate[2])):
            NI += -1
            ND +=  1

        NT  = max(NT, -T[-1][-1])
        NE  = max(NE, -E[-1][-1])
        NI  = max(NI, -I[-1][-1])
        
        T[0] = np.append(T[0], T[-1][-1] + NT)
        E[0] = np.append(E[0], E[-1][-1] + NE)
        I[0] = np.append(I[0], I[-1][-1] + NI)
        D[0] = np.append(D[0], D[-1][-1] + ND)
        V[0] = np.append(V[0], V[-1][-1] + tau*(Rate[3] - Rate[4]))

#        if (count < endtime):
#            if (np.floor(timetime) == savetime[count]):
#                with open('Data/Gillespie_Model/Deterministic/Gillespie_DV_'+str(i)+'.txt', 'a') as outfile:
#                        outfile.write(str(Time[-1][-1])+', '+str(T[-1][-1])+', '+str(E[-1][-1])+', '+str(I[-1][-1])+', '+str(D[-1][-1])+', '+str(V[-1][-1]/scale)+', '+'\n')
#                count += 1
        while count < len(savetime) and timetime >= savetime[count]:
            with open('Data/Gillespie_Model/Deterministic/Gillespie_DV_'+str(i)+'.txt', 'a') as outfile:
                outfile.write(f"{savetime[count]}, {T[-1][-1]}, {E[-1][-1]}, {I[-1][-1]}, {D[-1][-1]}, {V[-1][-1]/scale}\n")
            count += 1

        Time[0] = np.append(Time[0], Time[-1][-1] + tau)
        timetime = Time[-1][-1]

        if (T[-1][-1] + E[-1][-1] + I[-1][-1]) == 0.0:
            with open('Data/Gillespie_Model/Deterministic/Gillespie_DV_'+str(i)+'.txt', 'a') as outfile:
                outfile.write(str(Time[-1][-1])+', '+str(T[-1][-1])+', '+str(E[-1][-1])+', '+str(I[-1][-1])+', '+str(D[-1][-1])+', '+str(V[-1][-1]/scale)+', '+'\n')
            break

gc.collect()

