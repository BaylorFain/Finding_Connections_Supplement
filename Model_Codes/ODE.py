import numpy as np
import gc

np.random.seed(903)
###############################################

# Model parameters
Vo = 7.5*pow(10,-2)   #TCID_50/ml 
beta = 3.2*pow(10,-5) #((TCID_50/ml)*d)^-1
c = 5.2               #d^-1
p = 4.6*pow(10,-2)    #(TCID_50*ml)*d^-1
k = 4.0               #d^-1
delta = 5.2           #d^-1

MOI = Vo/(4.0*pow(10,8))

#Initial Conditions 
T0 = 1
E0 = 0
I0 = 0
D0 = 0
V0 = MOI*10**2

beta = beta/24
c = c/24
p = p/24
k = k/24
delta = delta/24

beta = beta*(4.0*pow(10,8))
c = c
p = p
k = k
delta = delta

n_simulations = 1
dt = 0.005
endtime = 24*20

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
    timetime = 0
    while timetime < endtime:

        Rate[0] = beta*V[i][-1][-1]*T[i][-1][-1]
        Rate[1] = k*E[i][-1][-1]
        Rate[2] = delta*I[i][-1][-1]
        Rate[3] = p*I[i][-1][-1]
        Rate[4] = c*V[i][-1][-1]

        T[i][0] = np.append(T[i][0], T[i][-1][-1] - dt*Rate[0])
        E[i][0] = np.append(E[i][0], E[i][-1][-1] + dt*Rate[0] - dt*Rate[1])
        I[i][0] = np.append(I[i][0], I[i][-1][-1] + dt*Rate[1] - dt*Rate[2])
        D[i][0] = np.append(D[i][0], D[i][-1][-1] + dt*Rate[2])
        V[i][0] = np.append(V[i][0], V[i][-1][-1] + dt*(Rate[3] - Rate[4]))
        Time[i][0] = np.append(Time[i][0], Time[i][-1][-1] + dt)
        timetime = Time[i][-1][-1]
        

##################################################
Time_ave = Time[0][0]
T_ave = T.mean(axis=0)[0]
E_ave = E.mean(axis=0)[0]
I_ave = I.mean(axis=0)[0]
D_ave = D.mean(axis=0)[0]
V_ave = V.mean(axis=0)[0]

with open('Data/ODE_Model/ODE.txt', 'w') as outfile:
    for i in range(len(Time_ave)):
        if i%int(1/dt) == 0: 
            outfile.write(str(T_ave[i])+', '+str(E_ave[i])+', '+str(I_ave[i])+', '+str(D_ave[i])+', '+str(V_ave[i])+', '+'\n')

gc.collect()



