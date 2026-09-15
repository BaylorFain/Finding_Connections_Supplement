import gc
import numpy as np
import scipy.stats as stats 

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
beta = beta*(4.0*pow(10,8))
c = c
p = p
k = k
delta = delta

tauE = 1/k
nE = 30

tauI = 1/delta
nI = 100

endtime = 24*20
dt =  0.025
save = int(1/dt)
t = np.arange(0,endtime,dt)

inU = np.zeros([int(endtime/dt)])
U = np.zeros([int(endtime/dt)])
inE = np.zeros([int(endtime/dt)])
E = np.zeros([int(endtime/dt)])
inI = np.zeros([int(endtime/dt)])
I = np.zeros([int(endtime/dt)])
inD = np.zeros([int(endtime/dt)])
D = np.zeros([int(endtime/dt)])
V = np.zeros([int(endtime/dt)])

U[0] = 1.0
E[0] = 0.0
I[0] = 0.0
D[0] = 0.0
V[0] = MOI*10**2

print(beta, 'beta')
print(c, 'c')
print(p, 'p')
print(k, 'k')
print(delta, 'delta')
print(V[0])

exit(0)

for i in range(len(U)):

    if (i+1) == len(U):
        break

    MoveAmount = 0.0
    if (U[i] == 0.0) or (V[i] == 0.0):
        ProbU = np.zeros([len(t)])
    else:
        ProbU = stats.expon.pdf((t-(i*dt)),scale=(1/(beta*U[i]*V[i])))
    inU = ProbU*dt
    outU = min(inU[i],U[i])
    U[i+1] = U[i] + MoveAmount - outU       

    MoveAmount = outU
    ProbE = stats.gamma.pdf((t-(i*dt)), a=nE, scale=(tauE/nE))
    inE = inE + MoveAmount*ProbE*dt

    outE = min(inE[i],E[i])
    E[i+1] = E[i] + MoveAmount - outE

    MoveAmount = outE
    ProbI = stats.gamma.pdf((t-(i*dt)), a=nI, scale=(tauI/nI))
    inI = inI + MoveAmount*ProbI*dt

    outI = min(inI[i],I[i])
    I[i+1] = I[i] + MoveAmount - outI

    MoveAmount = outI
    ProbD = 0.0
    inD = inD + MoveAmount*ProbD*dt
    outD = min(inD[i],D[i])
    D[i+1] = D[i] + MoveAmount - outD

    V[i+1] = V[i] + (p*I[i] - c*V[i])*dt

with open(r'Data/AS_Model/AGE.txt', 'w') as outfile:
    for i in range(len(U)):
        if i%save == 0.0:
            outfile.write(str(U[i])+', '+str(E[i])+', '+str(I[i])+', '+str(D[i])+', '+str(V[i])+', '+'\n')

gc.collect()
