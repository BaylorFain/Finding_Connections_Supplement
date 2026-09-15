
from scipy.integrate import odeint
from scipy.integrate import quad
import gc
import numpy as np
import datetime
import matplotlib
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from scipy.stats import shapiro
from scipy.stats import anderson
from scipy.stats import normaltest
from scipy.stats import ttest_1samp
import scipy.integrate as integrate
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from scipy.optimize import minimize
import scipy.stats as stats 
import warnings
import math
from scipy.special import gammainc
from scipy.special import gammaincc
from scipy.signal import deconvolve


startdex = 1

MOI = [0.0]

NumberOfRuns = 10

save = 1/1
xxspace = save
beta = 533.3
#beta = 0.0
c = 0.2
p = 0.0
D = 2.16*10**-8

Analysis = 1

#OuterPath = r"/home/baylor/Documents/Research/Papers/Current_Papers/REVIEW/Data/"
OuterPath = ""

BigFit = np.zeros([len(MOI), NumberOfRuns, 3],dtype=float )

BigPeakVirus = np.empty([len(MOI), NumberOfRuns],dtype=float )
BigTimeOfPeak = np.empty([len(MOI), NumberOfRuns],dtype=float )
BigUpSlope = np.empty([len(MOI), NumberOfRuns],dtype=float )
BigDownSlope = np.empty([len(MOI), NumberOfRuns],dtype=float )
BigEndTime = np.empty([len(MOI), NumberOfRuns],dtype=float )
BigAUC = np.empty([len(MOI), NumberOfRuns],dtype=float )
        
BigTimeArray = np.empty([NumberOfRuns],dtype=tuple )
BigHealthyArray = np.empty([NumberOfRuns],dtype=tuple )
BigEclipseArray = np.empty([NumberOfRuns],dtype=tuple )
BigInfectedArray = np.empty([NumberOfRuns],dtype=tuple )
BigDeadArray = np.empty([NumberOfRuns],dtype=tuple )
BigVirusArray = np.empty([NumberOfRuns],dtype=tuple )

ReallyBigTimeArray = np.empty([len(MOI),NumberOfRuns],dtype=tuple )
ReallyBigHealthyArray = np.empty([len(MOI),NumberOfRuns],dtype=tuple )
ReallyBigEclipseArray = np.empty([len(MOI),NumberOfRuns],dtype=tuple )
ReallyBigInfectedArray = np.empty([len(MOI),NumberOfRuns],dtype=tuple )
ReallyBigDeadArray = np.empty([len(MOI),NumberOfRuns],dtype=tuple )
ReallyBigVirusArray = np.empty([len(MOI),NumberOfRuns],dtype=tuple )

if Analysis == 1:
    for j in range(len(MOI)):

        lamb = np.array([1, 0.7, 0.8])
        amp = np.array([1, 0.01, 0.0001])

        LongestTime = 0.0
        for i in range(NumberOfRuns):    
            os.system('cls' if os.name=='nt' else 'clear')
            if((c != 0.0) and (p != 0.0)):
                print("c != 0.0, p != 0.0")
            elif(p != 0.0):
                print("p != 0.0")
            elif(c != 0.0):
                print("c != 0.0")
            else:
                print("beta != 0.0")       
            print(str(MOI[j])+"\n")
            print(i)

            InnerPath = r"Data/AB_Model/607_"+str(i)+"-CELLFREE_"+str(MOI[j])+"-MOI_"+str(beta)+"-beta_"+str(c)+"-c_"+str(p)+"-p_"
            
            INNERPATH = os.path.join(OuterPath+InnerPath).replace('\\','//')
            directory = os.path.join("".join(INNERPATH + "/Analysis"))
            if not os.path.exists(directory):
                os.makedirs(directory)
            Inner_Path_to_Folder = os.path.abspath(directory)

########################################################################

            TimeArray = []
            HealthyArray = []
            EclipseArray = []
            InfectedArray = []
            DeadArray = []
            VirusArray = []
            
            count = 0.0
            with open(os.path.join(INNERPATH,"PerTimeStep.txt")) as CellInFile:
                count1 = 0
                for cellline in CellInFile:
                    CellData = cellline
                    CellData = CellData.split(",")
                    del CellData[-1]
                    count = count + 1
                    if(LongestTime < count):
                        LongestTime = count
                        
                    TimeArray.append(float(count1*save))
                    HealthyArray.append(int(CellData[1]))
                    EclipseArray.append(int(CellData[2]))
                    InfectedArray.append(int(CellData[3]))
                    DeadArray.append(int(CellData[4]))
                    VirusArray.append(float(CellData[5]))
            
                    count1 = count1 + 1

            BigTimeArray[i] = np.array(TimeArray)
            BigHealthyArray[i] = np.array(HealthyArray)
            BigEclipseArray[i] = np.array(EclipseArray)
            BigInfectedArray[i] = np.array(InfectedArray)
            BigDeadArray[i] = np.array(DeadArray)
            BigVirusArray[i] = np.array(VirusArray)

            ReallyBigTimeArray[j][i] = BigTimeArray[i]
            ReallyBigHealthyArray[j][i] = BigHealthyArray[i]
            ReallyBigEclipseArray[j][i] = BigEclipseArray[i]
            ReallyBigInfectedArray[j][i] = BigInfectedArray[i]
            ReallyBigDeadArray[j][i] = BigDeadArray[i]
            ReallyBigVirusArray[j][i] = BigVirusArray[i]
    #################################################################################

    with open('beta-'+str(beta)+',c-'+str(c)+',p-'+str(p)+'.npy', 'wb') as outfile:
        np.save(outfile, ReallyBigTimeArray)
        np.save(outfile, ReallyBigHealthyArray)
        np.save(outfile, ReallyBigEclipseArray)
        np.save(outfile, ReallyBigInfectedArray)
        np.save(outfile, ReallyBigDeadArray)
        np.save(outfile, ReallyBigVirusArray)

#################################################################################
with open('beta-'+str(beta)+',c-'+str(c)+',p-'+str(p)+'.npy', 'rb') as infile:
    Time = np.load(infile, allow_pickle=True)
    Healthy = np.load(infile, allow_pickle=True)
    Eclipse = np.load(infile, allow_pickle=True)
    Infected = np.load(infile, allow_pickle=True)
    Dead = np.load(infile, allow_pickle=True)
    Virus = np.load(infile, allow_pickle=True)

################################################################################

for jj in range(len(MOI)):
    LongestTime = 0.0
    for ii in range(NumberOfRuns): 
        if len(Time[jj][ii]) > LongestTime:
            LongestTime = len(Time[jj][ii])

    MedianTimeArray = []
    MedianHealthyArray = []
    MedianEclipseArray = []
    MedianInfectedArray = []
    MedianDeadArray = []
    MedianVirusArray = []

    STDTimeArray = []
    STDHealthyArray = []
    STDEclipseArray = []
    STDInfectedArray = []
    STDDeadArray = []
    STDVirusArray = []

    InbetweenTimeArray = []
    InbetweenHealthyArray = []
    InbetweenEclipseArray = []
    InbetweenInfectedArray = []
    InbetweenDeadArray = []
    InbetweenVirusArray = []

    countt = 0
    for i in range(int(LongestTime)):
        for k in range(NumberOfRuns):      
            try:
                InbetweenTimeArray.append(Time[jj][k][i])
            except IndexError:
                InbetweenTimeArray.append(0.0)
                
            try:
                InbetweenHealthyArray.append(Healthy[jj][k][i])
            except IndexError:
                InbetweenHealthyArray.append(0.0)

            try:
                InbetweenEclipseArray.append(Eclipse[jj][k][i])
            except IndexError:
                InbetweenEclipseArray.append(0.0)

            try:
                InbetweenInfectedArray.append(Infected[jj][k][i])
            except IndexError:
                InbetweenInfectedArray.append(0.0)

            try:
                InbetweenDeadArray.append(Dead[jj][k][i])
            except IndexError:
                InbetweenDeadArray.append(0.0)

            try:
                InbetweenVirusArray.append(Virus[jj][k][i])
            except IndexError:
                InbetweenVirusArray.append(0.0)
        
        zeros = len(np.where(np.array(InbetweenVirusArray) == 0.0)[0])
        if (i > 5*24) and (zeros >= 3) and (countt == 0):
            index = i
            countt = countt + 1
#            print(i)
#            print(InbetweenVirusArray)
        MedianTimeArray.append(np.mean(InbetweenTimeArray))
        MedianHealthyArray.append(np.mean(InbetweenHealthyArray))
        MedianEclipseArray.append(np.mean(InbetweenEclipseArray))
        MedianInfectedArray.append(np.mean(InbetweenInfectedArray))
        MedianDeadArray.append(np.mean(InbetweenDeadArray))
        MedianVirusArray.append(np.mean(InbetweenVirusArray))

        STDTimeArray.append(np.std(InbetweenTimeArray))
        STDHealthyArray.append(np.std(InbetweenHealthyArray))
        STDEclipseArray.append(np.std(InbetweenEclipseArray))
        STDInfectedArray.append(np.std(InbetweenInfectedArray))
        STDDeadArray.append(np.std(InbetweenDeadArray))
        STDVirusArray.append(np.std(InbetweenVirusArray))
        
        InbetweenTimeArray.clear()
        InbetweenHealthyArray.clear() 
        InbetweenEclipseArray.clear() 
        InbetweenInfectedArray.clear() 
        InbetweenDeadArray.clear() 
        InbetweenVirusArray.clear()     

    endex = -1
    MedianTimeArray = np.array(MedianTimeArray)[startdex:endex]-startdex*save
    maxsoup = max(np.array(MedianHealthyArray)[startdex:endex])
    MedianHealthyArray = np.array(MedianHealthyArray)[startdex:endex]/1001365
    MedianEclipseArray = np.array(MedianEclipseArray)[startdex:endex]/1001365
    MedianInfectedArray = np.array(MedianInfectedArray)[startdex:endex]/1001365
    MedianDeadArray = np.array(MedianDeadArray)[startdex:endex]/1001365
    MedianVirusArray = np.array(MedianVirusArray)[startdex:endex]/1001365

#########################################################################
    TimeArray = MedianTimeArray
    HealthyArray = MedianHealthyArray
    EclipseArray = MedianEclipseArray
    InfectedArray = MedianInfectedArray
    DeadArray = MedianDeadArray
    VirusArray = MedianVirusArray

    STDTimeArray = np.array(STDTimeArray)/1001365
    STDHealthyArray = np.array(STDHealthyArray)/1001365
    STDEclipseArray = np.array(STDEclipseArray)/1001365
    STDInfectedArray = np.array(STDInfectedArray)/1001365
    STDDeadArray = np.array(STDDeadArray)/1001365
    STDVirusArray = np.array(STDVirusArray)/1001365

###############################################################################

    ODEHealthyArray = []
    ODEEclipseArray = []
    ODEInfectedArray = []
    ODEDeadArray = []
    ODEVirusArray = []
    
    count = 0.0
    with open(os.path.join(OuterPath,"Data/ODE_Model/ODE.txt")) as CellInFile:
        count1 = 0
        for cellline in CellInFile:
            CellData = cellline
            CellData = CellData.split(",")
            del CellData[-1]
            count = count + 1
            if(LongestTime < count):
                LongestTime = count
                
            ODEHealthyArray.append(float(CellData[0]))
            ODEEclipseArray.append(float(CellData[1]))
            ODEInfectedArray.append(float(CellData[2]))
            ODEDeadArray.append(float(CellData[3]))
            ODEVirusArray.append(float(CellData[4]))
    
            count1 = count1 + 1

    ODEHealthyArray = np.array(ODEHealthyArray)
    ODEEclipseArray = np.array(ODEEclipseArray)
    ODEInfectedArray = np.array(ODEInfectedArray)
    ODEDeadArray = np.array(ODEDeadArray)
    ODEVirusArray = np.array(ODEVirusArray)

###############################################################################
    
    IntHealthyArray = []
    IntEclipseArray = []
    IntInfectedArray = []
    IntDeadArray = []
    IntVirusArray = []
    
    count = 0.0
    with open(os.path.join(OuterPath,"Data/AS_Model/AGE.txt")) as CellInFile:
        count1 = 0
        for cellline in CellInFile:
            CellData = cellline
            CellData = CellData.split(",")
            del CellData[-1]
            count = count + 1
            if(LongestTime < count):
                LongestTime = count
                
            IntHealthyArray.append(float(CellData[0]))
            IntEclipseArray.append(float(CellData[1]))
            IntInfectedArray.append(float(CellData[2]))
            IntDeadArray.append(float(CellData[3]))
            IntVirusArray.append(float(CellData[4]))
    
            count1 = count1 + 1

    IntHealthyArray = np.array(IntHealthyArray)
    IntEclipseArray = np.array(IntEclipseArray)
    IntInfectedArray = np.array(IntInfectedArray)
    IntDeadArray = np.array(IntDeadArray)
    IntVirusArray = np.array(IntVirusArray)

#########################################################################        

    TauTimeArray = []
    TauHealthyArray = []
    TauEclipseArray = []
    TauInfectedArray = []
    TauDeadArray = []
    TauVirusArray = []
    
    count = 0.0
    with open(os.path.join(OuterPath,"Data/Tau_Leaping_Model/Deterministic/TauLeap_DV_average.txt")) as CellInFile:
        count1 = 0
        for cellline in CellInFile:
            CellData = cellline
            CellData = CellData.split(",")
            del CellData[-1]
            count = count + 1
            if(LongestTime < count):
                LongestTime = count

            TauTimeArray.append(float(CellData[0]))                
            TauHealthyArray.append(float(CellData[1]))
            TauEclipseArray.append(float(CellData[2]))
            TauInfectedArray.append(float(CellData[3]))
            TauDeadArray.append(float(CellData[4]))
            TauVirusArray.append(float(CellData[5]))
    
            count1 = count1 + 1

    TauTimeArray = np.array(TauTimeArray)
    TauHealthyArray = np.array(TauHealthyArray)
    TauEclipseArray = np.array(TauEclipseArray)
    TauInfectedArray = np.array(TauInfectedArray)
    TauDeadArray = np.array(TauDeadArray)
    TauVirusArray = np.array(TauVirusArray)

######################################################################### 

    BigGillespieTimeArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieHealthyArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieEclipseArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieInfectedArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieDeadArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieVirusArray = np.empty([NumberOfRuns],dtype=tuple )

    LongestTime = 0.0
    for i in range(NumberOfRuns):

        GillespieTimeArray = []
        GillespieHealthyArray = []
        GillespieEclipseArray = []
        GillespieInfectedArray = []
        GillespieDeadArray = []
        GillespieVirusArray = []
        
        count = 0.0
        with open(os.path.join(OuterPath,"Data/Gillespie_Model/Deterministic/Gillespie_"+str(i)+".txt")) as CellInFile:
#        with open(os.path.join(OuterPath,"GillespieModel/Tampered/Gillespie"+str(i)+".txt")) as CellInFile:
            count1 = 0
            for cellline in CellInFile:
                CellData = cellline
                CellData = CellData.split(",")
#                del CellData[-1]
                count = count + 1
#                if(LongestTime < count):
#                    LongestTime = count

                GillespieTimeArray.append(float(CellData[0]))
                GillespieHealthyArray.append(float(CellData[1]))
                GillespieEclipseArray.append(float(CellData[2]))
                GillespieInfectedArray.append(float(CellData[3]))
                GillespieDeadArray.append(float(CellData[4]))
                GillespieVirusArray.append(float(CellData[5]))
        
                if((LongestTime < float(CellData[0])) and (float(CellData[0]) < np.inf)):
                    LongestTime = float(CellData[0])

                count1 = count1 + 1

        BigGillespieTimeArray[i] = np.array(GillespieTimeArray)
        BigGillespieHealthyArray[i] = np.array(GillespieHealthyArray)
        BigGillespieEclipseArray[i] = np.array(GillespieEclipseArray)
        BigGillespieInfectedArray[i] = np.array(GillespieInfectedArray)
        BigGillespieDeadArray[i] = np.array(GillespieDeadArray)
        BigGillespieVirusArray[i] = np.array(GillespieVirusArray)

    BTGillespieTimeArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieHealthyArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieEclipseArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieInfectedArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieDeadArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieVirusArray = np.zeros([10,int(LongestTime)],dtype=float )

    for j in range(NumberOfRuns):
        for i in range(int(LongestTime)):
            index = np.where(np.floor(BigGillespieTimeArray[j]) == i)[0]
            if len(index) > 0:
                index = index[0]

                BTGillespieTimeArray[j][i] = BigGillespieTimeArray[j][index]
                BTGillespieHealthyArray[j][i] = BigGillespieHealthyArray[j][index]
                BTGillespieEclipseArray[j][i] = BigGillespieEclipseArray[j][index]
                BTGillespieInfectedArray[j][i] = BigGillespieInfectedArray[j][i]
                BTGillespieDeadArray[j][i] = BigGillespieDeadArray[j][index]
                BTGillespieVirusArray[j][i] = BigGillespieVirusArray[j][index]
            else:
                break

    GillespieTimeArray = np.arange(0,LongestTime-1,1)
    GillespieHealthyArray = np.mean(BTGillespieHealthyArray,axis=0)
    GillespieEclipseArray = np.mean(BTGillespieEclipseArray,axis=0)
    GillespieInfectedArray = np.mean(BTGillespieInfectedArray,axis=0)
    GillespieDeadArray = np.mean(BTGillespieDeadArray,axis=0)
    GillespieVirusArray = np.mean(BTGillespieVirusArray,axis=0)

#    for i in range(10):
#        print(BTGillespieVirusArray[i][int(4.5*24):int(5.5*24)])

#    BTGillespieVirusArray[BTGillespieVirusArray == 0.0] = np.nan
#    GillespieVirusArray = np.nanmean(BTGillespieVirusArray,axis=0)
#    print(GillespieVirusArray[int(4.5*24):int(5.5*24)])

######################################################################### 

    index = np.argwhere(TimeArray==max(TimeArray))[0][0]
#    index = 247

################################################################################
    size = 32
    capsizes = [50, 21, 50, 21, 50, 21, 50, 21, 50, 21]
    axessize = 2

    plt.rcParams['xtick.labelsize'] = size*0.9
    plt.rcParams['ytick.labelsize'] = size
    plt.rcParams['axes.labelsize'] = size
    plt.rc('legend', fontsize=size)
    plt.rc('legend', title_fontsize=size)
    plt.rcParams['figure.figsize'] = [12 , 8]
    plt.rcParams['xtick.major.size'] = axessize*4
    plt.rcParams['xtick.major.width'] = axessize
    plt.rcParams['xtick.minor.width'] = axessize
    plt.rcParams['ytick.major.width'] = axessize
    plt.rcParams['ytick.minor.width'] = axessize
    plt.rc('axes', linewidth=axessize)
    plt.rcParams['lines.marker'] = "."
    plt.rcParams['lines.linestyle'] = ""
    rotate = -45

    tautime = TauTimeArray
    Gillespietime = GillespieTimeArray

#ALL
#########################################################################3
    TauHealthyArray = TauHealthyArray#/1001365
    TauEclipseArray = TauEclipseArray#/1001365
    TauInfectedArray = TauInfectedArray#/1001365
    TauDeadArray = TauDeadArray#/1001365
    TauVirusArray = TauVirusArray#/1001365

    plt.plot(TimeArray[:index], HealthyArray[:index], '.-', label="AB")
    plt.plot(TimeArray[:index], ODEHealthyArray[:index], '.-', label="ODE")
    plt.plot(TimeArray[:index], IntHealthyArray[:index], '.-', label="Age Structured")
    plt.plot(tautime, TauHealthyArray, '.-', label="Tau-leap")
    plt.plot(Gillespietime, GillespieHealthyArray[:-1]/1001365, '.-', label="Gillespie")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Cells Over time")
    plt.ylabel("Fraction of Uninfected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/All_in_one_U.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], EclipseArray[:index], '.-', label="AB") 
    plt.plot(TimeArray[:index], ODEEclipseArray[:index], '.-', label="ODE")
    plt.plot(TimeArray[:index], IntEclipseArray[:index], '.-', label="Age Structured")
    plt.plot(tautime, TauEclipseArray, '.-', label="Tau-leap")
    plt.plot(Gillespietime, GillespieEclipseArray[:-1]/1001365, '.-', label="Gillespie")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Eclipse Over time")
    plt.ylabel("Fraction of Eclipse Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/All_in_one_E.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()


    plt.plot(TimeArray[:index], InfectedArray[:index], '.-', label="AB")
    plt.plot(TimeArray[:index], ODEInfectedArray[:index], '.-', label="ODE")
    plt.plot(TimeArray[:index], IntInfectedArray[:index], '.-', label="Age Structured")
    plt.plot(tautime, TauInfectedArray, '.-', label="Tau-leap")
    plt.plot(Gillespietime, GillespieInfectedArray[:-1]/1001365, '.-', label="Gillespie")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Infected Over time")
    plt.ylabel("Fraction of Infected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/All_in_one_I.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], DeadArray[:index], '.-', label="AB")
    plt.plot(TimeArray[:index], ODEDeadArray[:index], '.-', label="ODE")
    plt.plot(TimeArray[:index], IntDeadArray[:index], '.-', label="Age Structured")
    plt.plot(tautime, TauDeadArray, '.-', label="Tau-leap")
    plt.plot(Gillespietime, GillespieDeadArray[:-1]/1001365, '.-', label="Gillespie")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Dead Over time")
    plt.ylabel("Fraction of Dead of Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/All_in_one_D.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], VirusArray[:index], '.-', label="AB")
    plt.plot(TimeArray[:index], ODEVirusArray[:index], '.-', label="ODE")
    plt.plot(TimeArray[:index], IntVirusArray[:index], '.-', label="Age Structured")
    plt.plot(tautime, TauVirusArray, '.-', label="Tau-leap")
    plt.plot(Gillespietime, GillespieVirusArray[:-1]/1001365, '.-', label="Gillespie")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Virus Over time")
    plt.ylabel("Amount of Virus per Cell")
    plt.xlabel("Days")
    plt.yscale("log")
    plt.ylim( (10**-10,10**-2) )
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/All_in_one_V.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

#AB
################################################################
    plt.plot(TimeArray[:index], HealthyArray[:index], '.-', label="AB")
    for i in range(len(Time[0])):
        plt.plot(Time[0][i], Healthy[0][i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Cells Over time")
    plt.ylabel("Fraction of Uninfected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/AB_U.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], EclipseArray[:index], '.-', label="AB") 
    for i in range(len(Time[0])):
        plt.plot(Time[0][i], Eclipse[0][i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Eclipse Over time")
    plt.ylabel("Fraction of Eclipse Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/AB_E.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()


    plt.plot(TimeArray[:index], InfectedArray[:index], '.-', label="AB")
    for i in range(len(Time[0])):
        plt.plot(Time[0][i], Infected[0][i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Infected Over time")
    plt.ylabel("Fraction of Infected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/AB_I.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], DeadArray[:index], '.-', label="AB")
    for i in range(len(Time[0])):
        plt.plot(Time[0][i], Dead[0][i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Dead Over time")
    plt.ylabel("Fraction of Dead of Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/AB_D.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], VirusArray[:index], '.-', label="AB")
    for i in range(len(Time[0])):
        plt.plot(Time[0][i], Virus[0][i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Virus Over time")
    plt.ylabel("Amount of Virus per Cell")
    plt.xlabel("Days")
    plt.yscale("log")
    plt.ylim( (10**-10,10**-2) )
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/AB_V.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

#ODE
################################################################
    plt.plot(TimeArray[:index], ODEHealthyArray[:index], '.-', label="ODE")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Cells Over time")
    plt.ylabel("Fraction of Uninfected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/ODE_U.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()
 
    plt.plot(TimeArray[:index], ODEEclipseArray[:index], '.-', label="ODE")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Eclipse Over time")
    plt.ylabel("Fraction of Eclipse Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/ODE_E.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], ODEInfectedArray[:index], '.-', label="ODE")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Infected Over time")
    plt.ylabel("Fraction of Infected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/ODE_I.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], ODEDeadArray[:index], '.-', label="ODE")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Dead Over time")
    plt.ylabel("Fraction of Dead of Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/ODE_D.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], ODEVirusArray[:index], '.-', label="ODE")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Virus Over time")
    plt.ylabel("Amount of Virus per Cell")
    plt.xlabel("Days")
    plt.yscale("log")
    plt.ylim( (10**-10,10**-2) )
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/ODE_V.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

#Renewal
################################################################
    plt.plot(TimeArray[:index], IntHealthyArray[:index], '.-', label="Age Structured")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Cells Over time")
    plt.ylabel("Fraction of Uninfected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Renewal_U.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], IntEclipseArray[:index], '.-', label="Age Structured")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Eclipse Over time")
    plt.ylabel("Fraction of Eclipse Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Renewal_E.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], IntInfectedArray[:index], '.-', label="Age Structured")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Infected Over time")
    plt.ylabel("Fraction of Infected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Renewal_I.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], IntDeadArray[:index], '.-', label="Age Structured")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Dead Over time")
    plt.ylabel("Fraction of Dead of Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Renewal_D.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(TimeArray[:index], IntVirusArray[:index], '.-', label="Age Structured")
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Virus Over time")
    plt.ylabel("Amount of Virus per Cell")
    plt.xlabel("Days")
    plt.yscale("log")
    plt.ylim( (10**-10,10**-2) )
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Renewal_V.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

#Tau-leap
################################################################
    BigTauTimeArray = np.empty([NumberOfRuns],dtype=tuple )
    BigTauHealthyArray = np.empty([NumberOfRuns],dtype=tuple )
    BigTauEclipseArray = np.empty([NumberOfRuns],dtype=tuple )
    BigTauInfectedArray = np.empty([NumberOfRuns],dtype=tuple )
    BigTauDeadArray = np.empty([NumberOfRuns],dtype=tuple )
    BigTauVirusArray = np.empty([NumberOfRuns],dtype=tuple )
    
    for i in range(NumberOfRuns):

        TauTimeArray = []
        TauHealthyArray = []
        TauEclipseArray = []
        TauInfectedArray = []
        TauDeadArray = []
        TauVirusArray = []

        count = 0
        with open(os.path.join(OuterPath,"Data/Tau_Leaping_Model/Stochastic/TauLeap_SV_"+str(i)+".txt")) as CellInFile:
            count1 = 0
            for cellline in CellInFile:
                CellData = cellline
                CellData = CellData.split(",")
                del CellData[-1]
                count = count + 1
                if(LongestTime < count):
                    LongestTime = count

                TauTimeArray.append(float(CellData[0]))                
                TauHealthyArray.append(float(CellData[1]))
                TauEclipseArray.append(float(CellData[2]))
                TauInfectedArray.append(float(CellData[3]))
                TauDeadArray.append(float(CellData[4]))
                TauVirusArray.append(float(CellData[5]))
        
                count1 = count1 + 1

        BigTauTimeArray[i] = np.array(TauTimeArray)
        BigTauHealthyArray[i] = np.array(TauHealthyArray)
        BigTauEclipseArray[i] = np.array(TauEclipseArray)
        BigTauInfectedArray[i] = np.array(TauInfectedArray)
        BigTauDeadArray[i] = np.array(TauDeadArray)
        BigTauVirusArray[i] = np.array(TauVirusArray)

    plt.plot(np.mean(BigTauTimeArray), np.mean(BigTauHealthyArray)/1001365, '.-', label="Tau-leap")
    for i in range(NumberOfRuns):
        plt.plot(BigTauTimeArray[i], BigTauHealthyArray[i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Cells Over time")
    plt.ylabel("Fraction of Uninfected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/tau_U.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(np.mean(BigTauTimeArray), np.mean(BigTauEclipseArray)/1001365, '.-', label="Tau-leap")
    for i in range(NumberOfRuns):
        plt.plot(BigTauTimeArray[i], BigTauEclipseArray[i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Eclipse Over time")
    plt.ylabel("Fraction of Eclipse Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/tau_E.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(np.mean(BigTauTimeArray), np.mean(BigTauInfectedArray)/1001365, '.-', label="Tau-leap")
    for i in range(NumberOfRuns):
        plt.plot(BigTauTimeArray[i], BigTauInfectedArray[i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Infected Over time")
    plt.ylabel("Fraction of Infected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/tau_I.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(np.mean(BigTauTimeArray), np.mean(BigTauDeadArray)/1001365, '.-', label="Tau-leap")
    for i in range(NumberOfRuns):
        plt.plot(BigTauTimeArray[i], BigTauDeadArray[i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Dead Over time")
    plt.ylabel("Fraction of Dead of Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/tau_D.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(np.mean(BigTauTimeArray), np.mean(BigTauVirusArray)/1001365, '.-', label="Tau-leap")
    for i in range(NumberOfRuns):
        plt.plot(BigTauTimeArray[i], BigTauVirusArray[i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Virus Over time")
    plt.ylabel("Amount of Virus per Cell")
    plt.xlabel("Days")
    plt.yscale("log")
    plt.ylim( (10**-10,10**-2) )
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/tau_V.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

#Gillespie
################################################################
    BigGillespieTimeArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieHealthyArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieEclipseArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieInfectedArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieDeadArray = np.empty([NumberOfRuns],dtype=tuple )
    BigGillespieVirusArray = np.empty([NumberOfRuns],dtype=tuple )

    LongestTime = 0.0
    for i in range(NumberOfRuns):

        GillespieTimeArray = []
        GillespieHealthyArray = []
        GillespieEclipseArray = []
        GillespieInfectedArray = []
        GillespieDeadArray = []
        GillespieVirusArray = []
        
        count = 0.0
        with open(os.path.join(OuterPath,"Data/Gillespie_Model/Stochastic/Gillespie_"+str(i)+".txt")) as CellInFile:
#        with open(os.path.join(OuterPath,"GillespieModel/Tampered/Gillespie"+str(i)+".txt")) as CellInFile:
            count1 = 0
            for cellline in CellInFile:
                CellData = cellline
                CellData = CellData.replace("\n",",")                
                CellData = CellData.split(",")
                del CellData[-1]
                count = count + 1
#                if(LongestTime < count):
#                    LongestTime = count

                GillespieTimeArray.append(float(CellData[0]))
                GillespieHealthyArray.append(float(CellData[1]))
                GillespieEclipseArray.append(float(CellData[2]))
                GillespieInfectedArray.append(float(CellData[3]))
                GillespieDeadArray.append(float(CellData[4]))
                GillespieVirusArray.append(float(CellData[5]))
        
                if(LongestTime < float(CellData[0])):
                    LongestTime = float(CellData[0])

                count1 = count1 + 1

        BigGillespieTimeArray[i] = np.array(GillespieTimeArray)
        BigGillespieHealthyArray[i] = np.array(GillespieHealthyArray)
        BigGillespieEclipseArray[i] = np.array(GillespieEclipseArray)
        BigGillespieInfectedArray[i] = np.array(GillespieInfectedArray)
        BigGillespieDeadArray[i] = np.array(GillespieDeadArray)
        BigGillespieVirusArray[i] = np.array(GillespieVirusArray)

    BTGillespieTimeArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieHealthyArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieEclipseArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieInfectedArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieDeadArray = np.zeros([10,int(LongestTime)],dtype=float )
    BTGillespieVirusArray = np.zeros([10,int(LongestTime)],dtype=float )

    for j in range(NumberOfRuns):
        for i in range(int(LongestTime)):
            index = np.where(np.floor(BigGillespieTimeArray[j]) == i)[0]
            if len(index) > 0:
                index = index[0]

                BTGillespieTimeArray[j][i] = BigGillespieTimeArray[j][index]
                BTGillespieHealthyArray[j][i] = BigGillespieHealthyArray[j][index]
                BTGillespieEclipseArray[j][i] = BigGillespieEclipseArray[j][index]
                BTGillespieInfectedArray[j][i] = BigGillespieInfectedArray[j][i]
                BTGillespieDeadArray[j][i] = BigGillespieDeadArray[j][index]
                BTGillespieVirusArray[j][i] = BigGillespieVirusArray[j][index]
            else:
                break

    GillespieTimeArray = np.arange(0,LongestTime,1)
    GillespieHealthyArray = np.mean(BTGillespieHealthyArray,axis=0)
    GillespieEclipseArray = np.mean(BTGillespieEclipseArray,axis=0)
    GillespieInfectedArray = np.mean(BTGillespieInfectedArray,axis=0)
    GillespieDeadArray = np.mean(BTGillespieDeadArray,axis=0)
    GillespieVirusArray = np.mean(BTGillespieVirusArray,axis=0)

    index = np.argwhere(TimeArray==max(TimeArray))[0][0]

    plt.plot(GillespieTimeArray, GillespieHealthyArray/1001365, '.-', label="Gillespie")
    for i in range(NumberOfRuns):
        gtime = np. trim_zeros(BTGillespieTimeArray[i], trim='b')
        plt.plot(gtime, BTGillespieHealthyArray[i][:len(gtime)]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Cells Over time")
    plt.ylabel("Fraction of Uninfected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Gillespie_U.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(GillespieTimeArray, GillespieEclipseArray/1001365, '.-', label="Gillespie")
    for i in range(NumberOfRuns):
        plt.plot(BigGillespieTimeArray[i], BigGillespieEclipseArray[i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Eclipse Over time")
    plt.ylabel("Fraction of Eclipse Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Gillespie_E.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(GillespieTimeArray, GillespieInfectedArray/1001365, '.-', label="Gillespie")
    for i in range(NumberOfRuns):
        plt.plot(BigGillespieTimeArray[i], BigGillespieInfectedArray[i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Infected Over time")
    plt.ylabel("Fraction of Infected Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Gillespie_I.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(GillespieTimeArray, GillespieDeadArray/1001365, '.-', label="Gillespie")
    for i in range(NumberOfRuns):
        plt.plot(BigGillespieTimeArray[i], BigGillespieDeadArray[i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Dead Over time")
    plt.ylabel("Fraction of Dead of Cells")
    plt.xlabel("Days")
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Gillespie_D.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    plt.plot(GillespieTimeArray, GillespieVirusArray/1001365, '.-', label="Gillespie")
    for i in range(NumberOfRuns):
        plt.plot(BigGillespieTimeArray[i], BigGillespieVirusArray[i]/1001365, 'b-', alpha=0.2)
    plt.xticks(np.arange(0, index, 48),np.arange(0, index/24, 2))
#    plt.title("Virus Over time")
    plt.ylabel("Amount of Virus per Cell")
    plt.xlabel("Days")
    plt.yscale("log")
    plt.ylim( (10**-10,10**-2) )
    plt.xlim( (0,24*10) )
    plt.legend()
    plt.tick_params(axis='x', rotation=rotate)
    plt.savefig("Output_Plots/Gillespie_V.png", bbox_inches='tight', pad_inches=0.0)
    plt.close()

    gc.collect()
