// nvcc Gillespie_DV.cu -o program.out && ./program.out 

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <math.h>
#include <random>
#include <iostream>

// Model parameters
float Vo = 7.5*pow(10,-2);   //TCID_50/ml 
float beta = 3.2*pow(10,-5); //((TCID_50/ml)*d)^-1
float c = 5.2;               //d^-1
float p = 4.6*pow(10,-2);    //(TCID_50*ml)*d^-1
float k = 4.0;               //d^-1
float delta = 5.2;           //d^-1

float MOI = Vo/(4.0*pow(10,8));

//Initial Conditions 
int T0 = 1001365;
int E0 = 0;
int I0 = 0;
int D0 = 0;
float V0 = MOI*T0*pow(10,2);

int n_simulations = 10;
float endtime = 24*10; //maximum days of infection
int save = 1;
int num_saves = int(endtime/save)+1;
int* savetime = (int*) calloc(num_saves,sizeof(int));

char Path_to_Folder[100] = "";

float num = 5;
float scale = num/V0;

///////////////////////////////////////////////////////////////////
double PU(){
    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    return distribution(generator);
}

void convertVariables(){
    beta = beta*(4.0*pow(10,8))/T0; 

    beta = beta/24.0; 
    c = c/24.0;
    p = p/24.0;
    k = k/24.0;
    delta = delta/24.0;    
}

float Rate0;
float Rate1;
float Rate2;
float Rate3;
float Rate4;
float R;
float u;
float tau;

int NT;
int NE;
int NI;
int ND;
int NV;

float Time;
int T;
int E;
int I;
int D;
float V;

char Buffer[5];

int count;

float RandRateSum;


int main(void){

cout << beta[0] << endl;
cout << c << endl;
cout << rho << endl;
cout << k_paper/24.0 << endl;
cout << delta_paper/24.0 << endl;
cout << MOI[0] << endl;

convertVariables();

for (int i=0; i < num_saves; i++){
    savetime[i] = i;
}

for (int i = 0; i < n_simulations; i++){

    Time = 0.0;
    T = T0;
    E = E0;
    I = I0;
    D = D0;
    V = V0*scale;

    strcpy(Path_to_Folder, "");
    strcat(Path_to_Folder,"Data/Gillespie_Model/Deterministic/");
    strcat(Path_to_Folder,"Gillespie_");
    sprintf(Buffer,"%d",i);
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,".txt");

    char File1[100] = "";
    strcat(File1,Path_to_Folder);
    FILE *outfile1 = fopen(File1,"w");
    if (outfile1 == NULL){
        printf("Error opening file!\n");
        exit(0);
    }

    fclose(outfile1);

    printf("simulation = %d\n",i);

    count = 0;
    while (Time < endtime){
        Rate0 = beta*V/scale*T;
        Rate1 = k*E;
        Rate2 = delta*I;
        Rate3 = p*I*scale;
        Rate4 = c*V;

        R = Rate0 + Rate1 + Rate2;
        u = PU();
        tau = 1/R * log(1/u);

        NT = 0;
        NE = 0;
        NI = 0;
        ND = 0;
        NV = 0;

        u = PU();

        RandRateSum = u*R;

        if ((RandRateSum >= 0)&(RandRateSum < Rate0)){
            NT += -1;
            NE +=  1;
        }
        else if ((RandRateSum >= Rate0)&(RandRateSum < Rate0+Rate1)){
            NE += -1;
            NI +=  1;
        }
        else if ((RandRateSum >= Rate0+Rate1)&(RandRateSum < Rate0+Rate1+Rate2)){
            NI += -1;
            ND +=  1;
        }

        NT = MAX(NT, -T);
        NE = MAX(NE, -E);
        NI = MAX(NI, -I);

        Time += tau;
        T += NT;
        E += NE;
        I += NI;
        D += ND;
        V += tau*(Rate3 - Rate4);

        if (Time >= savetime[count]){

            FILE *outfile1 = fopen(File1,"a");
            if (outfile1 == NULL){
                printf("Error opening file!\n");
                exit(0);
            }

            fprintf(outfile1,"%.10f, %d, %d, %d, %d, %.10f", Time, T, E, I, D, V/scale);
            fprintf(outfile1,"\n");
            fclose(outfile1);
            count += 1;
        }
    }
}
}

