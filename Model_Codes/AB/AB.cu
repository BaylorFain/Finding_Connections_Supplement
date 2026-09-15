// nvcc AB.cu -o program.out && ./program.out

#define CODETESTINGCONDITIONS 0

//Simulation Parameters
int CELL2CELL = 0;
int FREECELL = 1;

//Change initial conditions
int INITIALCELL = 0;
int INITIALVIRUS = 1;

///////////////////////////////////////////////////////////////////////////////////
float timestep = 0.005;    //hours //0.005 hr = 18 sec, (1/3600) hr = 1 sec
float endtime = 24*10;   //hours
int Save = 1/(timestep); //the number of time the program saves to file, (1/timestep) results in 1 save every simulated hour; For shorter saves, must be wholw min, 1/2 = 30min, 1/4 = 15min, 1/20 = 3min
int NumberOfLayers = 607;//607; //607 is a million hexigon in a circle
int StartRuns = 0;
int NumberOfRuns = 10;

//Paper parameters
double Vo_paper = 7.5*pow(10,-2);   //TCID_50/ml 
float beta_paper = 3.2*pow(10,-5); //((TCID_50/ml)*d)^-1
float c_paper = 5.2;               //d^-1
float p_paper = 4.6*pow(10,-2);    //(TCID_50*ml)*d^-1
float k_paper = 4.0;               //d^-1
float delta_paper = 5.2;           //d^-1

double MOI = Vo_paper/(4.0*pow(10,8));

double T0 = 1001365;
double V0 = MOI*T0*pow(10,2);

//Physical Parameters////////////////////////////////////////////////
// ****SEE MAIN****  MOI; //pow(10,-5) to 1
// ****SEE MAIN**** double beta = 0.1;           //Infection rate, units: 1/(hour)       
double rho = p_paper/24.0;
// ****SEE MAIN**** double D = 2.16*pow(10,-8);  //Diffusion rate at 37 degrees celsius unit: m^2/h
double c = c_paper/24.0; //Clearance rate, units: 1/(hour) 
double deltx = 25.0*pow(10,-6);
double deltxprime = deltx*2;
double Dtsx2;// = D*timestep*pow(deltxprime,-2);

//Probability Constants////////////////////////////////////////////////////////////
float TauI = 1/(delta_paper/24.0); //Avg time for infection
float ni = 100.0;   //Number of compartments/partitions for time of infected cells
float TauE = 1/(k_paper/24.0); //Avg time for eclipse
float ne = 30.0;    //Number of compartments/partitions for time of eclipse cells
// ****SEE MAIN**** float probi = 0.2;  //Probability of cell to cell infection

#include "Functions/Includes.h"
#include "Functions/Variables.h"
#include "Functions/Functions.h"

int main(void){

beta_paper = beta_paper*(4.0*pow(10,8)); 
double MOI[1] = {V0/T0};
double beta[1] = {beta_paper/24.0};
double D[1] = {2.16*pow(10,-8)};
double probi[1] = {0.0};

cout << beta[0] << endl;
cout << c << endl;
cout << rho << endl;
cout << k_paper/24.0 << endl;
cout << delta_paper/24.0 << endl;
cout << MOI[0] << endl;

/////////////////////////////////////////////////////////////////////////
for(int q=0;q<(sizeof(MOI)/sizeof(MOI[0]));q++){
for(int k=0;k<(sizeof(probi)/sizeof(probi[0]));k++){
for(int m=0;m<(sizeof(D)/sizeof(D[0]));m++){
for(int n=0;n<(sizeof(beta)/sizeof(beta[0]));n++){
    Dtsx2 = D[m]*timestep*pow(deltxprime,-2);
    checksAndClear(D[m], timestep, deltxprime);
/////////////////////////////////////////////////////////////////////////
    //Loop For The number Of Simulations To Run Per Setting/////////////
    for(int BigIndex=0;BigIndex<NumberOfRuns;BigIndex++){

        //Creating Save Path///////////////////////////////////////////
        creatingPathToFolderAndDirectory(StartRuns+BigIndex, NumberOfLayers, MOI[q], probi[k], beta[n], c, rho, D[m]);
        //Creating placeholder variables for multipy runs/////////////
        int cell2cell = CELL2CELL;
        int freecell = FREECELL;

        //Building Cells/////////////////////////////////////////////
        creatingCellLocations();
        
        //Number of initial infected cells///////////////////////////
        int Ni = NumberOfCells*MOI[q]; 
        int Nx = (2*NumberOfLayers-1);      //Range of cells on x axis
        int Ny = (2*NumberOfLayers-1);      //Range of cells on y axis
        
        //Makeing empty matrices////////////////////////////////////
        allocateMemory(Nx, Ny);
        
        //Initial Conditions///////////////////////////////////////
        initailConditions(Nx, Ny, MOI[q]);

        //Initializing print files/////////////////////////////////
        printToFileCellAndVirusAnalysisInitial(Nx, Ny);

        //Infects a random cell////////////////////////////////////
        if (INITIALCELL == 1){
            if(Ni < 1){ printf("Use larger MOI"); exit(0);} 
            infectANumberOfCellsRandomly(Nx, Ny, Ni);
        }
        
        //Load constants to pass to GPU///////////////////////////
        loadConstants(MOI[q], probi[k], beta[n], D[m]);
        
        //Allocate memory on GPU/////////////////////////////////
        deviceSetupAndMemoryAllocation(Nx, Ny);

        //Copy memory to GPU/////////////////////////////////////
        memoryCopyToDevice(Nx, Ny);
        
        /////////////////////////////////////////////////////////
        curandState* state;
        cudaMalloc((void**)&state, Nx*Ny*sizeof(curandState));
        errorCheck("cudaMalloc Random Setup");
        cuRand_Setup<<<GridConfig,BlockConfig>>>(NumberOfLayers, seeds_GPU, state);
        errorCheck("Random Setup");
        ////////////////////////////////////////////////////////
        
        //Time step and counting///////////////////////////////
        int NumberofTimeSteps = endtime/timestep;
        int NumberofSavedTimeSteps = NumberofTimeSteps/Save;
        float timestepcount = 0.0; //equal to the number of timesteps elapsed
        while(timestepcount < (NumberofTimeSteps-1)){

            //Calls GPU code////////////////////////////////////
            kernel<<<GridConfig,BlockConfig>>>(timestepcount, cells_GPU, vtemp_GPU, ut_GPU, ecl_GPU, inf_GPU, th_GPU, EclipsePhaseLength_GPU, InfectionPhaseLength_GPU, TransmissionInfection_GPU, SystemConstants, cell2cell, freecell, NumberOfLayers, probi[k], seeds_GPU, state);

            //These actions are preformed when a measurement wants to be preformed. The time frame is dictated by the 'Save' variable////
            if((int(timestepcount)%Save) == 0){ 
                //Copys memory back to CPU/////////////////////
                cudaMemcpy( cells, cells_GPU, Nx*Ny*2*sizeof(char), cudaMemcpyDeviceToHost );
                errorCheck("cudaMemcpy cells DtoH");
                cudaMemcpy( vtemp, vtemp_GPU, Nx*Ny*2*sizeof(double), cudaMemcpyDeviceToHost );
                errorCheck("cudaMemcpy vtemp DtoH");
                //Prints last save analitics to file//Maybe move below analysis dish/////////////////////////////////////////////////////////

                //Analysis the dish//////////////////////////////////
                NumberDead = 0;
                NumberInfected = 0;
                NumberEclipse = 0;
                NumberHealthy = 0;
                AmountOfVirus = 0.0;
                for(int j=0; j<Ny; j++){
                    for(int i=0; i<Nx; i++){
                        AmountOfVirus = AmountOfVirus + vtemp[i+Nx*j+Nx*Ny*0];
                        
                        if(cells[i+Nx*j+Nx*Ny*0] == 'd'){
                            NumberDead = NumberDead + 1;
                        }
                        else if(cells[i+Nx*j+Nx*Ny*0] == 'i'){
                            NumberInfected = NumberInfected + 1;
                        }
                        else if(cells[i+Nx*j+Nx*Ny*0] == 'e'){
                            NumberEclipse = NumberEclipse +1;
                        }
                        else if(cells[i+Nx*j+Nx*Ny*0] == 'h'){
                            NumberHealthy = NumberHealthy + 1;
                        }
                    }
                }
                printToFileCellAndVirusAnalysis(timestepcount*timestep+1.0);
            }
            timestepcount = timestepcount+1.0;
        }

        //Copies and prints number of cell-to-cell and cell-free infections to file
        cudaMemcpy( TransmissionInfection, TransmissionInfection_GPU, Nx*Ny*2*sizeof(char), cudaMemcpyDeviceToHost );
        errorCheck("cudaMemcpy TransmissionInfection DtoH");
        Amountc2c = 0.0;
        Amountcf = 0.0;
        Amountboth = 0.0;
        for(int j=0; j<Ny; j++){
            for(int i=0; i<Nx; i++){
                if(TransmissionInfection[i+Nx*j+Nx*Ny*0] == 'c'){
                    Amountc2c = Amountc2c + 1.0;
                }
                else if(TransmissionInfection[i+Nx*j+Nx*Ny*0] == 'f'){
                    Amountcf = Amountcf + 1.0;
                }
                else if(TransmissionInfection[i+Nx*j+Nx*Ny*0] == 'b'){
                    Amountboth = Amountboth + 1.0;
                }
            }
        }
        char Fileti[100] = "";
        strcat(Fileti,Path_to_Folder);
        strcat(Fileti,"/TransmissionInfection.txt");
        FILE *outfileti = fopen(Fileti,"w");
        if (outfileti == NULL){
            printf("Error opening file!\n");
            exit(0);
        }
        fprintf(outfileti, "Amount via c2c = %f\n", Amountc2c);
        fprintf(outfileti, "Amount via cf = %f\n", Amountcf);
        fprintf(outfileti, "Amount via both = %f\n", Amountboth);
        fclose(outfileti);

        //Writes a file with all of our parameters/variables/////////
        createParameterFile(timestep, endtime, NumberofSavedTimeSteps, timestepcount, AmountOfVirus, MOI[q], beta[n], rho, D[m], c, deltxprime, TauI, ni, TauE, ne, probi[k]);

        //Frees memory//////////////////////////////////////////////
        freeMemory();
        cudaFree(state);
    }
}
}
}
}
printf("\nPROGRAM DONE\n");
}
