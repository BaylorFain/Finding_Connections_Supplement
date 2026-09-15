//For kernals
dim3 BlockConfig, GridConfig;

//For c
char Path_to_Folder[100] = "";
char Directroy[100] = "";
char** LocationData;
char* cells;
char* cells_GPU;
float* ecl;
float* ecl_GPU;
float* inf;
float* inf_GPU;
double* vtemp;
double* vtemp_GPU;
float* th;
float* th_GPU;
float* ut;
float* ut_GPU;
float* EclipsePhaseLength;
float* EclipsePhaseLength_GPU;
float* InfectionPhaseLength;
float* InfectionPhaseLength_GPU;
char* TransmissionInfection;
char* TransmissionInfection_GPU;
int NumberOfCells;

int NumberDead;
int NumberInfected;
int NumberEclipse;
int NumberHealthy;
double AmountOfVirus;
float Amountc2c;
float Amountcf;
float Amountboth;

double* seeds;
double* seeds_GPU;
