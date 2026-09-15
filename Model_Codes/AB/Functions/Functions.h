using namespace std;

float Te(float TauE, float ne){
//    Picks a random number from the gamma distribution
//    The number is to be used as a time step in the Eclipse Time Matrix
    random_device rd;
    default_random_engine generator(rd());
    gamma_distribution<double> distribution(ne, TauE/ne);

    return distribution(generator);
}

float Ti(float TauI, float ni){
//    Picks a random number from the gamma distribution
//    The number is to be used as a time step in the Infected Time Matrix
    random_device rd;
    default_random_engine generator(rd());
    gamma_distribution<double> distribution(ni, TauI/ni);

    return distribution(generator);
}

double PU1(){
//    Picks a random number from a uniform distribution
    random_device rd;
    default_random_engine generator(rd());
    uniform_real_distribution<double> distribution(0.0,1.0);

    return distribution(generator);
}

void checksAndClear(double D, double timestep, double deltxprime){
//Checks for Heisenberg status of viral diffusion
if(D*timestep/pow(deltxprime,2.0) > 0.25){
    printf("%.1f",D*timestep/pow(deltxprime,2.0));
    printf("CHANGE PARAMETERS TO FIT DIFFUSION LIMITS. VALUE MUST BE UNDER 0.5. VALUE SHOWN ABOVE");
    exit(0);
}
//Clear Terminal
system("clear");
}

void creatingPathToFolderAndDirectory(int BigIndex, int NumberOfLayers, double MOI, double probi, double beta, double c, double p, double D){
    char TransmissionType[10] = "";
    if (CELL2CELL == 1){
		if (FREECELL == 1){
            strcat(TransmissionType,"Both");
        }
        else {
            strcat(TransmissionType,"CELL2CELL");
        }
    }
    else if(CELL2CELL == 0){
	    if (FREECELL == 0){
            strcat(TransmissionType,"Neither");
        }
        else{
            strcat(TransmissionType,"CELLFREE");
        }
    }
    
    char Buffer[5]; //Buffer String For Conversion To Char
    time_t RawTime = time(NULL);
    tm* SpecificMoment = localtime(&RawTime);

    strcpy(Path_to_Folder, "");
    strcpy(Directroy, "");
    
/*    strcat(Path_to_Folder,"ViralModel/");*/
    strcat(Path_to_Folder,"../Data/AB_Model/");

    sprintf(Buffer,"%d",NumberOfLayers);
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,"_");
    sprintf(Buffer,"%d",BigIndex);
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,"-");
    strcat(Path_to_Folder,TransmissionType);
    strcat(Path_to_Folder,"_");

    sprintf(Buffer,"%.1f",MOI);
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,"-");
    strcat(Path_to_Folder,"MOI");
    strcat(Path_to_Folder,"_");
  
    sprintf(Buffer,"%.1f",beta);
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,"-");
    strcat(Path_to_Folder,"beta");
    strcat(Path_to_Folder,"_");
    
    sprintf(Buffer,"%.1f",c);
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,"-");
    strcat(Path_to_Folder,"c");
    strcat(Path_to_Folder,"_");
    
    sprintf(Buffer,"%.1f",rho);
    strcat(Path_to_Folder,Buffer);
    strcat(Path_to_Folder,"-");
    strcat(Path_to_Folder,"p");
    strcat(Path_to_Folder,"_");

    strcat(Directroy,"mkdir -p ");
    strcat(Directroy,Path_to_Folder);
    int check = system(strdup(Directroy));
    if(check != 0){
        exit(0);
    }
}

void creatingCellLocations(){
    float SideLenght = (2.0/3.0);
    int RadiusScale = 0;
    for(int i=0; i<NumberOfLayers; i++){
        if(i == 0){
            RadiusScale = RadiusScale + 1;
        }
        else{
            if((i)%2 == 1){
                RadiusScale = RadiusScale + 1;
            }
            else{
                RadiusScale = RadiusScale + 2;
            }
        }
    }
    float RadiusOfCircle = SideLenght*RadiusScale;

    int count = 0;
    for(int i=0; i<NumberOfLayers; i++){
        count = count + i;
    }
    int NumberOfHexagons=(count)*6+1;

    float** coord;
    int n = NumberOfHexagons;
    int m = 3;
    coord = (float**) calloc(n,sizeof(float*));  
    for (int i = 0; i < n; i++){
       coord[i] = (float*) calloc(m,sizeof(float));
    }

    float** percyclecoord;
    n = NumberOfHexagons;
    m = 3;
    percyclecoord = (float**) calloc(n,sizeof(float*));  
    for (int i = 0; i < n; i++){
       percyclecoord[i] = (float*) calloc(m,sizeof(float));  
    }

    int temp;
    for(int j=0; j<NumberOfLayers; j++){
        for(int i=0; i<(2*j); i++){
            if(i < j){
                temp = i;
            }
            percyclecoord[i+(j-1)*j+1][0] =  -temp-1;
            percyclecoord[i+(j-1)*j+1][1] =   temp+j-i;
            percyclecoord[i+(j-1)*j+1][2] =  -j+1+i;
            
        }
    }
    float c0[3] = {percyclecoord[0][0], percyclecoord[0][1], percyclecoord[0][2]};
    coord[0][2] = c0[2];
    coord[0][1] = c0[1];
    coord[0][0] = c0[0];

    count = 0;
    for(int j=0; j<(NumberOfHexagons/3); j++){
        for(int i=0; i<3; i++){
            coord[(i+0)%3+3*j+1][2] = percyclecoord[j+1][i]+c0[i];
            coord[(i+1)%3+3*j+1][1] = percyclecoord[j+1][i]+c0[i];
            coord[(i+2)%3+3*j+1][0] = percyclecoord[j+1][i]+c0[i];
        }
    }

    float hi = coord[0][0];
    float vi = coord[0][2];
    float xmin = INFINITY;
    float xcoord;
    float ycoord;
    double dist;
    for(int i=0; i<NumberOfHexagons; i++){
        xcoord = coord[i][0];
        if(coord[i][0] < xmin){
            xmin = coord[i][0];
        }
        ycoord = (2.0*sin(PI*(60.0/180.0))*(coord[i][1]-coord[i][2])/3.0)+vi;
        dist = sqrtf(pow(double(xcoord-hi),2.0)+pow(double(ycoord-vi),2.0));
        if(dist >= RadiusOfCircle){
            coord[i][0] = 5000.0;
            coord[i][1] = 0.0;
            coord[i][2] = 0.0;
        }
    }

    n = ((2*NumberOfLayers)-1);
    m = ((2*NumberOfLayers)-1);
    LocationData = (char**) malloc(n*sizeof(char*));  
    for(int j=0; j<n; j++){ 
        LocationData[j] = (char*) malloc(m*sizeof(char));  
        for(int i=0; i<m; i++){
            LocationData[j][i] = 'o';
       }
    }
    
    NumberOfCells = 0;
    for(int i=0; i<NumberOfHexagons; i++){
        if(coord[i][0] != 5000.0){
            LocationData[int(coord[i][2])-int(xmin)][int(coord[i][0])-int(xmin)] = 'h';
            NumberOfCells = NumberOfCells + 1;
        }
    }
    
    char File2[100] = "";
    strcat(File2,Path_to_Folder);
    strcat(File2,"/Parameters.txt");
    FILE *outfile2 = fopen(File2,"w");
    if (outfile2 == NULL){
        printf("Error opening file!\n");
        exit(0);
    }
    fprintf(outfile2, "Hexagon Side Lenght: %f\n", SideLenght);
    fprintf(outfile2, "Number of Layers: %d\n", NumberOfLayers);
    fprintf(outfile2, "Radius of Circle: %f\n", RadiusOfCircle);
    fprintf(outfile2, "Number of Cells: %d\n", NumberOfCells);
    fclose(outfile2);
    
    for (int i = 0; i < NumberOfHexagons; i++){  
       free(coord[i]);  
    }     
    free(coord);
    
    for (int i = 0; i < NumberOfHexagons; i++){  
       free(percyclecoord[i]);  
    }     
    free(percyclecoord);
}

void allocateMemory(int Nx, int Ny){    
    //Produces a matrix for the cells
    cells = (char*) malloc(Nx*Ny*2*sizeof(char));
    //Produces a matrix that will track the amount virus above each cell
    vtemp = (double*) malloc(Nx*Ny*2*sizeof(double));
    //Produces a univeral time matrix (ut)
    ut = (float*) malloc(Nx*Ny*sizeof(float));
    
    //Produces a time matrix for after eclipse phase (e)
    ecl = (float*) malloc(Nx*Ny*sizeof(float));
    //Produces a time matrix for after infection phase (i)
    inf = (float*) malloc(Nx*Ny*sizeof(float));
    //Produces a time matrix hor healthy cells (t)
    th = (float*) malloc(Nx*Ny*sizeof(float));
    
    //Produces an array of eclipse phase durations for cells
    EclipsePhaseLength = (float*) malloc(Nx*Ny*sizeof(float));
    //Produces an array of infection phase durations for cells
    InfectionPhaseLength = (float*) malloc(Nx*Ny*sizeof(float));

    //TransmissionInfection
    TransmissionInfection = (char*) malloc(Nx*Ny*2*sizeof(char));
    
    seeds = (double*) malloc(Nx*Ny*sizeof(double));
    
}

void initailConditions(int Nx, int Ny, double MOI){
    for(int j=0; j<Ny; j++){
        for(int i=0; i<Nx; i++){
            for(int k=0;k<2;k++){
                cells[i+Nx*j+Nx*Ny*k] = LocationData[i][j];
                vtemp[i+Nx*j+Nx*Ny*k] = 0.0;
                TransmissionInfection[i+Nx*j+Nx*Ny*k]  = 'o';
            }
            ut[i+Nx*j] = 0.0;
            ecl[i+Nx*j] = 0.0;
            inf[i+Nx*j] = 0.0;
            th[i+Nx*j] = 0.0;
            EclipsePhaseLength[i+Nx*j] = Te(TauE,ne);
            InfectionPhaseLength[i+Nx*j]  = Ti(TauI,ni);
            seeds[i+Nx*j] = PU1()*pow(10,12)+PU1()*pow(10,6);
       }
    }

    if (INITIALVIRUS == 1){
        for(int j=0; j<Ny; j++){
            for(int i=0; i<Nx; i++){
                for(int k=0;k<2;k++){
                    if(cells[i+Nx*j+Nx*Ny*k] != 'o'){
                        vtemp[i+Nx*j+Nx*Ny*k] = MOI;
                    }
                }
           }
        }
    }
}

void infectANumberOfCellsRandomly(int Nx, int Ny, int Ni){
    if(CODETESTINGCONDITIONS == 1){
        cells[(NumberOfLayers-1)+Nx*(NumberOfLayers-1)+Nx*Ny*0] = 'i';
        cells[(NumberOfLayers-1)+Nx*(NumberOfLayers-1)+Nx*Ny*1] = 'i'; //Only the center cell
    }
    else {
        srand(time(NULL));
        int randx;
        int randy; 
        int NumberOfInfectedCellsCount = 0;
        while(NumberOfInfectedCellsCount < Ni){
            randx = (rand()%Nx);
            randy = (rand()%Ny);
            if((cells[randx+Nx*randy+Nx*Ny*0] != 'o') && (cells[randx+Nx*randy+Nx*Ny*0] == 'h')){
                cells[randx+Nx*randy+Nx*Ny*0] = 'e';
                cells[randx+Nx*randy+Nx*Ny*1] = 'e';
                ecl[randx+Nx*randy] = Te(TauE,ne);
                NumberOfInfectedCellsCount = NumberOfInfectedCellsCount + 1;
            }
        }
    }
}

void printToFileCellAndVirusInitial(int Nx, int Ny, int NumberOfLayers){
    char File3[100] = "";
    strcat(File3,Path_to_Folder);
    strcat(File3,"/cells_over_time.txt");
    FILE *outfile3 = fopen(File3,"w");
    if (outfile3 == NULL){
        printf("Error opening file!\n");
        exit(0);
    }
    for(int i=0; i<((2*NumberOfLayers)-1); i++){
        for(int j=0; j<((2*NumberOfLayers)-1); j++){
            fprintf(outfile3,"%c,",LocationData[i][j]);
        }
            fprintf(outfile3,"\n");
    }
    fclose(outfile3);

    char File4[100] = "";
    strcat(File4,Path_to_Folder);
    strcat(File4,"/virus_over_time.txt");
    FILE *outfile4 = fopen(File4,"w");
    if (outfile4 == NULL){
        printf("Error opening file!\n");
        exit(0);
    }
    for(int i=0; i<((2*NumberOfLayers)-1); i++){
        for(int j=0; j<((2*NumberOfLayers)-1); j++){
            fprintf(outfile4,"%f,",0.0);
        }
            fprintf(outfile4,"\n");
    }
    fclose(outfile4);
}

void printToFileCellAndVirusAnalysisInitial(int Nx, int Ny){
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
    
    char File9[100] = "";
    strcat(File9,Path_to_Folder);
    strcat(File9,"/PerTimeStep.txt");
    FILE *outfile9 = fopen(File9,"w");
    if (outfile9 == NULL){
        printf("Error opening file!\n");
        exit(0);
    }
    
    fprintf(outfile9,"%.2f, %d, %d, %d, %d, %f,", 0.0, NumberHealthy, NumberEclipse, NumberInfected, NumberDead, AmountOfVirus);
    fprintf(outfile9,"\n");

    fclose(outfile9);
}

void printToFileCellAndVirus(int Nx, int Ny, int NumberOfLayers){
    char File5[100] = "";
    strcat(File5,Path_to_Folder);
    strcat(File5,"/cells_over_time.txt");
    FILE *outfile5 = fopen(File5,"a");
    if (outfile5 == NULL){
        printf("Error opening file!\n");
        exit(0);
    }
    for(int i=0; i<((2*NumberOfLayers)-1); i++){
        for(int j=0; j<((2*NumberOfLayers)-1); j++){
            fprintf(outfile5,"%c,",cells[i+Nx*j+Nx*Ny*0]);
        }
            fprintf(outfile5,"\n");
    }
    fclose(outfile5);
   
    char File6[100] = "";
    strcat(File6,Path_to_Folder);
    strcat(File6,"/virus_over_time.txt");
    FILE *outfile6 = fopen(File6,"a");
    if (outfile6 == NULL){
        printf("Error opening file!\n");
        exit(0);
    }
    for(int i=0; i<((2*NumberOfLayers)-1); i++){
        for(int j=0; j<((2*NumberOfLayers)-1); j++){
            fprintf(outfile6,"%f,",vtemp[i+Nx*j+Nx*Ny*1]);
        }
            fprintf(outfile6,"\n");
    }
    fclose(outfile6);
}

void printToFileCellAndVirusAnalysis(float timestep){
    char File8[100] = "";
    strcat(File8,Path_to_Folder);
    strcat(File8,"/PerTimeStep.txt");
    FILE *outfile8 = fopen(File8,"a");
    if (outfile8 == NULL){
        printf("Error opening file!\n");
        exit(0);
    }
    
    fprintf(outfile8,"%.2f, %d, %d, %d, %d, %f,", timestep, NumberHealthy, NumberEclipse, NumberInfected, NumberDead, AmountOfVirus);
    fprintf(outfile8,"\n");

    fclose(outfile8);
}

void createParameterFile(float timestep, float endtime, int NumberofSavedTimeSteps, float timestepcount, double AmountOfVirus, double MOI, double beta, double rho, double D, double c, double deltxprime, float TauI, float ni, float TauE, float ne, double probi){
    char File7[100] = "";
    strcat(File7,Path_to_Folder);
    strcat(File7,"/Parameters.txt");
    FILE *outfile7 = fopen(File7,"a");
    if (outfile7 == NULL){
        printf("Error opening file!\n");
        exit(0);
    }
    fprintf(outfile7, "Time Step: %f\n", timestep);
    fprintf(outfile7, "Initial End Time: %f\n", endtime);
    fprintf(outfile7, "Number of Saved Time Steps: %d\n", NumberofSavedTimeSteps);
    fprintf(outfile7, "Actual Hours Simulated: %f\n", timestepcount*timestep);
    fprintf(outfile7, "Final Amount of Virus: %f\n", AmountOfVirus);
    fprintf(outfile7, "MOI: %f\n", MOI);    
    fprintf(outfile7, "beta: %f\n", beta);
    fprintf(outfile7, "log10(beta): %f\n", log10(beta));
    fprintf(outfile7, "p: %f\n", rho);
    fprintf(outfile7, "D: %f\n", D);
    fprintf(outfile7, "log10(D): %f\n", log10(D));
    fprintf(outfile7, "c: %f\n", c);
    fprintf(outfile7, "delta x: %f\n", deltxprime);
    fprintf(outfile7, "TauI: %f\n", TauI);
    fprintf(outfile7, "ni: %f\n", ni);
    fprintf(outfile7, "TauE: %f\n", TauE);
    fprintf(outfile7, "ne: %f\n", ne);
    fprintf(outfile7, "Probability of cell to cell infection: %f\n", probi);

    fclose(outfile7);
}

void freeMemory(){ 
    for(int i=0; i<((2*NumberOfLayers)-1); i++){
        free(LocationData[i]);
    }
    free(LocationData);  
    free(cells);   
    free(ecl);
    free(inf);  
    free(vtemp);
    free(th);
    free(ut);
    free(EclipsePhaseLength);
    free(InfectionPhaseLength);
    free(TransmissionInfection);
    free(seeds);
    
    cudaFree(cells_GPU);
    cudaFree(ecl_GPU);
    cudaFree(inf_GPU);
    cudaFree(vtemp_GPU);
    cudaFree(th_GPU);
    cudaFree(ut_GPU);
    cudaFree(EclipsePhaseLength_GPU);
    cudaFree(InfectionPhaseLength_GPU);
    cudaFree(TransmissionInfection_GPU);
    cudaFree(seeds_GPU);

}

void errorCheck(const char *message){
  cudaError_t  error;
  error = cudaGetLastError();

  if(error != cudaSuccess)
  {
    printf("\n CUDA ERROR: %s = %s\n", message, cudaGetErrorString(error));
    exit(0);
  }
}

struct systemConstantsStruct
{
    double MOI;
    double beta;
    double rho;
    double D;
    double c;
    double deltx;
    double deltxprime;
    double Dtsx2;

    float TauI;
    float TauE;
    float ne;
    float ni;
    double probi;
    
    float timestep;
};

systemConstantsStruct SystemConstants;

void loadConstants(double MOI, double probi, double beta, double D){
    SystemConstants.MOI = MOI;
    SystemConstants.beta = beta;
    SystemConstants.rho = rho;
    SystemConstants.D = D;
    SystemConstants.c = c;
    SystemConstants.deltx = deltx;
    SystemConstants.deltxprime = deltxprime;
    SystemConstants.Dtsx2 = Dtsx2;

    SystemConstants.TauI = TauI;
    SystemConstants.TauE = TauE;
    SystemConstants.ne = ne;
    SystemConstants.ni = ni;
    SystemConstants.probi = probi;
    
    SystemConstants.timestep = timestep;
}

void deviceSetupAndMemoryAllocation(int Nx, int Ny){

	BlockConfig.x = 16;
	BlockConfig.y = 16;
	BlockConfig.z = 1;
	
	GridConfig.x = (Nx-1)/BlockConfig.x + 1;
	GridConfig.y = (Ny-1)/BlockConfig.y + 1;
	GridConfig.z = 1;
	
	cudaMalloc((void**)&cells_GPU, Nx*Ny*2*sizeof(char));
	errorCheck("cudaMalloc cells Mem");
	cudaMalloc((void**)&vtemp_GPU, Nx*Ny*2*sizeof(double));
	errorCheck("cudaMalloc vtemp Mem");
	cudaMalloc((void**)&ut_GPU, Nx*Ny*sizeof(float));
	errorCheck("cudaMalloc ut Mem");
	
	cudaMalloc((void**)&ecl_GPU, Nx*Ny*sizeof(float));
	errorCheck("cudaMalloc ecl Mem");
	cudaMalloc((void**)&inf_GPU, Nx*Ny*sizeof(float));
	errorCheck("cudaMalloc inf Mem");
	cudaMalloc((void**)&th_GPU, Nx*Ny*sizeof(float));
	errorCheck("cudaMalloc th Mem");
	
	cudaMalloc((void**)&EclipsePhaseLength_GPU, Nx*Ny*sizeof(float));
	errorCheck("cudaMalloc EclipsePhaseLength Mem");
	cudaMalloc((void**)&InfectionPhaseLength_GPU, Nx*Ny*sizeof(float));
	errorCheck("cudaMalloc InfectionPhaseLength Mem");

	cudaMalloc((void**)&TransmissionInfection_GPU, Nx*Ny*2*sizeof(char));
	errorCheck("cudaMalloc TransmissionInfection Mem");

	cudaMalloc((void**)&seeds_GPU, Nx*Ny*sizeof(double));
	errorCheck("cudaMalloc seeds Mem");

}

void memoryCopyToDevice(int Nx, int Ny){
    cudaMemcpy( cells_GPU, cells, Nx*Ny*2*sizeof(char), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy cells HtoD");
    cudaMemcpy( vtemp_GPU, vtemp, Nx*Ny*2*sizeof(double), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy vtemp HtoD");
    cudaMemcpy( ut_GPU, ut, Nx*Ny*sizeof(float), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy ut HtoD");

    cudaMemcpy( ecl_GPU, ecl, Nx*Ny*sizeof(float), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy ecl HtoD");
    cudaMemcpy( inf_GPU, inf, Nx*Ny*sizeof(float), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy inf HtoD");
    cudaMemcpy( th_GPU, th, Nx*Ny*sizeof(float), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy th HtoD");

    cudaMemcpy( EclipsePhaseLength_GPU, EclipsePhaseLength, Nx*Ny*sizeof(float), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy EclipsePhaseLength HtoD");
    cudaMemcpy( InfectionPhaseLength_GPU, InfectionPhaseLength, Nx*Ny*sizeof(float), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy InfectionPhaseLength HtoD");

    cudaMemcpy( TransmissionInfection_GPU, TransmissionInfection, Nx*Ny*2*sizeof(char), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy TransmissionInfection HtoD");
    
    cudaMemcpy( seeds_GPU, seeds, Nx*Ny*sizeof(double), cudaMemcpyHostToDevice );
    errorCheck("cudaMemcpy seeds HtoD");
    
}

__global__ void cuRand_Setup(int NumberOfLayers, double *seeds, curandState *state){
    int Row = threadIdx.x + blockIdx.x * blockDim.x;
    int Column =  threadIdx.y + blockIdx.y * blockDim.y;

    int NX = (2*NumberOfLayers-1);
    
    if(Row+NX*Column < NX*NX){
        long long seed = seeds[Row+NX*Column];
        curand_init (seed, Row+NX*Column, 0, &state[Row+NX*Column]);
    }


}

__device__ double PU_GPU(curandState *state){
//    Picks a random number from a uniform distribution

    double Random = curand_uniform(state);

    return Random;
}

__global__ void kernel(float timestepcount, char *cells, double *vtemp, float *ut, float *ecl, float *inf, float *th,  float *epl, float *ipl, char *ti, systemConstantsStruct constant, int cell2cell, int freecell, int NumberOfLayers, double probi, double *seeds, curandState *state){

    int Row = threadIdx.x + blockIdx.x * blockDim.x;
    int Column =  threadIdx.y + blockIdx.y * blockDim.y;

    int NX = (2*NumberOfLayers-1);
    int NY = (2*NumberOfLayers-1);
    int NXNY = NX*NY;

    if((cells[Row+NX*Column+NXNY*0] != 'o') && (Row+NX*Column < NXNY)){
        //Virus Spreads
        int AboveRow = Row-1;   //row coordinate above cell
        int LeftColumn = Column-1;   //column coordinate left of cell
        int BelowRow = Row+1;   //row coordinate below cell
        int RightColumn = Column+1;   //column coordinate right of cell

//          if the cell one row up doesn't exist, it's taken out of the equation
        if(AboveRow < 0){
            AboveRow = Row;
        }
//          if the cell one column to the left doesn't exist, it's taken out of the equation
        if(LeftColumn < 0){
            LeftColumn = Column;
        }
//          if the cell one row down doesn't exist, it's taken out of the equation
        if(BelowRow > (NY-1)){
            BelowRow = Row;
        }
//          if the cell one column to the right doesn't exist, it's taken out of the equation
        if(RightColumn > (NX-1)){
            RightColumn = Column;
        }

        if(cells[AboveRow+NX*Column+NXNY*0] == 'o'){
            AboveRow = Row;
        }
        if(cells[AboveRow+NX*RightColumn+NXNY*0] == 'o'){
            AboveRow = Row;
            RightColumn = Column;
        }
        if(cells[Row+NX*RightColumn+NXNY*0] == 'o'){
            RightColumn = Column;
        }
        if(cells[BelowRow+NX*Column+NXNY*0] == 'o'){
            BelowRow = Row;
        }
        if(cells[Row+NX*LeftColumn+NXNY*0] == 'o'){
            LeftColumn = Column;
        }
        if(cells[BelowRow+NX*LeftColumn+NXNY*0] == 'o'){
            BelowRow = Row;
            LeftColumn = Column;
        }
        
        double NNN = (vtemp[AboveRow+NX*Column+NXNY*0] + vtemp[AboveRow+NX*RightColumn+NXNY*0] + vtemp[Row+NX*RightColumn+NXNY*0] + vtemp[BelowRow+NX*Column+NXNY*0] + vtemp[Row+NX*LeftColumn+NXNY*0] + vtemp[BelowRow+NX*LeftColumn+NXNY*0]);
        
        double p = 0.0;
        if(cells[Row+NX*Column+NXNY*0] == 'i'){
            p = constant.rho;
        }
        double VirusProduced = p*constant.timestep;
        double VirusDecay = constant.c*vtemp[Row+NX*Column+NXNY*0]*constant.timestep;
        double VirusOut = 4.0*constant.Dtsx2*vtemp[Row+NX*Column+NXNY*0];
        double VirusIn = 2.0*constant.Dtsx2*NNN/3.0;

        __syncthreads();
        
        vtemp[Row+NX*Column+NXNY*1] = vtemp[Row+NX*Column+NXNY*0] + VirusProduced - VirusOut + VirusIn - VirusDecay;
        if(vtemp[Row+NX*Column+NXNY*1] < pow(10.0,-10.0)){
            vtemp[Row+NX*Column+NXNY*1] = 0.0;
        }
        
        //The Cell behavior
        if(cells[Row+NX*Column+NXNY*0] == 'i'){
            // Infectied
            if(ut[Row+NX*Column] > (inf[Row+NX*Column] + ecl[Row+NX*Column] + th[Row+NX*Column])){
                cells[Row+NX*Column+NXNY*1] = 'd';
                if(CODETESTINGCONDITIONS == 1){
                    cells[Row+NX*Column+NXNY*1] = 'i';
                }
            }
        }
        else if(cells[Row+NX*Column+NXNY*0] == 'e'){
            // Eclipse
            if(ut[Row+NX*Column] > (ecl[Row+NX*Column] + th[Row+NX*Column])){
                cells[Row+NX*Column+NXNY*1] = 'i';
                inf[Row+NX*Column] = inf[Row+NX*Column] + ipl[Row+NX*Column];
            }
        }
        else if(cells[Row+NX*Column+NXNY*0] == 'h'){
            // Healthy
            th[Row+NX*Column] = th[Row+NX*Column] + constant.timestep;
            
            if(cell2cell == 1){  
                // Cell to cell transmission
                int AboveRow = Row-1;   //row coordinate above cell
                int LeftColumn = Column-1;   //column coordinate left of cell
                int BelowRow = Row+1;   //row coordinate below cell
                int RightColumn = Column+1;   //column coordinate right of cell
                
        //        if the cell one row up doesn't exist, it's taken out of the equation
                if(AboveRow < 0){         
                    AboveRow = 0;
                }
        //        if the cell one column to the left doesn't exist, it's taken out of the equation
                if(LeftColumn < 0){
                    LeftColumn = 0;
                }
        //        if the cell one row down doesn't exist, it's taken out of the equation
                if(BelowRow > NY-1){
                    BelowRow = 0;
                }
        //        if the cell one column to the right doesn't exist, it's taken out of the equation
                if(RightColumn > NX-1){
                    RightColumn = 0;
                }

                if(cells[Row+NX*LeftColumn+NXNY*0] == 'i'){
                    if(PU_GPU(&state[Row+NX*Column]) < constant.probi){
                        cells[Row+NX*Column+NXNY*1] = 'e';
                        ecl[Row+NX*Column] = epl[Row+NX*Column];
                        ti[Row+NX*Column+NXNY*0] = 'c';
                    }
                }
                if(cells[Row+NX*RightColumn+NXNY*0] == 'i'){
                    if(PU_GPU(&state[Row+NX*Column]) < constant.probi){
                        cells[Row+NX*Column+NXNY*1] = 'e';
                        ecl[Row+NX*Column] = epl[Row+NX*Column];
                        ti[Row+NX*Column+NXNY*0] = 'c';
                    }
                }
                if(cells[AboveRow+NX*Column+NXNY*0] == 'i'){
                    if(PU_GPU(&state[Row+NX*Column]) < constant.probi){
                        cells[Row+NX*Column+NXNY*1] = 'e';
                        ecl[Row+NX*Column] = epl[Row+NX*Column];
                        ti[Row+NX*Column+NXNY*0] = 'c';
                    }
                }
                if(cells[BelowRow+NX*Column+NXNY*0] == 'i'){
                    if(PU_GPU(&state[Row+NX*Column]) < constant.probi){
                        cells[Row+NX*Column+NXNY*1] = 'e';
                        ecl[Row+NX*Column] = epl[Row+NX*Column];
                        ti[Row+NX*Column+NXNY*0] = 'c';
                    }
                }
                if(cells[AboveRow+NX*RightColumn+NXNY*0] == 'i'){
                    if(PU_GPU(&state[Row+NX*Column]) < constant.probi){
                        cells[Row+NX*Column+NXNY*1] = 'e';
                        ecl[Row+NX*Column] = epl[Row+NX*Column];
                        ti[Row+NX*Column+NXNY*0] = 'c';
                    }
                }
                if(cells[BelowRow+NX*LeftColumn+NXNY*0] == 'i'){
                    if(PU_GPU(&state[Row+NX*Column]) < constant.probi){
                        cells[Row+NX*Column+NXNY*1] = 'e';
                        ecl[Row+NX*Column] = epl[Row+NX*Column];
                        ti[Row+NX*Column+NXNY*0] = 'c';
                    }
                }
            }
            
            if(freecell == 1){
                // Cell free transmission
                double probablity = PU_GPU(&state[Row+NX*Column]);
                double pnoinfect = exp(-vtemp[Row+NX*Column+NXNY*0]*constant.beta*constant.timestep);
                if(probablity > pnoinfect){
	                cells[Row+NX*Column+NXNY*1] = 'e';
	                ecl[Row+NX*Column] = epl[Row+NX*Column];
			        if(ti[Row+NX*Column] == 'c'){
				        ti[Row+NX*Column] = 'b';
			        }
			        else{
				        ti[Row+NX*Column] = 'f';
			        }
                }
            }
        }

        //The Universal Time for the cells is kept here (ut)
        ut[Row+NX*Column] = ut[Row+NX*Column] + constant.timestep;
        vtemp[Row+NX*Column+NXNY*0] = vtemp[Row+NX*Column+NXNY*1];
        cells[Row+NX*Column+NXNY*0] = cells[Row+NX*Column+NXNY*1];
    }
}
