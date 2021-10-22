#include "cuda_Init.cuh"
#define THREADS_PER_BLOCK 512   

void Set_Device_Parameter(){
    int grid,block;
    int mingrid;
    int numBlocksPerSm;
    int numThreads;
    int numSms;
    int size;
    // Find good grid and block size
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)SetSeed,0,nsp*Gsize); 
    grid = (nsp*Gsize + block - 1) / block;
    cudaMalloc((void**) &devStates, nsp*Gsize*sizeof(curandState));
    SetSeed<<<grid,block>>>(devStates,seed,nsp*Gsize); // Each thread gets same seed
    printf(" Find good grids and blocks. \n");
    cudaMalloc((void**) &dev_vsave, h_nvel * sizeof(float));
    cudaMemcpy(dev_vsave, vsave, h_nvel * sizeof(float),cudaMemcpyHostToDevice);
    // Field Solver 
    sMemSize = sizeof(double) * THREADS_PER_BLOCK;
    numBlocksPerSm = 0;
    numThreads = THREADS_PER_BLOCK;
    checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, PCG, numThreads, sMemSize));
    numSms = prop.multiProcessorCount;
    FIELD_GRID = dim3(numSms*numBlocksPerSm, 1, 1);
    FIELD_BLOCK = dim3(THREADS_PER_BLOCK, 1, 1);
    printf(" - Field Solver : [%d][%d]\n",numSms*numBlocksPerSm,THREADS_PER_BLOCK);
    printf("   Cooperative_groups = [%d]\n",numSms*numBlocksPerSm*THREADS_PER_BLOCK);
    // Field2
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)Cond_Sigma_Lap,0,Gsize*nsp); 
    grid = (Gsize*nsp + block - 1) / block;
    printf(" - Field module : [%d][%d]\n",grid,block);
    FIELD_GRID2 = dim3(grid, 1, 1);
    FIELD_BLOCK2 = dim3(block, 1, 1);
    // Deposit 
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)DepositAtom,0,Gsize*nsp); 
    grid = (Gsize*nsp + block - 1) / block;
    printf(" - Deposit module : [%d][%d]\n",grid,block);
    DEPOSIT_GRID = dim3(grid, 1, 1);
    DEPOSIT_BLOCK = dim3(block, 1, 1);
    // Efield
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)GCondAInit,0,CondNUMR); 
    grid = (CondNUMR*nsp + block - 1) / block;
    printf(" - Efield module : [%d][%d]\n",grid,block);
    EFIELD_GRID = dim3(grid, 1, 1);
    EFIELD_BLOCK = dim3(block, 1, 1);
    // Move 
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)MoveE_Basic,0,Gsize*nsp); 
    grid = (Gsize*nsp + block - 1) / block;
    printf(" - Move module : [%d][%d]\n",grid,block);
    MOVE_GRID = dim3(grid, 1, 1);
    MOVE_BLOCK = dim3(block, 1, 1);
    // SortBoundary
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)SortBoundary_Basic,0,Gsize*nsp); 
    grid = (Gsize*nsp + block - 1) / block;
    printf(" - SORT module : [%d][%d]\n",grid,block);
    SORT_GRID = dim3(grid, 1, 1);
    SORT_BLOCK = dim3(block, 1, 1);
    // MCC
    size = nsp * Gsize;
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)MCC_Ar_Basic,0,size); 
    grid = (size + block - 1) / block;
    printf(" - MCC module : [%d][%d]\n",grid,block);
    MCC_GRID = dim3(grid, 1, 1);
    MCC_BLOCK = dim3(block, 1, 1);
    // time setting
    gputime_field	=0.0;	gputime_efield	=0.0;	gputime_diag	=0.0;	gputime_move	=0.0;
	gputime_mcc		=0.0;	gputime_continue=0.0;	gputime_deposit	=0.0;	gputime_sort	=0.0;
	gputime_Tec	    =0.0;	gputime_dump	=0.0;	totaltime		=0.0;	TotalT_D		=0;
	TotalT_H		=0;	TotalT_M		=0;	TotalT_S		=0;
}
__global__ void SetSeed(curandState *state,long int seed,int num)
{
	int TID = blockDim.x * blockIdx.x + threadIdx.x;
	if(TID>=num) return;
    /* Each thread gets same seed, a different sequence number, no offset */
    curand_init(seed, TID, 0, &state[TID]);
}
