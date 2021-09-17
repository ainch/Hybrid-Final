#include "cuda_Init.cuh"
#define THREADS_PER_BLOCK 512   


void Set_DiagParameter_cuda(){
    // Host BUF VECTOR
    Host_G_buf = VFMalloc(Gsize);
    Host_C_buf = VFMalloc(Csize);
    VFInit(Host_G_buf,0.0,Gsize);
    VFInit(Host_C_buf,0.0,Csize);
}
void Set_Device_Parameter(){
    int grid,block;
    int mingrid;
    int numBlocksPerSm;
    int numThreads;
    int numSms;
    // Find good grid and block size
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)SetSeed,0,Gsize); 
    grid = (Gsize + block - 1) / block;
    cudaMalloc((void**) &devStates, Gsize * sizeof(curandState));
    SetSeed<<<grid,block>>>(devStates,seed,Gsize); // Each thread gets same seed
    //
    printf(" Find good grids and blocks. \n");
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

    // Example : Find good grid and block size
    int Search_Occupancy_Flag = 0;
    if(Search_Occupancy_Flag){
        //Example
        int *array;
        int blockSize;      // The launch configurator returned block size
        int minGridSize;    // The minimum grid size needed to achieve the
                            // maximum occupancy for a full device
                            // launch
        int gridSize;       // The actual grid size needed, based on input
                            // size
        cudaDeviceProp prop;
        int numBlocks;       // Occupancy in terms of active blocks
        int activeWarps;
        int maxWarps;
        cudaGetDevice(&device_num);
        cudaGetDeviceProperties(&prop, device_num);
        cudaMalloc((void**) &array, Gsize * sizeof(int));
        cudaMemset((void *) array, 0, Gsize * sizeof(int));
        // Get MiniGridSize and blockSize
        cudaOccupancyMaxPotentialBlockSize(&minGridSize,&blockSize,(void*)MyKernel,0,Gsize); 
        // Get Occupancy
        cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocks,MyKernel,blockSize,0);
        activeWarps = numBlocks * blockSize / prop.warpSize;
        maxWarps = prop.maxThreadsPerMultiProcessor / prop.warpSize;
        printf("Occupancy: = %3.0f %\n",(double)activeWarps / maxWarps * 100);
        printf("maxWarps = %d\n",maxWarps);
        // Round up according to array size
        gridSize = (Gsize + blockSize - 1) / blockSize;
        printf("minGridSize = %d\n",minGridSize);
        printf("blockSize = %d\n",blockSize);
        printf("gridSize = %d\n",gridSize);
        MyKernel<<<gridSize, blockSize>>>(array, Gsize);
        cudaDeviceSynchronize();
    }
}
__global__ void SetSeed(curandState *state,long int seed,int num)
{
		int TID = blockDim.x * blockIdx.x + threadIdx.x;
		if(TID>num) return;
    /* Each thread gets same seed, a different sequence number, no offset */
    curand_init(seed, TID, 0, &state[TID]);
}
__global__ void MyKernel(int *array, int arrayCount)
{
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    if (idx < arrayCount) {
        array[idx] *= array[idx];
    }
}