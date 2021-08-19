#include "cuda_Init.cuh"

void Set_NullCollisionTime_cuda(){
    
}
void Set_DiagParameter_cuda(){
    
}
void Set_Particle_cuda(){
    
}
void Set_Device_Parameter(){
    int grid,block;
    int mingrid;
    
    int   *vec_cond_Garray;
    int   *vec_boundary_Garray;
    int   *vec_face_Garray;
    float *vec_area_Garray;
    float *vec_eps_Carray;
    float *dev_Sigma;
    int   *dev_face_Garray;
    float *dev_area_Garray;
    float *dev_eps_Carray;
    // Find good grid and block size
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)SetSeed,0,Gsize); 
    grid = (Gsize + block - 1) / block;
    

    cudaMalloc((void**) &devStates, Gsize * sizeof(curandState));
    SetSeed<<<grid,block>>>(devStates,seed,Gsize); // Each thread gets same seed

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
		int TID=blockDim.x*blockIdx.x+threadIdx.x;
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