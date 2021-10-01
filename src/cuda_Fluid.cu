#include "cuda_Fluid.cuh"

void Set_Fluid_cuda(){

    checkCudaErrors(cudaMalloc((void**)&dev_FG, nfsp * sizeof(Fluid)));
    checkCudaErrors(cudaMemcpy(dev_FG, FG, nfsp * sizeof(Fluid), cudaMemcpyHostToDevice));
    
    checkCudaErrors(cudaMalloc((void**)&dev_C_F, nfsp * Csize *  sizeof(GFC)));
    checkCudaErrors(cudaMemcpy(dev_C_F, Host_C_F, nfsp * Csize *  sizeof(GFC), cudaMemcpyHostToDevice));
}