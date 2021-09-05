#include "cuda_Field.cuh"
// FOR Field method 4
#define THREADS_PER_BLOCK 512   
//
void PCG_SOLVER_Laplace(){
    // OUTPUT
    // Lap_TEMP_Sol[Gsize] : Temperature Profile
    // Lap_PHI_Sol[CondNUMR][Gsize] : Each of conductor Phi Profile, This is Device value
    // Lap_SIG_Sol[CondNUMR][CondNUMR] : Each of conductor Sigma Profile for external circuit

    //Make a Field DATA set  
    //Host_PCG_DATA = (DPS_Data*)malloc(A_size*sizeof(DPS_Data));
    //checkCudaErrors(cudaMalloc((void**)&dev_PCG_DATA, A_size*sizeof(DPS_Data)));
    //Make_PCG_DATA_Init<<<A_size/4,4>>>(dev_PCG_DATA,A_size,dev_M);
    //checkCudaErrors(cudaMemcpy(Host_PCG_DATA, dev_PCG_DATA, A_size*sizeof(DPS_Data), cudaMemcpyDeviceToHost));
    // Laplace Solution
    //cudaMallocPitch(&Lap_PHI_Sol, &pitch, Gsize * sizeof(float), CondNUMR); // for Laplace Solution
    //cudaMalloc((void**) &Lap_TEMP_Sol, Gsize * sizeof(int));
    // cudaMemset((void *) array, 0, Gsize * sizeof(int));
    

           // Find good grid and block size
        //cudaOccupancyMaxPotentialBlockSize(&mingrid,&FIELD_BLOCK,(void*)cusparseSpMV,0,Gsize); 
        //FIELD_GRID = (Gsize + FIELD_BLOCK - 1) / FIELD_BLOCK;
        //printf("blockSize = %d\n",FIELD_BLOCK);
        //printf("gridSize = %d\n",FIELD_GRID);
    // For test

    //if (matA       ) { checkCudaErrors(cusparseDestroySpMat(matA)); }
    //if (vecx       ) { checkCudaErrors(cusparseDestroyDnVec(vecx)); }
    //if (vecAP      ) { checkCudaErrors(cusparseDestroyDnVec(vecAP)); }
    //if (vecp       ) { checkCudaErrors(cusparseDestroyDnVec(vecp)); }
}
void Set_MatrixPCG_cuda(){

}
