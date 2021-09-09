#include "xypic.cuh"

extern long int seed;
extern int ngx,ngy;
extern int Gsize;
extern int Csize;
extern int device_num;
extern cudaDeviceProp prop;
extern __global__ void PCG(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, double *d_result);
#ifndef __CUDA_INIT_CUH__
#define __CUDA_INIT_CUH__
curandState *devStates;
int REDUCTION_GRID, EFIELD_GRID, MOVE_GRID, SORT_GRID, MCC_GRID,DEPOSIT_GRID;
int REDUCTION_BLOCK, EFIELD_BLOCK, MOVE_BLOCK, SORT_BLOCK,MCC_BLOCK, DEPOSIT_BLOCK;
int sMemSize;
dim3 FIELD_GRID,FIELD_BLOCK;
float *Host_G_buf, *Host_C_buf;
void Set_Device_Parameter();
void Set_Particle_cuda();
void Set_NullCollisionTime_cuda();
void Set_DiagParameter_cuda();
__global__ void SetSeed(curandState *state,long int seed,int num);
__global__ void MyKernel(int *array, int arrayCount);
#endif