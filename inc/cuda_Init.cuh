#include "xypic.cuh"

extern long int seed;
extern int ngx,ngy;
extern int Gsize;
extern int Csize;
extern int nsp;
extern int device_num;
extern cudaDeviceProp prop;
extern __global__ void PCG(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, double *d_result);
extern __global__ void DepositAtom(int Gsize, int ngy, Species *info, GCP *sp, GPG *data);
#ifndef __CUDA_INIT_CUH__
#define __CUDA_INIT_CUH__
curandState *devStates;
int REDUCTION_GRID, EFIELD_GRID, MOVE_GRID, SORT_GRID, MCC_GRID;
int REDUCTION_BLOCK, EFIELD_BLOCK, MOVE_BLOCK, SORT_BLOCK,MCC_BLOCK;
int sMemSize;
dim3 FIELD_GRID,FIELD_BLOCK;
dim3 DEPOSIT_GRID,DEPOSIT_BLOCK;
float *Host_G_buf, *Host_C_buf;
void Set_Device_Parameter();
void Set_DiagParameter_cuda();
__global__ void SetSeed(curandState *state,long int seed,int num);
__global__ void MyKernel(int *array, int arrayCount);
#endif