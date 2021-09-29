#include "xypic.cuh"

extern long int seed;
extern int ngx,ngy;
extern int Gsize;
extern int Csize;
extern int nsp;
extern int device_num;
extern int CondNUMR;
extern cudaDeviceProp prop;
extern __global__ void PCG(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, double *d_result);
extern __global__ void DepositAtom(int Gsize, int ngy, Species *info, GCP *sp, GPG *data);
extern __global__ void Cond_Sigma_Lap(int ngx, int ngy, float dx, float dy, float zlength, GGA *vecG, GCA *vecC, float *Phi, float *Sigma);
extern __global__ void GCondAInit(int CondNUMR, int nsp, float Value, GCondA *data);
extern __global__ void MoveE_Basic(int Gsize,int ngy,float dt_dx,float dt_dy, Species *info, GCP *sp, GPG *data, GGA *Field);
extern __global__ void SortBoundary_Basic(int Gsize,int ngy,float dt_dx,float dt_dy,int *StructureIndex, Species *info, GCP *sp, GPG *data, GGA *Field, GCondA *Cond, int *ReArgFlag);
#ifndef __CUDA_INIT_CUH__
#define __CUDA_INIT_CUH__
curandState *devStates;
int REDUCTION_GRID, MCC_GRID;
int REDUCTION_BLOCK, MCC_BLOCK;
int sMemSize;
dim3 FIELD_GRID,FIELD_BLOCK;
dim3 FIELD_GRID2,FIELD_BLOCK2;
dim3 DEPOSIT_GRID,DEPOSIT_BLOCK;
dim3 EFIELD_GRID,EFIELD_BLOCK;
dim3 MOVE_GRID, MOVE_BLOCK;
dim3 SORT_GRID, SORT_BLOCK;
float *Host_G_buf, *Host_C_buf;
void Set_Device_Parameter();
void Set_DiagParameter_cuda();
__global__ void SetSeed(curandState *state,long int seed,int num);
__global__ void MyKernel(int *array, int arrayCount);
#endif