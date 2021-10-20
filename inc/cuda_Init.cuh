#include "xypic.cuh"

extern long int seed;
extern int MainGas;
extern int ngx,ngy;
extern int Gsize;
extern int Csize;
extern int nsp;
extern int device_num;
extern int CondNUMR;
extern cudaDeviceProp prop;
extern int h_nvel;
extern float *vsave;
extern __global__ void PCG(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, double *d_result);
extern __global__ void DepositAtom(int Gsize, int ngy, Species *info, GCP *sp, GPG *data);
extern __global__ void Cond_Sigma_Lap(int ngx, int ngy, float dx, float dy, float zlength, GGA *vecG, GCA *vecC, float *Phi, float *Sigma);
extern __global__ void GCondAInit(int CondNUMR, int nsp, float Value, GCondA *data);
extern __global__ void MoveE_Basic(int Gsize,int ngy,float dt_dx,float dt_dy, Species *info, GCP *sp, GPG *data, GGA *Field);
extern __global__ void SortBoundary_Basic(int Gsize,int ngy,float dt_dx,float dt_dy,int *StructureIndex, Species *info, GCP *sp, GPG *data, GGA *Field, GCondA *Cond, int *ReArgFlag);
extern __global__ void MCC_Ar_Basic(int Gsize, int Csize, int ngy, int nsp, float dt, int MCCn, float dtm,float idx,float idy, int nvel, float *vsave,
											curandState *states, int N_LOGX, float idLOGX, MCC_sigmav *sigv, CollF *CollP, ArCollD *CX,  int TnRct, float*MCCR,
											Fluid *infoF, GFC *Fluid, GGA *BG, Species *info, GPG *data, GCP *sp);
namespace cg = cooperative_groups;
#ifndef __CUDA_INIT_CUH__
#define __CUDA_INIT_CUH__
curandState *devStates;
float *dev_vsave;
int sMemSize;
int sMemSize_MCC;
dim3 FIELD_GRID,FIELD_BLOCK;
dim3 FIELD_GRID2,FIELD_BLOCK2;
dim3 DEPOSIT_GRID,DEPOSIT_BLOCK;
dim3 EFIELD_GRID,EFIELD_BLOCK;
dim3 MOVE_GRID, MOVE_BLOCK;
dim3 SORT_GRID, SORT_BLOCK;
dim3 MCC_GRID, MCC_BLOCK;
float *Host_G_buf, *Host_C_buf;
void Set_Device_Parameter();
void Set_DiagParameter_cuda();
__global__ void SetSeed(curandState *state,long int seed,int num);
__global__ void MyKernel(int *array, int arrayCount);
#endif