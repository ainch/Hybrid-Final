#include "xypic.cuh"

extern long int seed;
extern int MainGas;
extern int ngx,ngy;
extern int Gsize;
extern int Csize;
extern int nsp,nfsp;
extern int device_num;
extern int CondNUMR;
extern cudaDeviceProp prop;
extern int h_nvel;
extern float *vsave;
extern __global__ void PCG(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, double *d_result);
extern __global__ void PCG_float(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, float *d_result);
extern __global__ void DepositAtom(int Gsize, int ngy, Species *info, GCP *sp, GPG *data);
extern __global__ void Cond_Sigma_Lap(int ngx, int ngy, float dx, float dy, float zlength, GGA *vecG, GCA *vecC, float *Phi, float *Sigma);
extern __global__ void MoveE_Basic(int Gsize,int ngy,float dt_dx,float dt_dy, Species *info, GCP *sp, GPG *data, GGA *Field);
extern __global__ void SortBoundary_Basic(int Gsize,int ngy,int CondNum, float dt_dx,float dt_dy,int *StructureIndex, Species *info, GCP *sp, GPG *data, GGA *Field, float *Cond, int *ReArgFlag);
extern __global__ void MCC_Ar_Basic(int Gsize, int Csize, int ngy, int nsp, float dt, int MCCn, float dtm,float idx,float idy, int nvel, float *vsave,
											curandState *states, int N_LOGX, float idLOGX, MCC_sigmav *sigv, CollF *CollP, ArCollD *CX,  int TnRct, float *MCCR, float *Stack, 
											Fluid *infoF, GFG *Fluid, GGA *BG, Species *info, GPG *data, GCP *sp);
extern __global__ void Average_Particle_Density(int nsp, int Gsize, int N_ave, Species *info, GPG *data);
extern __global__ void Average_Argon_MCC_rate(int Gsize, int TnRct, int N_ave, float dt, float *MCCR, float *ave_MCCR, Species *info);
namespace cg = cooperative_groups;
extern float 	totaltime,gputime_field,gputime_efield;
extern float 	gputime_move,gputime_mcc,gputime_deposit;
extern float 	gputime_diag,gputime_sort,gputime_Tec;
extern float 	gputime_continue,gputime_dump;
extern int		TotalT_D;
extern int		TotalT_H;
extern int		TotalT_M;
extern int		TotalT_S;
#ifndef __CUDA_INIT_CUH__
#define __CUDA_INIT_CUH__
extern curandState *devStates;
extern float *dev_vsave;
extern int sMemSize;
extern int sMemSize_MCC;
extern dim3 FIELD_GRID,FIELD_BLOCK;
extern dim3 FIELD_GRID2,FIELD_BLOCK2;
extern dim3 DEPOSIT_GRID,DEPOSIT_BLOCK;
extern dim3 EFIELD_GRID,EFIELD_BLOCK;
extern dim3 MOVE_GRID, MOVE_BLOCK;
extern dim3 SORT_GRID, SORT_BLOCK;
extern dim3 MCC_GRID, MCC_BLOCK;
extern dim3 DIAG_G_GRID, DIAG_G_BLOCK;
extern dim3 DIAG_NSPG_GRID, DIAG_NSPG_BLOCK;
void Set_Device_Parameter();
__global__ void SetSeed(curandState *state,long int seed,int num);
__global__ void MyKernel(int *array, int arrayCount);
#endif