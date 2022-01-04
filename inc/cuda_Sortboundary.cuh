#include "xypic.cuh"
extern int ncx,ncy;
extern int Gsize,ngy;
extern int CondNUMR;
extern int nsp;
extern float dt_dx,dt_dy;
extern Species *dev_info_sp;// particle species
extern GCP *dev_sp;
extern GPG *dev_G_sp;
extern GGA *dev_GvecSet;
extern float *dev_CondCharge;
extern int *vec_StructureIndex;
extern dim3 SORT_GRID, SORT_BLOCK;
extern curandState *devStates;
extern void Deposit_cuda();
#ifndef __CUDA_MOVE_CUH__
#define __CUDA_MOVE_CUH__
extern int *dev_StructureIndex;
extern int *ReArgFlag;
void Set_SortBoundary_cuda();
void SortBounndary_cuda();
__global__ void ReArgPtcls(int Gsize, int nsp, int ngy, curandState *states, Species *info, GCP *sp, GPG *data);
__global__ void SortBoundary_Basic(int Gsize,int ngy,int CondNum, float dt_dx,float dt_dy,int *StructureIndex, Species *info, GCP *sp, GPG *data, GGA *Field, float *Cond, int *ReArgFlag);
#endif