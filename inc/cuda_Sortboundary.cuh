#include "xypic.cuh"
extern int ncx,ncy;
extern int Gsize,ngy;
extern float dt_dx,dt_dy;
extern Species *SP;
extern GCP *dev_sp;
extern GPG *dev_G_sp;
extern GGA *dev_GvecSet;
extern GCondA *dev_CondVec;
extern int *vec_StructureIndex;
extern dim3 SORT_GRID, SORT_BLOCK;
#ifndef __CUDA_MOVE_CUH__
#define __CUDA_MOVE_CUH__
int *dev_StructureIndex;
int *dev_ReArgFlag;
void Set_SortBoundary_cuda();
void SortBounndary_cuda();
__global__ void SortBoundary_Basic(int Gsize,int ngy,float dt_dx,float dt_dy,int *StructureIndex, Species *info, GCP *sp, GPG *data, GGA *Field, GCondA *Cond, int *ReArgFlag);
#endif