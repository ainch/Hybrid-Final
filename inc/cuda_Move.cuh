#include "xypic.cuh"
extern int Gsize,ngy;
extern float dt_dx, dt_dy;
extern Species *SP;
extern GCP *dev_sp;
extern GPG *dev_G_sp;
extern GGA *dev_GvecSet;
extern dim3 MOVE_GRID, MOVE_BLOCK;
#ifndef __CUDA_MOVE_CUH__
#define __CUDA_MOVE_CUH__
void Move_cuda();
__global__ void MoveE_Basic(int Gsize,int ngy,float dt_dx,float dt_dy, Species *info, GCP *sp, GPG *data, GGA *Field);
#endif