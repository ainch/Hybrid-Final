#include "xypic.cuh"
extern int Gsize,ngy,ngx;
extern int nsp;
extern Species *dev_info_sp;// particle species
extern GCP *dev_sp;
extern GPG *dev_G_sp;
extern GPG *Host_G_sp;
extern GGA *dev_GvecSet;
extern dim3 DEPOSIT_GRID,DEPOSIT_BLOCK;
#ifndef __CUDA_MCC_CUH__
#define __CUDA_MCC_CUH__
void Deposit_cuda();
__global__ void DepositAtom(int Gsize, int ngy, Species *info, GCP *sp, GPG *data);
__global__ void DepositBoundary(int Gsize, int ngy, int nsp, GGA *vecSet, GPG *data);
#endif