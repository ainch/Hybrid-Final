#include "xypic.cuh"
extern int Gsize,ngy,ngx;
extern int nsp;
extern Species *dev_info_sp;// particle species
extern GCP *dev_sp;
extern GPG *dev_G_sp;
extern GPG *Host_G_sp;
extern GGA *dev_GvecSet;
extern dim3 DEPOSIT_GRID,DEPOSIT_BLOCK;
extern int N_smt;
extern float *dev_Sigma, *dev_Source;
extern int   *dev_A_idx;
extern float *dev_R;
extern float *Host_G_buf, *Host_C_buf;
#ifndef __CUDA_MCC_CUH__
#define __CUDA_MCC_CUH__
void Deposit_cuda();
__global__ void DepositInitDensity(int Gsize, Species *info, GPG *data);
__global__ void DepositAtom(int Gsize, int ngy, Species *info, GCP *sp, GPG *data);
__global__ void DepositBoundary(int Gsize, int ngy, int nsp, GGA *vecSet, GPG *data);
__global__ void Smooth_121_A(int ngx, int ngy, int nsp, GGA *vecSet, GPG *data);
__global__ void Smooth_121_B(int ngx, int ngy, int nsp, GGA *vecSet, GPG *data);
__global__ void SumSource(int ngx, int ngy, Species *info, GGA *vecSet, GPG *data, float *Sigma, float *Source);
__global__ void PCG_Set(int Gsize, int *IDX, GGA *vecSet, float *Source, float *B);
#endif