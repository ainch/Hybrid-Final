#include "xypic.cuh"
extern float dt;
extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern int N_ave;
extern int Gsize,Csize;
extern int ngx,ngy,ncx,ncy;
extern int Conti_Flag;
extern int MainGas;
extern int nsp,nfsp;
extern Species *SP;// particle species
extern Species *dev_info_sp;// particle species
extern GPG *dev_G_sp;
extern int HISTMAX;
extern int dHIST;
extern int Hcount;
extern int hist_count;
extern Hist *HistPt,*HistFG;
extern float *t_array;  
extern float *iter_array;  
extern int *FIter;
extern int nave_count;
extern dim3 MCC_GRID, MCC_BLOCK;
#ifndef __CUDA_DIAGNOSTIC_CUH__
#define __CUDA_DIAGNOSTIC_CUH__
float *Host_G_buf, *Host_C_buf;
void Diagnostic();
void Set_Diagnostic_cuda();
__global__ void Accomulate_Particle_Density(int nsp, int Gsize, GPG *data);
__global__ void Average_Particle_Density(int nsp, int Gsize, int N_ave, Species *info, GPG *data);
#endif