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
extern dim3 FIELD_GRID,FIELD_BLOCK;
extern dim3 FIELD_GRID2,FIELD_BLOCK2;
extern dim3 DEPOSIT_GRID,DEPOSIT_BLOCK;
extern dim3 EFIELD_GRID,EFIELD_BLOCK;
extern dim3 MOVE_GRID, MOVE_BLOCK;
extern dim3 SORT_GRID, SORT_BLOCK;
extern dim3 MCC_GRID, MCC_BLOCK;
extern dim3 CONTI_GRID,CONTI_BLOCK;
extern float *TotPotential;
extern float *dev_Source;
extern float *dev_Sigma;
extern float *vec_Potential, *sum_Potential, *ave_Potential; // [Gsize] potential
extern float *vec_Source, *sum_Source, *ave_Source; // [Gsize] Charge density 
extern float *vec_Sigma, *sum_Sigma, *ave_Sigma; // [Gsize] [Dielectric] surface charge density, [Conductor] Surface current
extern float *sum_Ex, *ave_Ex, *sum_Ey, *ave_Ey;
extern GGA *vec_G,*dev_GvecSet;
extern GCA *vec_C,*dev_CvecSet;
#ifndef __CUDA_DIAGNOSTIC_CUH__
#define __CUDA_DIAGNOSTIC_CUH__
float *Host_G_buf, *Host_C_buf;
void Diagnostic();
void Set_Diagnostic_cuda();
__global__ void Accomulate_Particle_Density(int nsp, int Gsize, GPG *data);
__global__ void Average_Particle_Density(int nsp, int Gsize, int N_ave, Species *info, GPG *data);
#endif