#include "xypic.cuh"

extern int ngx,ngy,Gsize;
extern int A_size;
extern float PCGtol;
extern HGA *vec_G;
extern int CondNUMR;
extern int PCG_Method;
extern int FieldIter;
extern int A_size;
extern float *A_val,*MatTA;
extern int *Ai,*Aj;
extern int **A_idx;
extern float *MatM,**cond_b,*temp_b;
extern float **phi_dw,**phi_u;
extern "C" void CPU_PCG_Laplace_Solution_Save(float **Sol);
extern "C" float *VFMalloc(int size);
extern "C" void VFInit(float *V,float C,int size);
extern "C" void VFCopy(float *V,float *C,int size);

#ifndef __CUDA_FIELD_CUH__
#define __CUDA_FIELD_CUH__
float *Lap_TEMP_Sol; // Temperature Solution of Laplace Equation
float **Lap_PHI_Sol; // Each of conductor Phi Solution of Laplace Equation, This is Device value
float **Lap_SIG_Sol; // Each of conductor Sigma Solution of Laplace Equation for external circuit
size_t pitch;
float rsold,rnew,Temp;
float alpha,beta;
// Single CPU PCG parameter
float *AX,*X,*B,*R0,*Z0,*P0,*AP,*PAP;
// ORIGIN GPU PCG parameter
float *dev_A;	
int *dev_Aj,*dev_Ai;
int   *vec_A_idx;
float *dev_M;
float *dev_AP,*dev_X,*dev_b,*dev_R,*dev_Z,*dev_P;
DPS_Const *dev_PCG_const;
DPS_Const *Host_PCG_const;
DPS_Data *dev_PCG_DATA;
DPS_Data *Host_PCG_DATA;
void Set_MatrixPCG_cuda();
void PCG_SOLVER_Laplace();
int PCG_SINGLECPU();
__global__ void Make_PCG_DATA_Init(DPS_Data *p, int size,float *MatrixM);
__global__ void Make_PCG_Const_Init(DPS_Const *p,int Asize, float tol);
__global__ void PCG_LAP(int *Iter,int Gsize,int Asize,float *A,int *Ai,int *Aj,float *M,float *AP,float *R,float *Z,float *P,float *X,float *b);
__global__ void SaveAT2D(float *A, size_t pitch, int height, float *PHI, int n);
__global__ void LoadAT2D(float *A, size_t pitch, int height, float *PHI, int n);
#endif