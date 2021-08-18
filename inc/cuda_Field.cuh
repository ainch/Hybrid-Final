#include "xypic.cuh"

extern int ngx,ngy,Gsize;
extern int A_size;
extern float PCGtol;
extern int CondNUMR;
extern int PCG_Method;
extern int FieldIter;
extern int A_size;
extern float *A_val,*MatTA;
extern int *Ai,*Aj;
extern int **A_idx;
extern float *MatM,**cond_b,*temp_b;
extern float **phi_dw,**phi_u;
extern float *VFMalloc(int size);
extern void VFInit(float *V,float C,int size);
extern void VFCopy(float *V,float *C,int size);

#ifndef __CUDA_FIELD_CUH__
#define __CUDA_FIELD_CUH__
float *Lap_TEMP_Sol; // Temperature Solution of Laplace Equation
float **Lap_PHI_Sol; // Each of conductor Phi Solution of Laplace Equation, This is Device value
float **Lap_SIG_Sol; // Each of conductor Sigma Solution of Laplace Equation for external circuit
size_t pitch;
void Set_MatrixPCG_cuda();
void PCG_SOLVER_Laplace();
__global__ void PCG(int Iter,int Gsize,int Asize,float *A,int *Ai,int *Aj,float *X,float *b);
__global__ void SaveAT2D(float *A, size_t pitch, int height, float *PHI, int n);
__global__ void LoadAT2D(float *A, size_t pitch, int height, float *PHI, int n);
#endif