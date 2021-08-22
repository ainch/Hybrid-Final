#include "xypic.cuh"

extern int ngx,ngy,Gsize;
extern int A_size;
extern float PCGtol;
extern float PCGtol2;
extern HGA *vec_G;
extern int CondNUMR;
extern int Field_Solver_Flag;
extern int FieldIter;
extern int A_size;
extern float *A_val,*MatTA;
extern int *Ai,*Aj;
extern int **A_idx;
extern float *MatM,**cond_b,*temp_b;
extern float **phi_dw,**phi_u;
extern int FIELD_GRID, FIELD_BLOCK;
extern "C" void Field_Laplace_Solution_Save(char *Filename,float **Sol);
extern "C" float *VFMalloc(int size);
extern "C" void VFInit(float *V,float C,int size);
extern "C" void VFCopy(float *V,float *C,int size);

#ifndef __CUDA_FIELD_CUH__
#define __CUDA_FIELD_CUH__
cublasHandle_t cublasHandle;
cublasStatus_t cublasStatus;
cusparseHandle_t cusparseHandle;
 cusparseStatus_t cusparseStatus;
cusparseSpMatDescr_t matA;
cusparseDnVecDescr_t vecx;
cusparseDnVecDescr_t vecp;
cusparseDnVecDescr_t vecAP;
cudaStream_t stream1, streamForGraph;
cudaGraph_t initGraph;
cudaGraphExec_t graphExec;
float *Lap_TEMP_Sol; // Temperature Solution of Laplace Equation
float **Lap_PHI_Sol; // Each of conductor Phi Solution of Laplace Equation, This is Device value
float **Lap_SIG_Sol; // Each of conductor Sigma Solution of Laplace Equation for external circuit
size_t pitch;
// Single CPU PCG parameter
float *AX,*X,*B,*R0,*Z0,*P0,*AP,*PAP;
// ORIGIN GPU PCG parameter
float *dev_A;	
int *dev_Aj,*dev_Ai;
int   *vec_A_idx;
float *dev_M;
float *dev_AP,*dev_X,*dev_b,*dev_R,*dev_Z,*dev_P;
float *d_r1, *d_r0, *d_dot, *d_a, *d_na, *d_b;
DPS_Const *dev_PCG_const;
DPS_Const *Host_PCG_const;
DPS_Data *dev_PCG_DATA;
DPS_Data *Host_PCG_DATA;
void Field_Method0_Initial();
int CG_CPU();
void Field_Method1_Initial();
int PCG_CPU();
void Field_Method2_Initial();
int CG_GPU();
void Field_Method3_Initial();
int CG_GPU_CudaGraphs();
void Field_Method4_Initial();
void Field_Method5_Initial();
void Field_Method6_Initial();
void Set_MatrixPCG_cuda();
void PCG_SOLVER_Laplace();
int PCG_LAP_Divide(int grid,int block);
__global__ void Make_PCG_DATA_Init(DPS_Data *p, int size,float *MatrixM);
__global__ void Make_PCG_Const_Init(DPS_Const *p,int Asize, float tol);
__global__ void PCG_LAP(float *A,int *Ai,int *Aj, DPS_Const *PCG_C, DPS_Data *PCG_D,float *X,float *b);
__global__ void initVectors(float *rhs, float *x, int N); 
__global__ void r1_div_x(float *r1, float *r0, float *b);
__global__ void a_minus(float *a, float *na);
__global__ void SaveAT2D(float *A, size_t pitch, int height, float *PHI, int n);
__global__ void LoadAT2D(float *A, size_t pitch, int height, float *PHI, int n);
#endif