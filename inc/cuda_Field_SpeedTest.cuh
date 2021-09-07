#include "xypic.cuh"

extern int device_num;
extern int ngx,ngy,Gsize;
extern int A_size;
extern float PCGtol;
extern float PCGtol2;
extern HGA *vec_G;
extern int CondNUMR;
extern int Lap_Field_Solver_Test,Lap_Field_Solver_Flag,Lap_Field_Solver_Save;
extern int FieldIter;
extern int A_size;
extern float *A_val,*TA_val;
extern int *Ai,*Aj;
extern int **A_idx;
extern float *MatM,**cond_b,*temp_b;
extern float **phi_dw,**phi_u;
extern int FIELD_GRID, FIELD_BLOCK;
extern float *dev_A;	
extern int *dev_Aj,*dev_Ai;
extern int   *vec_A_idx;
extern float *dev_M;
extern float *dev_AP,*dev_X,*dev_b,*dev_R,*dev_Z,*dev_P;
extern "C" void Field_Laplace_Solution_Save(char *Filename,float **Sol);
extern "C" float *VFMalloc(int size);
extern "C" void VFInit(float *V,float C,int size);
extern "C" void VFCopy(float *V,float *C,int size);

#ifndef __CUDA_FIELD_TEST_CUH__
#define __CUDA_FIELD_TEST_CUH__
cusparseSpMatDescr_t matA;
cusparseDnVecDescr_t vecx;
cusparseDnVecDescr_t vecp;
cusparseDnVecDescr_t vecAP;

cudaGraph_t initGraph;
cudaGraphExec_t graphExec;

size_t pitch;
// ORIGIN GPU PCG parameter
DPS_Const *dev_PCG_const;
DPS_Const *Host_PCG_const;
DPS_Data *dev_PCG_DATA;
DPS_Data *Host_PCG_DATA;
//
void PCG_Laplace_TEST();
void Field_Method0_Initial();
int CG_CPU();
float *AX,*X,*B,*R0,*Z0,*P0,*AP,*PAP;
// FOR Field method 1
void Field_Method1_Initial();
int PCG_CPU();
// FOR Field method 2
void Field_Method2_Initial();
int CG_GPU();

// FOR Field method 3
void Field_Method3_Initial();
cublasHandle_t cublasHandle;
cublasStatus_t cublasStatus;
cusparseHandle_t cusparseHandle;
cusparseStatus_t cusparseStatus;
cudaStream_t stream1, streamForGraph;
int CG_GPU_CudaGraphs();
__global__ void initVectors(float *rhs, float *x, int N); 
__global__ void r1_div_x(float *r1, float *r0, float *b);
__global__ void a_minus(float *a, float *na);
// FOR Field method 4
void Field_Method4_Initial();
namespace cg = cooperative_groups;
__device__ void gpuSaxpy(float *x, float *y, float a, int size, const cg::grid_group &grid);
__device__ void gpuSpMV(int *I, int *J, float *val, int nnz, int num_rows, float alpha, float *inputVecX, 
                        float *outputVecY, cg::thread_block &cta, const cg::grid_group &grid);
__device__ void gpuSaxpy(float *x, float *y, float a, int size, const cg::grid_group &grid);
__device__ void gpuRSaxpy(float *x, float *y, float a, int size, const cg::grid_group &grid);
__device__ void gpuDotProduct(float *vecA, float *vecB, double *result, int size, const cg::thread_block &cta, const cg::grid_group &grid);
__device__ void gpuCopyVector(float *srcA, float *destB, int size, const cg::grid_group &grid);
__device__ void gpuScaleVector(float *vec, float alpha, int size, const cg::grid_group &grid);
__global__ void gpuConjugateGradient(int *I, int *J, float *val, float *x,  float *Ax, float *p, float *r,DPS_Const *result,double *d_result);
// FOR Field method 5
void Field_Method5_Initial();
// namespace cg = cooperative_groups;
__device__ void gpuProductVector(float *vecA, float *vecB, float *vecC, int size, const cg::thread_block &cta, const cg::grid_group &grid);
__global__ void gpuPreConjugateGradient(int *I, int *J, float *val, float *M, float *x,  float *Ax, float *p, float *r, float *Z, 
            DPS_Const *result,double *d_result);
// FOR Field method 6
void Field_Method6_Initial();
// namespace cg = cooperative_groups;
cudaLaunchParams *launchParamsList;
int GPUn,*deviceN;
int *man_I,*man_J;
float *man_A;
float *man_R, *man_P, *man_AP, *man_Z, *man_M, *man_X;
__device__ double grid_dot_result = 0.0;
__device__ void setDotResultToZero(double *dot_result);
__device__ void MultigpuProductVector(float *vecA, float *vecB, float *vecC, int size, const cg::thread_block &cta, const cg::multi_grid_group &multi_grid);
__device__ void MultigpuScaleVector(float *vec, float alpha, int size, const cg::multi_grid_group &multi_grid);
__device__ void MultigpuCopyVector(float *srcA, float *destB, int size, const cg::multi_grid_group &multi_grid);
__device__ void MultigpuDotProduct(float *vecA, float *vecB, int size, const cg::thread_block &cta, const cg::multi_grid_group &multi_grid);
__device__ void MultigpuRSaxpy(float *x, float *y, float a, int size, const cg::multi_grid_group &multi_grid);
__device__ void MultigpuSaxpy(float *x, float *y, float a, int size, const cg::multi_grid_group &multi_grid);
__device__ void MultigpuSpMV(int *I, int *J, float *val, int nnz, int num_rows, float alpha, float *inputVecX, 
                        float *outputVecY, cg::thread_block &cta, const cg::multi_grid_group &multi_grid);
__global__ void multiGpuPreConjugateGradient(int *I, int *J, float *val, float *M, float *x,  float *Ax, float *p, float *r, float *Z, 
            int nnz, int N,float tol,double *d_result);
//
void Field_Method7_Initial();

//
void Field_Method8_Initial();
float *dev_Y;
float *dev_L,*dev_U;
//
int PCG_LAP_Divide(int grid,int block);
__global__ void Make_PCG_DATA_Init(DPS_Data *p, int size,float *MatrixM);
__global__ void Make_PCG_Const_Init(DPS_Const *p,int Asize, float tol);
__global__ void PCG_LAP(float *A,int *Ai,int *Aj, DPS_Const *PCG_C, DPS_Data *PCG_D,float *X,float *b);
#endif