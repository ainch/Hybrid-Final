#include "xypic.cuh"
extern int device_num;
extern int ngx,ngy,Gsize;
extern int Csize;
extern float dx, dy, zlength;
extern int A_size;
extern float PCGtol;
extern float PCGtol2;
extern GGA *vec_G;
extern GCA *vec_C;
extern int CondNUMR;
extern int Lap_Field_Solver_Test,Lap_Field_Solver_Flag,Lap_Field_Solver_Save;
extern int Preconditioner_Flag;
extern int FieldIter;
extern int A_size;
extern float *A_val,*TA_val;
extern int *Ai,*Aj;
extern int **A_idx;
extern float *MatM,**cond_b,*temp_b;
extern float *Host_G_buf, *Host_C_buf;
extern int sMemSize;
extern dim3 FIELD_GRID,FIELD_BLOCK;

extern "C" void Field_Laplace_Solution_Save(char *Filename,float **Sol);
extern "C" float *VFMalloc(int size);
extern "C" void VFInit(float *V,float C,int size);
extern "C" void VFCopy(float *V,float *C,int size);

#ifndef __CUDA_FIELD_CUH__
#define __CUDA_FIELD_CUH__
namespace cg = cooperative_groups;
int nz = 5*A_size;
int N = A_size;
float *dev_A,*dev_TA;	
int *dev_Aj,*dev_Ai;
float *dev_b,*dev_Tb;
int   *vec_A_idx;
int   *dev_A_idx;
float *dev_M;
float *dev_AP,*dev_X,*dev_R,*dev_Z,*dev_P;
double *dot_result;
int *FIter;
GGA *dev_GvecSet;
GCA *dev_CvecSet;
float *dev_phi;     // PCG Solution
float *dev_phi_buf; // Sigma or buf
size_t pitch;
float *Lap_TEMP_Sol; // Temperature Solution of Laplace Equation
float *Lap_PHI_Sol; // Each of conductor Phi Solution of Laplace Equation, This is Device value
float **Lap_SIG_Sol; // Each of conductor Sigma Solution of Laplace Equation for external circuit
void Set_MatrixPCG_cuda();
void PCG_SOLVER_Laplace();
__device__ void Mat_x_Vec(int *I, int *J, float *val, int nnz, int num_rows, float alpha, float *inputVecX, 
                        float *outputVecY, cg::thread_block &cta, const cg::grid_group &grid);
__device__ void A_x_X_p_Y(float a, float *x, float *y, int size, const cg::grid_group &grid);
__device__ void A_x_Y_p_X(float a, float *x, float *y, int size, const cg::grid_group &grid);
__device__ void Vec_Dot_Sum(float *vecA, float *vecB, double *result, int size, const cg::thread_block &cta, const cg::grid_group &grid);
__device__ void CopyVector(float *srcA, float *destB, int size, const cg::grid_group &grid);
__device__ void Vec_x_Vec(float *vecA, float *vecB, float *vecC, int size, const cg::thread_block &cta, const cg::grid_group &grid);
__global__ void PCG(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, double *d_result);
__global__ void PCG_Deposit_Lap(int Gsize, int *IDX, GGA *vecG, int k, float *X, float *PHI);
__global__ void PCG_Deposit_Temp(int Gsize, int *IDX, float *X, GGA *vecG);
__global__ void Cond_Sigma_Lap(int ngx, int ngy, float dx, float dy, float zlength, GGA *vecG, GCA *vecC, float *Phi, float *Sigma);
__global__ void SaveAT2D(float *A, size_t pitch, int height, float *PHI, int n);
__global__ void LoadAT2D(float *A, size_t pitch, int height, float *PHI, int n);
#endif