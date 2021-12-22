#include "xypic.cuh"
extern double t;
extern int DumpFlag;
extern int device_num;
extern int ngx,ngy,Gsize;
extern int Csize;
extern float dx, dy, zlength;
extern float idxy,idx,idy,dx2,dy2,dxdy2,hdx,hdy,r_eps0;
extern int A_size;
extern float PCGtol;
extern float PCGtol2;
extern GGA *vec_G;
extern GCA *vec_C;
extern int CondNUMR;
extern int External_Flag; // 0 : Voltage driven, 1: Power driven
extern int SrcNUM;
extern int *SrcM_ID;
extern float *SrcDC, *SrcPOWER, *SrcAC, *SrcFREQ, *Src2piFREQ, *SrcPHASE, *SrcRPHASE;
extern int *Efield_Flag, *Cond_Source_num, *Cond_count,**Cond_Power_ID;
extern float *Cond_Power;
extern float **CC_a; //Circuit_Const a
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
extern dim3 FIELD_GRID2,FIELD_BLOCK2;
extern dim3 EFIELD_GRID,EFIELD_BLOCK;
extern int MainGas;
extern int nsp;
extern BackG *BG;	// background species
extern Species *SP;
extern Species *dev_info_sp;// particle species
extern GPG *dev_G_sp;
extern float dt;
extern float gputime;
extern float *phi_cond;
extern float **AM,*V_t,*b_t,*extq,*extq_1,*extq_2,*extq_3;
extern float *CondCharge;
extern float *Surf_charge,*Old_Surf_charge,*Old2_Surf_charge;
extern cudaEvent_t start, stop;
extern "C" void Main_Variable_printorSave();
extern "C" float frand();
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
float *dot_result2;
int *FIter;
GGA *dev_GvecSet;
GCA *dev_CvecSet;
float *dev_CondCharge;
float *dev_phi;     // PCG Solution
float *dev_phi_buf; // Sigma or buf
size_t pitch;
float *Lap_TEMP_Sol; // Temperature Solution of Laplace Equation
float *Lap_PHI_Sol; // Each of conductor Phi Solution of Laplace Equation, This is Device value
float **Lap_SIG_Sol; // Each of conductor Sigma Solution of Laplace Equation for external circuit
float *Pois_SIG_Sol; // Each of conductor Sigma Solution of Poisson Equation for external circuit
float *TotPotential;
float *LapPotential;
float *dev_Sigma, *dev_Source;
void Set_MatrixPCG_cuda();
void PCG_SOLVER_Laplace();
void PCG_SOLVER();
void Efield_cuda();
void Efield_cuda_Basic();
void cofactor(float **matrix, float size);
float determinant(float **matrix, float size);
void transpose(float **matrix, float matrix_cofactor[N_MAX][N_MAX],float size);
__device__ void Mat_x_Vec(int *I, int *J, float *val, int nnz, int num_rows, float alpha, float *inputVecX, 
                        float *outputVecY, cg::thread_block &cta, const cg::grid_group &grid);
__device__ void A_x_X_p_Y(float a, float *x, float *y, int size, const cg::grid_group &grid);
__device__ void A_x_Y_p_X(float a, float *x, float *y, int size, const cg::grid_group &grid);
__device__ void Vec_Dot_Sum(float *vecA, float *vecB, double *result, int size, const cg::thread_block &cta, const cg::grid_group &grid);
__device__ void CopyVector(float *srcA, float *destB, int size, const cg::grid_group &grid);
__device__ void Vec_x_Vec(float *vecA, float *vecB, float *vecC, int size, const cg::thread_block &cta, const cg::grid_group &grid);
__global__ void PCG(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, double *d_result);
__global__ void PCG_Deposit(int Gsize, int *IDX, GGA *vecG, float *X, float *PHI);
__global__ void PCG_Deposit_Lap(int Gsize, int *IDX, GGA *vecG, int k, float *X, float *PHI);
__global__ void PCG_Deposit_Temp(int Gsize, int *IDX, float *X, GGA *vecG);
__global__ void Cond_Sigma_Lap(int ngx, int ngy, float dx, float dy, float zlength, GGA *vecG, GCA *vecC, float *Phi, float *Sigma);
__global__ void Cond_Sigma(int ngx, int ngy, float hdx, float hdy, float idxy, float zlength, GGA *vecG, GCA *vecC, Species *info, GPG *data, float *Phi, float *Sigma);
__global__ void Cond_Sigma_v2(int ngx, int ngy, float hdx, float hdy, float idxy, float zlength, GGA *vecG, GCA *vecC, Species *info, GPG *data, float *Phi, float *Sigma);
__global__ void Calculate_1GasPara(int Gsize, float mass, float press, GGA *vecG);
__global__ void Calculate_2GasPara(int Gsize, float mass1, float press1, float mass2, float press2, GGA *vecG);
__global__ void SaveAT2D(float *A, size_t pitch, int height, float *PHI, int n);
__global__ void LoadAT2D(float *A, size_t pitch, int height, float *PHI, int n);
__global__ void VectorSum(int Gsize,float *TotPhi,float V,float *Phi);
__global__ void GGACopy_Potential(int Gsize, GGA *vecG, float *V1, float *V2);
__global__ void VtoEfield(int ngx,int ngy,float dx,float dy,float hdx,float hdy,float idx,
                            float idy, float *Sigma, float *Phi, GPG *data, GCA *vecC, GGA *vecG);
__device__ void Vec_Dot_Sum_F(float *vecA, float *vecB, float *result, int size, const cg::thread_block &cta, const cg::grid_group &grid);
__global__ void PCG_float(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, float *d_result);
#endif