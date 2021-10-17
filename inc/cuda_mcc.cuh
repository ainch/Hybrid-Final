#include "xypic.cuh"
extern long int seed;
extern float dt;
extern int Gsize,ngx,ngy;
extern int Csize,ncx,ncy;
extern float idx,idy;
extern int nsp, nfsp, nBG;
extern int TnRct; // Total Number of reaction 
extern int MainGas;
extern int DT_MCCn; // mcc count for each step
extern float dt_mcc; // timestsep for MCC Module
extern Species *SP;// particle species
extern Fluid *FG;	// fluid species
extern BackG *BG;	// background species
extern float LOGX_MIN,LOGX_MAX;
extern float dLOGX,idLOGX;
extern int N_LOGX;
extern CollF *Coll_Flag;
extern ArCollD *Ar_Data;
extern O2CollD *O2_Data;
extern ArO2CollD *ArO2_Data;
extern int sMemSize_MCC;
extern dim3 MCC_GRID, MCC_BLOCK;
extern dim3 MCC_GRID2, MCC_BLOCK2;
extern Species *dev_info_sp;// particle species
extern GCP *dev_sp;
extern GPG *dev_G_sp;
extern GGA *dev_GvecSet;
extern GFC *dev_C_F;
extern Fluid *dev_FG;
extern curandState *devStates;
extern int h_nvel;
extern float *dev_vsave;
extern int N_LOGX;
extern float idLOGX;
namespace cg = cooperative_groups;
extern __device__ void Direct_Argon_Electron(int Gsize, int ngy, int ID, int MCCn, float dtm, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, GGA *BG, GFC *Fluid);
extern __device__ void Direct_Argon_ArIon(int Gsize, int ngy, int ID, int MCCn, float dt, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, GGA *BG, GFC *Fluid);
extern __device__ void ArCollision_Check(int Gsize, int Csize, int ngy, int ID, int isp, float dt, int MCCn, float dtm, curandState *states, 
											Species *info, GPG *data, GCP *sp, MCC_sigmav *sigv, GGA *BG, GFC *Fluid);
extern __device__ void Argon_E_Collision(int Gsize, int ngy, int TID, int MCCn, float dtm, int nvel,float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX,
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, GGA *BG, GFC *Fluid);
extern __device__ void Argon_Ar_Collision(int Gsize, int TID, float dt, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, GGA *BG, GFC *Fluid);
extern __device__ float Argon_CrossSection(int R, float engy, int N_LOGX, float idLOGX, ArCollD *data);
extern __device__ void  O2Collision_Check(int Gsize, int Csize, int ngy, int TID, float dt, int MCCn, float dtm, float dx, float dy,
                                        curandState *states, Species *info, GPG *data, GCP *sp, MCC_sigmav *sigv, GGA *BG, GFC *Fluid);
extern __device__ void O2_Electron(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, O2CollD *CX, GGA *BG);
extern __device__ void O2_O_negative(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, O2CollD *CX, GGA *BG);
extern __device__ void O2_O2_ion(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, O2CollD *CX, GGA *BG);
extern __device__ void O2_O_ion(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, O2CollD *CX, GGA *BG);
#ifndef __CUDA_MCC_CUH__
#define __CUDA_MCC_CUH__
CollF *dev_Coll_Flag;
ArCollD *dev_ArCX;
O2CollD *dev_O2CX;
ArO2CollD *dev_ArO2CX;
MCC_sigmav *Host_SigmaV;
MCC_sigmav *dev_SigmaV;
void MCC_Ar_cuda();
void MCC_O2_cuda();
void MCC_ArO2_cuda();
void Set_NullCollisionTime_cuda();
__global__ void MCC_ArO2_Basic(int Gsize, int Csize, int ngy, int nsp, float dt, int MCCn, float dtm,float idx,float idy, int nvel, float *vsave,
											curandState *states, int N_LOGX, float idLOGX, MCC_sigmav *sigv, CollF *CollP, ArCollD *CX, 
											Fluid *infoF, GFC *Fluid, GGA *BG, Species *info, GPG *data, GCP *sp);
__global__ void MCC_O2_Basic(int Gsize, int Csize, int ngy, int nsp, float dt, int MCCn, float dtm,float idx,float idy, int nvel, float *vsave,
											curandState *states, int N_LOGX, float idLOGX, MCC_sigmav *sigv, CollF *CollP, O2CollD *CX, 
											Fluid *infoF, GFC *Fluid, GGA *BG, Species *info, GPG *data, GCP *sp);
__global__ void MCC_Ar_Basic(int Gsize, int Csize, int ngy, int nsp, float dt, int MCCn, float dtm,float idx,float idy, int nvel, float *vsave,
											curandState *states, int N_LOGX, float idLOGX, MCC_sigmav *sigv, CollF *CollP, ArCollD *CX, 
											Fluid *infoF, GFC *Fluid, GGA *BG, Species *info, GPG *data, GCP *sp);
__device__ void dev_maxwellv(float *vx_local,float *vy_local,float *vz_local,float vsaven,float vti,float Rphi,float Rthe);
__device__ void dev_newvel_IONSC(float *vx_sc,float *vy_sc,float *vz_sc,float vel,float rand1,float rand2);		
__device__ void dev_anewvel(float energy,float vel,float* n_vx,float* n_vy,float* n_vz,int e_flag,float massrate,float rand1,float rand2);
#endif