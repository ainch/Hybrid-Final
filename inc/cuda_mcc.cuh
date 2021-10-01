#include "xypic.cuh"
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
namespace cg = cooperative_groups;
#ifndef __CUDA_MCC_CUH__
#define __CUDA_MCC_CUH__
CollF *dev_Coll_Flag;
ArCollD *dev_ArCX;
O2CollD *dev_O2CX;
ArO2CollD *dev_ArO2CX;
MCC_sigmav *Host_SigmaV;
MCC_sigmav *dev_SigmaV;
void MCC_Ar_cuda();
void Set_NullCollisionTime_cuda();
__global__ void MCC_Ar_cooper(int Gsize, int ngy, float dt, int MCCn, float dtm,float idx,float idy, int nvel, float *vsave,
											curandState *states, MCC_sigmav *sigv, CollF *CollP, ArCollD *CX, 
											Fluid *infoF, GFC *Fluid, GGA *BG, Species *info, GPG *data, GCP *sp);
__global__ void MCC_Ar_Basic(int Gsize, int ngy, curandState *states, Species *info, GCP *sp, GPG *data, GGA *Field);
#endif