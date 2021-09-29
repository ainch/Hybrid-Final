#include "xypic.cuh"
extern float dt;
extern int nsp, nfsp, nBG;
extern int TnRct; // Total Number of reaction 
extern int MainGas;
extern int DT_MCC; // mcc count for each step
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
#ifndef __CUDA_MCC_CUH__
#define __CUDA_MCC_CUH__
CollF *dev_Coll_Flag;
ArCollD *dev_ArCX;
O2CollD *dev_O2CX;
ArO2CollD *dev_ArO2CX;
MCC_sigmav *Host_SigmaV;
MCC_sigmav *dev_SigmaV;
void Set_NullCollisionTime_cuda();
#endif