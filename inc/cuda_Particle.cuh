#include "xypic.cuh"
extern int Gsize;
extern int nsp;
extern int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen
extern int nsp, nfsp, nBG;
extern Species *SP;// particle species
extern Fluid *FG;	// fluid species
extern BackG *BG;	// background species
extern HCP *PtD;
extern GPG *Host_G_sp;
#ifndef __CUDA_PARTICLE_CUH__
#define __CUDA_PARTICLE_CUH__
GCP *Host_sp, *dev_sp;
GPG *dev_G_sp;
Species *dev_info_sp;
void Set_Particle_cuda();
void GCPInit(int size, GCP *A);
void Copy_HCPtoGCP(Species *info, HCP *A,GCP *B, GPG *C);
void Copy_GCPtoHCP(Species *info, GCP *A,HCP *B, GPG *C);
#endif