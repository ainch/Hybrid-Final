#include "xypic.cuh"
extern int nfsp;
extern int Csize;
extern Fluid *FG;	// fluid species
extern GFC *Host_C_F;
#ifndef __CUDA_FLUID_CUH__
#define __CUDA_FLUID_CUH__
Fluid *dev_FG;
GFC *dev_C_F;
void Set_Fluid_cuda();
#endif