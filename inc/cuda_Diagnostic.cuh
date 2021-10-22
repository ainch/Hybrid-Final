#include "xypic.cuh"
extern int Gsize,Csize;
#ifndef __CUDA_DIAGNOSTIC_CUH__
#define __CUDA_DIAGNOSTIC_CUH__
float *Host_G_buf, *Host_C_buf;
void Diagnostic();
void Set_Diagnostic_cuda();
#endif