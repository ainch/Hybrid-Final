#include "xypic.cuh"
extern int ngx,ngy;
extern void info_Device(); 					// Select GPU
extern void start_cuda();					// Select Module

#ifndef __CUDA_MAIN_H__
#define __CUDA_MAIN_H__

#endif

#ifndef __TEST__
#define __TEST__
int test(void);
__global__ void testKernel(point *p);
__global__ void MakeVectorForMoveKernel(int ngx,int ngy,point *p);
#endif

