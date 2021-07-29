#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "def.h"
#include <omp.h> //OpenMP
//CUDA
#include "cuda.h"
#include "cublas_v2.h"
#include "cusparse_v2.h"
#include "cusolverDn.h"
#include "cusolverSp.h"
#include "cusolverRf.h"
#include "cusolver_common.h"
#include "curand_kernel.h"
#include "cuda_runtime.h"
#include "cufft.h"
#include "helper_cuda.h"
#include "interop.cuh"

#ifndef __CUDS_XYPIC__H
#define __CUDS_XYPIC__H
// extern C Function declaration
extern "C" int ***TIMalloc(int sizeX,int sizeY,int sizeZ);
extern "C" float ***TFMalloc(int sizeX,int sizeY,int sizeZ);
extern "C" int **MIMalloc(int sizeX,int sizeY);
extern "C" float **MFMalloc(int sizeX,int sizeY);
extern "C" int *VIMalloc(int size);
extern "C" float *VFMalloc(int size);
extern "C" void TFFree(float ***T,int sizeX,int sizeY);
extern "C" void TIFree(int ***T,int sizeX,int sizeY);
extern "C" void MFFree(float **M,int sizeX);
extern "C" void MIFree(int **M,int sizeX);
extern "C" void TFInit(float ***T,float C,int sizeX,int sizeY,int sizeZ);
extern "C" void TIInit(int ***T,int C,int sizeX,int sizeY,int sizeZ);
extern "C" void MFInit(float **M,float C,int sizeX,int sizeY);
extern "C" void MIInit(int **M,int C,int sizeX,int sizeY);
extern "C" void VFInit(float *V,float C,int size);
extern "C" void VIInit(int *V,int C,int size);
extern "C" void MFCopy(float **M,float **C,int sizeX,int sizeY);
extern "C" void MICopy(int **M,int **C,int sizeX,int sizeY);
extern "C" void VFCopy(float *V,float *C,int size);
extern "C" void VICopy(int *V,int *C,int size);
// CU Function declaration

#endif