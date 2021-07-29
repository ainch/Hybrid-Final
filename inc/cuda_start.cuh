#include "xypic.cuh"

extern int device_num;

#ifndef __CUDA_START_H__
#define __CUDA_START_H__
void info_Device(); 	// Select GPU
void start_cuda();		// Select Module
#endif