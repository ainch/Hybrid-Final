#include "xypic.cuh"

extern long int seed;
extern int ngx,ngy;
extern int Gsize;
extern int device_num;
#ifndef __CUDA_INIT_CUH__
#define __CUDA_INIT_CUH__
curandState *devStates;
int FIELD_GRID, REDUCTION_GRID, EFIELD_GRID, MOVE_GRID, SORT_GRID, MCC_GRID,DEPOSIT_GRID;
int FIELD_BLOCK, REDUCTION_BLOCK, EFIELD_BLOCK, MOVE_BLOCK, SORT_BLOCK,MCC_BLOCK, DEPOSIT_BLOCK;

void Set_Device_Parameter();
void Set_Particle_cuda();
void Set_NullCollisionTime_cuda();
void Set_DiagParameter_cuda();
__global__ void SetSeed(curandState *state,long int seed,int num);
__global__ void MyKernel(int *array, int arrayCount);
#endif