#include "xypic.cuh"

extern int device_num;
extern int ngx,ngy,Gsize;
extern int A_size;
extern float PCGtol;
extern float PCGtol2;
extern HGA *vec_G;
extern int CondNUMR;
extern int Lap_Field_Solver_Test,Lap_Field_Solver_Flag,Lap_Field_Solver_Save;
extern int Preconditioner_Flag;
extern int FieldIter;
extern int A_size;
extern float *A_val,*MatTA;
extern int *Ai,*Aj;
extern int **A_idx;
extern float *MatM,**cond_b,*temp_b;
extern float **phi_dw,**phi_u;
extern int FIELD_GRID, FIELD_BLOCK;
extern "C" void Field_Laplace_Solution_Save(char *Filename,float **Sol);
extern "C" float *VFMalloc(int size);
extern "C" void VFInit(float *V,float C,int size);
extern "C" void VFCopy(float *V,float *C,int size);

#ifndef __CUDA_FIELD_CUH__
#define __CUDA_FIELD_CUH__

void Set_MatrixPCG_cuda();
void PCG_SOLVER_Laplace();

#endif