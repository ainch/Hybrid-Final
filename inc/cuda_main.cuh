#include "xypic.cuh"
extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern float dt;   // timestsep for PIC
extern int CYCLE_NUM; // Minimum frequency number of cycle
extern int ngx,ngy,Gsize;
extern int nsp;
extern Species *SP;// particle species
extern int Lap_Field_Solver_Test,Lap_Field_Solver_Flag,Lap_Field_Solver_Save;
extern void (*FieldSolver)();
extern void (*MOVE)();
extern void (*SORT_BOUNDARY)();
extern void (*MCC)();
extern void (*DEPOSIT)();
extern void (*CONTIEQ)();
extern void (*DIAG)();
// In cuda_*.cu
extern void info_Device(); 					// Select GPU
extern void start_cuda();					// Select Module
extern void Set_Device_Parameter();
extern void Set_MatrixPCG_cuda();
extern void Set_Particle_cuda();
extern void Set_NullCollisionTime_cuda();
extern void Set_DiagParameter_cuda();
extern void PCG_Laplace_TEST();
extern void PCG_SOLVER_Laplace();

#ifndef __CUDA_MAIN_CUH__
#define __CUDA_MAIN_CUH__
float time_sum;
float gputime;
float 	totaltime,gputime_field,gputime_efield;
float 	gputime_move,gputime_mcc,gputime_deposit;
float 	gputime_diag,gputime_sort,gputime_trace;
float 	gputime_continue,gputime_dump;
int		TotalT_D;
int		TotalT_H;
int		TotalT_M;
int		TotalT_S;
#endif

#ifndef __TEST__
#define __TEST__
int test(void);
__global__ void testKernel(point *p);
__global__ void MakeVectorForMoveKernel(int ngx,int ngy,point *p);
#endif

