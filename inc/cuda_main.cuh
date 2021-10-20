#include "xypic.cuh"
extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern float dt;   // timestsep for PIC
extern int CYCLE_NUM; // Minimum frequency number of cycle
extern int ngx,ngy,Gsize;
extern int MainGas;
extern int TnRct; // Total Number of reaction 
extern int Msize;
extern float *MCC_rate;
extern float *dev_MCC_rate;
extern float dx,dy;
extern int nsp;
extern float PCGtol2;
extern Species *SP;// particle species
extern GPG *Host_G_sp;
extern GPG *dev_G_sp;
extern GGA *vec_G;
extern int Lap_Field_Solver_Test,Lap_Field_Solver_Flag,Lap_Field_Solver_Save;
extern double *dot_result;
extern int *FIter;
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
extern void Set_Fluid_cuda();
extern void Set_NullCollisionTime_cuda();
extern void Set_DiagParameter_cuda();
extern void Set_SortBoundary_cuda();
extern void PCG_Laplace_TEST();
extern void PCG_SOLVER_Laplace();
extern void Deposit_cuda();
extern void PCG_SOLVER();
extern void Efield_cuda();
extern void Move_Sort_cuda();
extern void Move_cuda();
extern void SortBounndary_cuda();
extern void MCC_Ar_cuda();
extern void MCC_O2_cuda();
extern void MCC_ArO2_cuda();
extern void Tecplot_save();
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