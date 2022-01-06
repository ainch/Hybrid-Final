#include "xypic.cuh"
extern int Gsize;
extern int device_num;
extern int ConstB_Flag; // Magnetic field 
extern int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen
extern int Lap_Field_Solver_Test,Lap_Field_Solver_Flag,Lap_Field_Solver_Save;
extern int nsp;
extern Species *SP;// particle species
extern int External_Flag; // 0 : Voltage driven, 1: Power driven
extern int CSS_Flag;
extern void (*EFIELD)();
extern void Efield_cuda();
extern void Efield_cuda_Basic();
extern void (*MOVE)();
extern void Move_cuda();
extern void (*SORT_BOUNDARY)();
extern void SortBounndary_cuda();
extern void (*MCC)();
extern void (*MCC_Basic)();
extern void MCC_Ar_cuda();
extern void MCC_O2_cuda();
extern void MCC_ArO2_cuda();
extern void (*DEPOSIT)();
extern void Deposit_cuda();
extern void Deposit_Basic();
extern void (*CONTIEQ)();
extern void Solve_Continuity_eqn();
extern void Solve_Continuity_eqn_check();
extern void (*DIAG)();
extern void Diagnostic();
extern void Diagnostic_Basic();
#ifndef __CUDA_START_CUH__
#define __CUDA_START_CUH__
extern cudaDeviceProp prop;
void info_Device(); 	// Select GPU
void start_cuda();		// Select Module
#endif