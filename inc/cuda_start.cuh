#include "xypic.cuh"

extern int device_num;
extern int ConstB_Flag; // Magnetic field 
extern int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen
extern int Lap_Field_Solver_Test,Lap_Field_Solver_Flag,Lap_Field_Solver_Save;
extern void (*FieldSolver)();
extern void (*MOVE)();
extern void (*SORT_BOUNDARY)();
extern void (*MCC)();
extern void (*DEPOSIT)();
extern void (*CONTIEQ)();
extern void (*DIAG)();
#ifndef __CUDA_START_CUH__
#define __CUDA_START_CUH__
void info_Device(); 	// Select GPU
void start_cuda();		// Select Module
#endif