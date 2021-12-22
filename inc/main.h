#include "xypic.h"
extern int TecplotS_CX_Flag;
extern int TecplotS_Gsize_Flag;
extern int TecplotS_Particle_Flag;
extern int TecplotS_Particle_Num;
//
extern HCP *PtD;
extern void InputRead();
extern void InputFileMake();
extern void Source_setting();
extern void Geometry_setting();
extern void FieldSolverSetting();
extern void GasSetting();
extern void Main_Variable_printorSave();
extern void Cross_Section_data_Save();
extern void Initial_Particle_Save(int size,HCP *PtData);
extern void Fluid_Setting();
extern void Diagnostic_Setting();
extern void LoadDumpFile();
extern void main_cuda();

#ifndef __MAIN_H__
#define __MAIN_H__
long int seed; // related to Random number
double t;  // real time
int tstep; // number of time step
int cstep; // Number of Cycle step
int DumpFlag; // Dump File ON,OFF
char InputFile[80]; // INPUT FILE NAME
char DumpFile[80];  // DUMP FILE NAME
char *ConstBFile;  // DUMP FILE NAME
void display_title();
void display_finish();
void (*EFIELD)();
void (*MOVE)();
void (*SORT_BOUNDARY)();
void (*MCC)();
void (*MCC_Basic)();
void (*DEPOSIT)();
void (*CONTIEQ)();
void (*DIAG)();
#endif