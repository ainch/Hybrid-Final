#include "xypic.h"

extern void InputRead(int argc, char *argv[]);
extern void start();
extern void DumpRead(int argc, char *argv[]);
extern void main_cuda();

#ifndef __MAIN_H__
#define __MAIN_H__
double t;  // real time
int tstep; // number of time step
int cstep; // Number of Cycle step
int DumpFlag; // Dump File ON,OFF
char InputFile[80]; // INPUT FILE NAME
char DumpFile[80];  // DUMP FILE NAME
char *ConstBFile;  // DUMP FILE NAME

void display_title();
void display_finish();
#endif