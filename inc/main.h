#include "xypic.h"
// C Variable declaration
long int seed; // related to Random number
double t;  // real time
int tstep; // number of time step
int cstep; // Number of Cycle step
int DumpFlag; // Dump File ON,OFF
char InputFile[80]; // INPUT FILE NAME
char DumpFile[80];  // DUMP FILE NAME
int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen

