#include "xypic.h"

#ifndef __CPUvalue_H__
#define __CPUvalue_H__
// C Variable declaration
long int seed; // related to Random number
double t;  // real time
int tstep; // number of time step
int cstep; // Number of Cycle step
int DumpFlag; // Dump File ON,OFF
char InputFile[80]; // INPUT FILE NAME
char DumpFile[80];  // DUMP FILE NAME
int device_num;
float xlength,ylength,zlength;
int ngx,ngy,Gsize;
int ncx,ncy,Csize;
float dx,dy;
float idx,idy,dx2,dy2,dxdy2,hdx,hdy,r_eps0;
float fncx,fncy,fngx,fngy;
float *x_Garray,*y_Garray;
int BoundaryNUN;
int *BoundaryX0,*BoundaryY0,*BoundaryX1,*BoundaryY1,*BoundaryBC;
float *BoundaryTEMP;
int CondNUN;
int*CondM_ID,*CondX0,*CondX1,*CondY0,*CondY1;
float *CondTEMP;
int SrcNUN;
int *SrcM_ID;
float *SrcDC, *SrcPOWER, *SrcAC, *SrcFREQ, *SrcPHASE, *SrcR, *SrcL, *SrcC;
float Min_FREQ;
int DielNUN;
int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
float *DielEPS;
int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen
#endif