#include "xypic.h"

extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern int DumpFlag; // Dump File ON,OFF
extern char InputFile[80]; // INPUT FILE NAME
extern char DumpFile[80];  // DUMP FILE NAME

#ifndef __START_H__
#define __START_H__
// C Variable declaration
long int seed; // related to Random number
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
int nsp, nfsp, nBG;
Species *SP;// particle species
Fluid *FG;	// fluid species
BackG *BG;	// background species
float Total_Pressure;

void InputRead(int argc, char *argv[]);
void start();
void DumpRead(int argc, char *argv[]);
int ***TIMalloc(int sizeX,int sizeY,int sizeZ);
float ***TFMalloc(int sizeX,int sizeY,int sizeZ);
int **MIMalloc(int sizeX,int sizeY);
float **MFMalloc(int sizeX,int sizeY);
int *VIMalloc(int size);
float *VFMalloc(int size);
void TFFree(float ***T,int sizeX,int sizeY);
void TIFree(int ***T,int sizeX,int sizeY);
void MFFree(float **M,int sizeX);
void MIFree(int **M,int sizeX);
void TFInit(float ***T,float C,int sizeX,int sizeY,int sizeZ);
void TIInit(int ***T,int C,int sizeX,int sizeY,int sizeZ);
void MFInit(float **M,float C,int sizeX,int sizeY);
void MIInit(int **M,int C,int sizeX,int sizeY);
void VFInit(float *V,float C,int size);
void VIInit(int *V,int C,int size);
void MFCopy(float **M,float **C,int sizeX,int sizeY);
void MICopy(int **M,int **C,int sizeX,int sizeY);
void VFCopy(float *V,float *C,int size);
void VICopy(int *V,int *C,int size);
#endif