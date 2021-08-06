#include "xypic.h"

extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern int DumpFlag; // Dump File ON,OFF
extern char InputFile[80]; // INPUT FILE NAME
extern char DumpFile[80];  // DUMP FILE NAME
extern char *ConstBFile;  // DUMP FILE NAME

extern void Argon_CrossSectionSET(CollF *CF);
extern void Oygen_CrossSectionSET(CollF *CF);
extern void ArO2_CrossSectionSET(CollF *CF);

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
int BoundaryNUM;
int *BoundaryX0,*BoundaryY0,*BoundaryX1,*BoundaryY1,*BoundaryBC;
float *BoundaryTEMP;
int CondNUM;
int*CondM_ID,*CondX0,*CondX1,*CondY0,*CondY1;
float *CondTEMP;
int SrcNUM;
int *SrcM_ID;
float *SrcDC, *SrcPOWER, *SrcAC, *SrcFREQ, *SrcPHASE, *SrcR, *SrcL, *SrcC;
float Min_FREQ, Max_FREQ;
int DielNUM;
int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
float *DielEPS;
HGA *vec_G;
HCA *vec_C;
int **StructureIndex;
int *vec_StructureIndex;
//
int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen
int nsp, nfsp, nBG;
Species *SP;// particle species
Fluid *FG;	// fluid species
BackG *BG;	// background species
float Total_Pressure;
//
int DT_PIC;
int DT_CONTI;
float dt;   // timestsep for PIC
float dtc; // time step for continuity equation
float PCGtol;
int HISTMAX;
int dHIST;
int np_lim; // total particle limit
int N_ave;
int N_smt;  // Number of smoothing every timestep
int ConstB_Flag; // Magnetic field 
//
int Basic_Flag; // 0 : Basic, 1: OTHERS
//
int nRct_cx,nRct_rc; // Number of reaction _ cross section or Reaction rate
int TnRct; // Total Number of reaction 
int mMnum;
int CX_TEC_Flag;
CollF *Coll_Flag;
//
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
int IVnZC(char A[50],int Value); // Int Value non Zero CHECK
float FIVnZC(char A[50],float Value);// float Value non Zero CHECK
#endif