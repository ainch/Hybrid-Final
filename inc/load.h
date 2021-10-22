#include "xypic.h"
extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern char DumpFile[80];  // DUMP FILE NAME
extern int ngx,ngy;
extern int ncx,ncy;
extern float dx,dy;
extern float xlength,ylength,zlength;
extern float *x_Garray,*y_Garray;
extern int CondNUM;
extern int *CondM_ID,*CondX0,*CondX1,*CondY0,*CondY1;
extern int DielNUM;
extern int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
extern GGA *vec_G;
extern GCA *vec_C;
extern int **A_idx;
extern int MainGas;
extern int nsp;
extern Species *SP;// particle species
extern HCP *PtD;
extern int **MIMalloc(int sizeX,int sizeY);
extern void MIInit(int **M,int C,int sizeX,int sizeY);
#ifndef __START_H__
#define __START_H__
void V000_LoadDUMP(FILE *LF);
void V001_LoadDUMP(FILE *LF);
void LoadDumpFile();
void CG_Matrix_Setting(float *A, int *Ai, int *Aj, float **b, float *M, float *Atemp, float *btemp);
#endif