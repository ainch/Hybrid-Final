#include "xypic.h"

extern int ngx,ngy;
extern int ncx,ncy;
extern float dx,dy;
extern float xlength,ylength,zlength;
extern float *x_Garray,*y_Garray;
extern int CondNUM;
extern int *CondM_ID,*CondX0,*CondX1,*CondY0,*CondY1;
extern int DielNUM;
extern int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
extern HGA *vec_G;
extern HCA *vec_C;
extern int **A_idx;
extern int **MIMalloc(int sizeX,int sizeY);
extern void MIInit(int **M,int C,int sizeX,int sizeY);
#ifndef __START_H__
#define __START_H__
void CG_Matrix_Setting(float *A, int *Ai, int *Aj, float **b, float *M, float *Atemp, float *btemp);
#endif