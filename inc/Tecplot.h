#include "xypic.h"

extern float xlength,ylength,zlength;
extern int ngx,ngy,Gsize;
extern float dx,dy;
extern float *x_Garray,*y_Garray;
extern int BoundaryNUM;
extern int *BoundaryX0,*BoundaryY0,*BoundaryX1,*BoundaryY1;
extern int CondNUM,CondNUMR;
extern int*CondM_ID, *CondX0,*CondX1,*CondY0,*CondY1;
extern int DielNUM,DielNUMR;
extern int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
extern HGA *vec_G;
extern int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen
extern float LOGX_MIN,LOGX_MAX;
extern float dLOGX,idLOGX;
extern int N_LOGX;
extern ArCollD *Ar_Data;
extern O2CollD *O2_Data;
extern ArO2CollD *ArO2_Data;
#ifndef __Tecplot_H__
#define __Tecplot_H__
void Cross_Section_data_Save();
void Main_Variable_printorSave();
#endif