#include "xypic.cuh"
extern char InputFile[80]; // INPUT FILE NAME
extern int TecplotS_2D_Flag;
extern int TecplotS_2D_Ncycle;
extern int TecplotS_Movie_Flag;
extern int TecplotS_Movie_Ncycle;
extern int TecplotS_Movie_Frame;
extern int TecplotS_Movie_SCYCLE;
extern int TecplotS_Movie_Count;
extern int TecplotS_PT_Movie_Flag;
extern int TecplotS_PT_Movie_Ncycle;
extern int TecplotS_PT_Movie_Frame;
extern int TecplotS_PT_Movie_SCYCLE;
extern int TecplotS_PT_Movie_Count;
extern int Total_maxnp;
extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern float dt;   // timestsep for PIC
extern int CYCLE_NUM; // Minimum frequency number of cycle
extern int ngx,ngy,Gsize;
extern int nsp;
extern float dx,dy;
extern float xlength,ylength,zlength;
extern float *x_Garray,*y_Garray;
extern int BoundaryNUM;
extern int *BoundaryX0,*BoundaryY0,*BoundaryX1,*BoundaryY1;
extern int CondNUM,CondNUMR;
extern int*CondM_ID, *CondX0,*CondX1,*CondY0,*CondY1;
extern int DielNUM,DielNUMR;
extern int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
extern Species *SP;
extern Species *dev_info_sp;
extern GGA *vec_G;
extern GGA *dev_GvecSet;
extern GPG *Host_G_sp;
extern GPG *dev_G_sp;
extern HCP *PtD;
extern GCP *Host_sp, *dev_sp;
extern void Copy_GCPtoHCP(Species *info, GCP *A,HCP *B, GPG *C);
#ifndef __CUDA_TECPLOT_CUH__
#define __CUDA_TECPLOT_CUH__
int PT_Movie_S_count;
void Tecplot_save();
void Tecplot_Gsize_Movie(int Init);
void Tecplot_PT_Movie(int Init,int isp);
void Tecplot_2D();
#endif