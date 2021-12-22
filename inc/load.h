#include "xypic.h"
extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern char DumpFile[80];  // DUMP FILE NAME
extern int Gsize,ngx,ngy;
extern int Csize,ncx,ncy;
extern float dx,dy;
extern float xlength,ylength,zlength;
extern float *x_Garray,*y_Garray;
extern int CondNUM, CondNUMR;
extern int *CondM_ID,*CondX0,*CondX1,*CondY0,*CondY1;
extern int DielNUM;
extern int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
extern GGA *vec_G;
extern GCA *vec_C;
extern GPG *Host_G_sp;
extern GFC *Fluid_sp;
extern GFG *Fluid_Den, *dev_FG_Den;
extern int **A_idx;
extern int MainGas;
extern int nsp;
extern Species *SP;// particle species
extern HCP *PtD;
extern int HISTMAX;
extern int dHIST;
extern int Hcount;
extern int hist_count, hist_ave_count;
extern Hist *HistPt,*HistFG;
extern float **Current_hist,**Surf_charge_hist,**Volt_hist,**Volt_cond_hist;
extern float ***SP_current_hist;
extern float *t_array;  
extern float *iter_array;  
extern float *t_ave_array;
extern Hist *Hist_ave_Pt,*Hist_ave_Pt_stack;
extern int nfsp;
extern int Conti_Flag;
extern GFC *Fluid_sp;
extern int DumpFlag;
extern float *MCC_rate, *ave_MCC_rate;
extern float *vec_Potential, *ave_Potential; // [Gsize] potential
extern float *vec_Source, *ave_Source; // [Gsize] Charge density 
extern float *vec_Sigma, *ave_Sigma; // [Gsize] [Dielectric] surface charge density, [Conductor] Surface current
extern float *ave_Ex, *ave_Ey;
extern int TnRct;
extern int **MIMalloc(int sizeX,int sizeY);
extern void MIInit(int **M,int C,int sizeX,int sizeY);
#ifndef __START_H__
#define __START_H__
void V000_LoadDUMP(FILE *LF);
void V001_LoadDUMP(FILE *LF);
void LoadDumpFile();
void CG_Matrix_Setting(float *A, int *Ai, int *Aj, float **b, float *M, float *Atemp, float *btemp);
#endif