#include "xypic.h"
extern float dt;
extern int nsp,nfsp;
extern int HISTMAX;
extern int dHIST;
extern int Gsize;
extern int CondNUMR;
extern Species *SP;// particle species
#ifndef __DIAGNOSTICS_H__
#define __DIAGNOSTICS_H__
int Hcount;
int hist_count;
int nave_count;
float *t_array;
float *iter_array;    
Hist *HistPt,*HistFG;
float **Current_hist,**Surf_charge_hist,**Volt_hist,**Volt_cond_hist;
float ***SP_current_hist;
float *Current_Now;
//
float *vec_Potential, *sum_Potential, *ave_Potential; // [Gsize] potential
float *vec_Source, *sum_Source, *ave_Source; // [Gsize] Charge density 
float *vec_Sigma, *sum_Sigma, *ave_Sigma; // [Gsize] [Dielectric] surface charge density, [Conductor] Surface current
float *sum_Ex, *ave_Ex, *sum_Ey, *ave_Ey;
//
int *Stack_ave_np, *ave_np;
float *new_ave_np, *old_ave_np;
void Diagnostic_Setting();
#endif

