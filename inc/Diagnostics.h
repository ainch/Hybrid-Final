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
extern int Hcount;
extern int hist_count, hist_ave_count;
extern int nave_count;
extern float *t_array;
extern float *iter_array;    
extern Hist *HistPt,*HistFG;
extern float *t_ave_array;
extern Hist *Hist_ave_Pt,*Hist_ave_Pt_stack;
extern float **Current_hist,**Surf_charge_hist,**Volt_hist,**Volt_cond_hist;
extern float ***SP_current_hist;
extern float *Current_Now;
//
extern float *vec_Potential, *sum_Potential, *ave_Potential; // [Gsize] potential
extern float *vec_Source, *sum_Source, *ave_Source; // [Gsize] Charge density 
extern float *vec_Sigma, *sum_Sigma, *ave_Sigma; // [Gsize] [Dielectric] surface charge density, [Conductor] Surface current
extern float *sum_Ex, *ave_Ex, *sum_Ey, *ave_Ey;
//
extern int *Stack_ave_np, *ave_np;
extern float *new_ave_np, *old_ave_np;
void Diagnostic_Setting();
#endif

