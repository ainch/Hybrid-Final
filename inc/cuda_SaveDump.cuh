#include "xypic.cuh"
extern char InputFile[80]; // INPUT FILE NAME
extern char DumpFile[80]; // DUMP FILE NAME
extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern int DT_PIC;  // Number of 1 cycle step
extern int init_dump_num;
extern int OVER_dump_order;
extern int dump_order;
extern int dump_num;
extern float *dump_cycle;
extern char InputFile[80]; // INPUT FILE NAME
extern int Total_maxnp;
extern Species *SP; // particle species
extern Species *dev_info_sp;// particle species
extern HCP *PtD;
extern GCP *Host_sp, *dev_sp;
extern GPG *Host_G_sp, *dev_G_sp;
extern GFC *Fluid_sp;
extern GFG *Fluid_Den, *dev_FG_Den;
extern GGA *vec_G, *dev_GvecSet;
extern float *MCC_rate, *ave_MCC_rate, *dev_ave_MCC_rate;
extern float *RCstack, *dev_RCstack;
extern float *TotPotential;
extern float *vec_Potential, *ave_Potential; // [Gsize] potential
extern float *vec_Source, *ave_Source; // [Gsize] Charge density 
extern float *vec_Sigma, *ave_Sigma; // [Gsize] [Dielectric] surface charge density, [Conductor] Surface current
extern float *ave_Ex, *ave_Ey;
extern float *dev_ave_Potential;
extern float *dev_ave_Source;
extern float *dev_ave_Sigma;
extern float *dev_ave_Ex;
extern float *dev_ave_Ey;
extern int MainGas;
extern int nsp,nfsp;
extern int Gsize,ngx,ngy;
extern int Csize,ncx,ncy;
extern int CondNUMR;
extern int TnRct;
extern int HISTMAX;
extern int dHIST;
extern int Hcount;
extern int hist_count, hist_ave_count;
extern Hist *HistPt,*HistFG;
extern float *t_ave_array;
extern Hist *Hist_ave_Pt,*Hist_ave_Pt_stack;
extern float **Current_hist,**Surf_charge_hist,**Volt_hist,**Volt_cond_hist;
extern float ***SP_current_hist;
extern float *t_array;  
extern float *iter_array;  
extern int Basic_Flag;
extern float *Host_G_buf, *Host_C_buf;
extern float 	totaltime,gputime_field,gputime_efield;
extern float 	gputime_move,gputime_mcc,gputime_deposit;
extern float 	gputime_diag,gputime_sort,gputime_Tec;
extern float 	gputime_continue,gputime_dump;
extern int		TotalT_D;
extern int		TotalT_H;
extern int		TotalT_M;
extern int		TotalT_S;
extern float EndTime;
extern float Margin_ave_np;
extern int Flag_ave_np, Same_ave_np, Stack_ave_np;
extern int Basic_Flag; // 0 : Basic, 1: OTHERS
extern int Conti_Flag;
extern GFC *Fluid_sp;
extern void Copy_GCPtoHCP(Species *info, GCP *A, HCP *B, GPG *C);
#ifndef __CUDA_SAVEDUMP_CUH__
#define __CUDA_SAVEDUMP_CUH__
void SaveDumpFile(int KEY2,int KEY1,int KEY0);
void V000_DUMP(FILE *SF);
void V001_DUMP(FILE *SF);
#endif