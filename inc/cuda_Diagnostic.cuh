#include "xypic.cuh"
extern float dt;
extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern int DT_PIC;
extern int N_ave;
extern int Gsize,Csize;
extern int ngx,ngy,ncx,ncy;
extern int Conti_Flag;
extern int MainGas;
extern int nsp,nfsp;
extern int TnRct;
extern int CondNUMR;
extern int *SrcM_ID;
extern float *SrcDC, *SrcPOWER, *SrcAC, *SrcFREQ, *Src2piFREQ, *SrcPHASE, *SrcRPHASE;
extern float *dev_CondCharge;
extern float *CondCharge;
extern Species *SP;// particle species
extern Species *dev_info_sp;// particle species
extern GPG *dev_G_sp;
extern int HISTMAX;
extern int dHIST;
extern int Hcount;
extern int hist_count, hist_ave_count;
extern Hist *HistPt,*HistFG;
extern float *t_ave_array;
extern Hist *Hist_ave_Pt,*Hist_ave_Pt_stack;
extern float **Current_hist,**Surf_charge_hist,**Volt_hist,**Volt_cond_hist;
extern float ***SP_current_hist;
extern float *Current_Now;
extern float *t_array;  
extern float *iter_array;  
extern int *FIter;
extern int nave_count;
extern int PD_intv;
extern float PD_Ratio;
extern float *phi_cond;
extern float **AM,*V_t,*b_t,*extq,*extq_1,*extq_2,*extq_3;
extern float *CondCharge;
extern float *Surf_charge,*Old_Surf_charge,*Old2_Surf_charge;
extern int *Efield_Flag, *Cond_Source_num, *Cond_count,**Cond_Power_ID;
extern float *Cond_Power;
extern dim3 FIELD_GRID,FIELD_BLOCK;
extern dim3 FIELD_GRID2,FIELD_BLOCK2;
extern dim3 DEPOSIT_GRID,DEPOSIT_BLOCK;
extern dim3 EFIELD_GRID,EFIELD_BLOCK;
extern dim3 MOVE_GRID, MOVE_BLOCK;
extern dim3 SORT_GRID, SORT_BLOCK;
extern dim3 MCC_GRID, MCC_BLOCK;
extern dim3 CONTI_GRID,CONTI_BLOCK;
extern dim3 DIAG_G_GRID, DIAG_G_BLOCK;
extern dim3 DIAG_NSPG_GRID, DIAG_NSPG_BLOCK;
extern float *TotPotential;
extern float *dev_Source;
extern float *dev_Sigma;
extern float *vec_Potential, *ave_Potential; // [Gsize] potential
extern float *vec_Source, *ave_Source; // [Gsize] Charge density 
extern float *vec_Sigma, *ave_Sigma; // [Gsize] [Dielectric] surface charge density, [Conductor] Surface current
extern float *ave_Ex, *ave_Ey;
extern float *MCC_rate, *ave_MCC_rate;
extern float *dev_MCC_rate,*dev_ave_MCC_rate;
extern GGA *vec_G,*dev_GvecSet;
extern GCA *vec_C,*dev_CvecSet;
extern float Margin_ave_np;
extern int Flag_ave_np, Same_ave_np;
extern int Basic_Flag; // 0 : Basic, 1: OTHERS
extern int *Stack_ave_np, *ave_np;
extern float *new_ave_np, *old_ave_np;
extern CollF *dev_Coll_Flag;
extern Fluid *FG;
extern GFG *dev_FG_Den, *dev_FG_Src;
extern GFG *Fluid_Den, *Fluid_Src;
extern GFC *Fluid_sp;
extern int CSS_Flag;
extern void Sync_Fluid_GFCtoGFG_forDen(GFC *A, GFG *B);
extern void Sync_Fluid_GFGtoGFC_forSource(GFG *A, GFC *B);
#ifndef __CUDA_DIAGNOSTIC_CUH__
#define __CUDA_DIAGNOSTIC_CUH__
extern float *Host_G_buf, *Host_C_buf;
extern float *dev_sum_Potential, *dev_ave_Potential;
extern float *dev_sum_Source, *dev_ave_Source;
extern float *dev_sum_Sigma, *dev_ave_Sigma;
extern float *dev_sum_Ex, *dev_ave_Ex;
extern float *dev_sum_Ey, *dev_ave_Ey;
extern float *dev_phi_buf;
void Diagnostic();
void Diagnostic_Basic();
void Set_Diagnostic_cuda();
__global__ void SumReductionINT1024All(int n, int isp, GPG *g_data, int *g_max);
__device__ void warpSumReduceINT(volatile int* sdata,int TID);
__global__ void Accomulate_Field_Data(int Gsize, float *TotPot, float *Source, float *Sigma, GGA *vecG
                        ,float *sum_Potential, float *sum_Source, float *sum_Sigma, float *sum_Ex, float *sum_Ey);
__global__ void Average_Field_Data(int Gsize, int N_ave, float *TotPot, float *Source, float *Sigma, GGA *vecG
                        ,float *sum_Potential, float *sum_Source, float *sum_Sigma, float *sum_Ex, float *sum_Ey
                        ,float *ave_Potential, float *ave_Source, float *ave_Sigma, float *ave_Ex, float *ave_Ey);
__global__ void Accomulate_Particle_Density(int nsp, int Gsize, GPG *data);
__global__ void Average_Particle_Density(int nsp, int Gsize, int N_ave, Species *info, GPG *data);
__global__ void Average_Argon_MCC_rate(int Gsize, int TnRct, int N_ave, float dt, float *MCCR, float *ave_MCCR, Species *info);
__global__ void Average_Oxygen_MCC_rate(int Gsize, int TnRct, int N_ave, float dt, float *MCCR, float *ave_MCCR, Species *info);
__global__ void Average_ArO2_MCC_rate(int Gsize, int TnRct, int N_ave, float dt, float *MCCR, float *ave_MCCR, Species *info);
__global__ void Argon_Update_Source(int Gsize, int ngy, int TnRct, float*MCCR, CollF *CollP, GPG *SP, GFG *FG, GFG *FG_S, GGA *BG);
__global__ void Oxygen_Update_Source(int Gsize, int TnRct, float*MCCR, CollF *CollP, GPG *SP, GFG *FG, GFG *FG_S, GGA *BG);
__global__ void ArO2_Update_Source(int Gsize, int ngy, int TnRct, float*MCCR, CollF *CollP, GPG *SP, GFG *FG, GFG *FG_S, GGA *BG);
#endif