#include "xypic.h"
extern double t;  // real time
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern char DumpFile[80];  // DUMP FILE NAME
extern int Gsize,ngx,ngy;
extern int Csize,ncx,ncy;
extern float dx,dy,dV;
extern float xlength,ylength,zlength;
extern float *x_Garray,*y_Garray;
extern float *x_Carray,*y_Carray;
extern int CondNUM,CondNUMR;
extern int *CondM_ID,*CondX0,*CondX1,*CondY0,*CondY1;
extern int DielNUM;
extern int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
extern GGA *vec_G;
extern GCA *vec_C;
extern GCP *Host_sp;
extern GPG *Host_G_sp;
extern GFC *Host_C_F;
extern GFG *Host_G_F;
extern float *MCC_rate, *ave_MCC_rate;
extern float *vec_Potential, *ave_Potential; // [Gsize] potential
extern float *vec_Source, *ave_Source; // [Gsize] Charge density 
extern float *vec_Sigma, *ave_Sigma; // [Gsize] [Dielectric] surface charge density, [Conductor] Surface current
extern float *ave_Ex, *ave_Ey;
extern int **A_idx;
extern int MainGas;
extern int nsp,nfsp;
extern int TnRct;
extern Species *SP;// particle species
extern Fluid *FG;	// fluid species
extern BackG *BG;	// background species
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
extern float frand();
extern int ***TIMalloc(int sizeX,int sizeY,int sizeZ);
extern float ***TFMalloc(int sizeX,int sizeY,int sizeZ);
extern int **MIMalloc(int sizeX,int sizeY);
extern float **MFMalloc(int sizeX,int sizeY);
extern int *VIMalloc(int size);
extern float *VFMalloc(int size);
extern void TFFree(float ***T,int sizeX,int sizeY);
extern void TIFree(int ***T,int sizeX,int sizeY);
extern void MFFree(float **M,int sizeX);
extern void MIFree(int **M,int sizeX);
extern void TFInit(float ***T,float C,int sizeX,int sizeY,int sizeZ);
extern void TIInit(int ***T,int C,int sizeX,int sizeY,int sizeZ);
extern void MFInit(float **M,float C,int sizeX,int sizeY);
extern void MIInit(int **M,int C,int sizeX,int sizeY);
extern void VFInit(float *V,float C,int size);
extern void VIInit(int *V,int C,int size);
extern void MFCopy(float **M,float **C,int sizeX,int sizeY);
extern void MICopy(int **M,int **C,int sizeX,int sizeY);
extern void VFCopy(float *V,float *C,int size);
extern void VICopy(int *V,int *C,int size);
extern int IVnZC(char A[50],int Value); // Int Value non Zero CHECK
extern float FIVnZC(char A[50],float Value);// float Value non Zero CHECK
extern void MFDigonal(float **M,float C,float D,int sizeX,int sizeY);
#ifndef __VIEWER_H__
#define __VIEWER_H__
int *np;
float **sp_x,**sp_y,**sp_vx,**sp_vy,**sp_vz;
float ***sp_den, ***sp_ave_den, ***sp_sigma;
float ***FG_den, ***FG_ave_den, ***FG_source, ***FG_flux_x, ***FG_flux_y;
float **GasTemp, **GasDens, **GasVel,**GasDens2, **GasVel2, **LAP_pot, **POIS_pot, **Ex, **Ey;
float **Potential, **AvePotential,**AveSource,**AveSigma,**AveEx,**AveEy;
float ***MCC_Diag;
float *MCC_Value;
void Make_Value();
void Make_Plot();
void Set_PhaseSpace();
void Set_Time_History();
void Set_Particle_2D();
void Set_Fluid_2D();
void Set_Background_2D();
void Set_FIELD_2D();
void Set_MCC_RATE();
void Dump();
void Quit();
#endif