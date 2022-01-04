#include "def.h"
#include "cuda.h"
#include "cublas_v2.h"
#include "cusparse_v2.h"
#include "cusolverDn.h"
#include "cusolverSp.h"
#include "cusolverRf.h"
#include "cusolver_common.h"
#include "curand_kernel.h"
#include "cuda_runtime.h"
#include "cufft.h"
#include "helper_cuda.h"
#include "interop.cuh"

int Conti_Flag;
GFC *Fluid_sp;
GFG *Fluid_Den, *Fluid_Src;
int Conti_xnum, Conti_ynum;
Con_RegionX *Conti_x;
Con_RegionY *Conti_y;
int *dev_StructureIndex;
int *ReArgFlag;

int h_nvel;
float *vsave;
float time_sum;
float gputime;
cudaEvent_t start, stop;
float 	totaltime,gputime_field,gputime_efield;
float 	gputime_move,gputime_mcc,gputime_deposit;
float 	gputime_diag,gputime_sort,gputime_Tec;
float 	gputime_continue,gputime_dump;
int		TotalT_D;
int		TotalT_H;
int		TotalT_M;
int		TotalT_S;

int PRINT_Flag;
int device_num;
float xlength,ylength,zlength;
int ngx,ngy,Gsize;
int ncx,ncy,Csize;
float dx,dy;
float idxy,idx,idy,dx2,dy2,dxdy2,hdx,hdy,r_eps0;
float fncx,fncy,fngx,fngy;
float *x_Garray,*y_Garray;
float *x_Carray,*y_Carray;
int BoundaryNUM;
int *BoundaryX0,*BoundaryY0,*BoundaryX1,*BoundaryY1,*BoundaryBC;
float *BoundaryTEMP;
int CondNUM,CondNUMR;  // CondNUMR is Real Conductor number
int*CondM_ID, *CondX0,*CondX1,*CondY0,*CondY1;
float *CondTEMP,*CondR,*CondL,*CondC;
int SrcNUM;
int *SrcM_ID;
float *SrcDC, *SrcPOWER, *SrcAC, *SrcFREQ, *Src2piFREQ, *SrcPHASE, *SrcRPHASE;
int External_Flag; // 0 : Voltage driven, 1: Power driven
float Min_FREQ, Max_FREQ;
int DielNUM,DielNUMR;
int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
float *DielEPS;
GGA *vec_G;
GCA *vec_C;
int **StructureIndex;
int *vec_StructureIndex;
int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen
int nsp, nfsp, nBG;
Species *SP;// particle species
Fluid *FG;	// fluid species
BackG *BG;	// background species
HCP *PtD;
GPG *Host_G_sp;
float Total_Pressure;
int DT_PIC;  // Number of 1 cycle step
int DT_CONTI; // How many times PIC dt?
int Lap_Field_Solver_Test,Lap_Field_Solver_Flag,Lap_Field_Solver_Save;
float PCGtol;
float PCGtol2;
int HISTMAX;
int dHIST;
int NP_LIMIT; //Each of particle limit
int N_ave;
int N_smt;  // Number of smoothing every timestep
int ConstB_Flag; // Magnetic field 
int PD_intv;
float PD_Ratio;
int CSS_Flag;
float EndTime;
float Margin_ave_np;
int Flag_ave_np, Same_ave_np;
int Basic_Flag; // 0 : Basic, 1: OTHERS
int nRct_cx,nRct_rc; // Number of reaction _ cross section or Reaction rate
int TnRct; // Total Number of reaction 
int mMnum;
CollF *Coll_Flag;
int Msize;
float *MCC_rate,*ave_MCC_rate;
int FieldIter;
int A_size;
float *A_val,*TA_val;
int *Ai,*Aj;
int **A_idx;
float *MatM,**cond_b,*temp_b;
float *phi_cond;
float **AM,*V_t,*b_t,*extq,*extq_1,*extq_2,*extq_3;
float *CondCharge;
float *Surf_charge,*Old_Surf_charge,*Old2_Surf_charge;
int init_dump_num;
int OVER_dump_order;
int dump_order;
int dump_num;
float *dump_cycle;
int TecplotS_CX_Flag;
int TecplotS_Gsize_Flag;
int TecplotS_Particle_Flag;
int TecplotS_Particle_Num;
int TecplotS_2D_Flag;
int TecplotS_2D_Ncycle;
int TecplotS_Movie_Flag;
int TecplotS_Movie_Ncycle;
int TecplotS_Movie_Frame;
int TecplotS_Movie_SCYCLE;
int TecplotS_Movie_Count;
int TecplotS_PT_Movie_Flag;
int TecplotS_PT_Movie_Ncycle;
int TecplotS_PT_Movie_Frame;
int TecplotS_PT_Movie_SCYCLE;
int TecplotS_PT_Movie_Count;

int nz;
int N;
float *dev_A,*dev_TA;	
int *dev_Aj,*dev_Ai;
float *dev_b,*dev_Tb;
int   *vec_A_idx;
int   *dev_A_idx;
float *dev_M;
float *dev_AP,*dev_X,*dev_R,*dev_Z,*dev_P;
double *dot_result;
float *dot_result2;
int *FIter;
GGA *dev_GvecSet;
GCA *dev_CvecSet;
float *dev_CondCharge;
float *dev_phi;     // PCG Solution
float *dev_phi_buf; // Sigma or buf
size_t pitch;
float *Lap_TEMP_Sol; // Temperature Solution of Laplace Equation
float *Lap_PHI_Sol; // Each of conductor Phi Solution of Laplace Equation, This is Device value
float *TotPotential;
float *LapPotential;
float *dev_Sigma, *dev_Source;
											   
											   
											   


float LOGX_MIN,LOGX_MAX;
float dLOGX,idLOGX;
int N_LOGX;
ArCollD *Ar_Data;
O2CollD *O2_Data;
ArO2CollD *ArO2_Data;
float *Host_G_buf, *Host_C_buf;
float *dev_sum_Potential, *dev_ave_Potential;
float *dev_sum_Source, *dev_ave_Source;
float *dev_sum_Sigma, *dev_ave_Sigma;
float *dev_sum_Ex, *dev_ave_Ex;
float *dev_sum_Ey, *dev_ave_Ey;





curandState *devStates;
float *dev_vsave;
int sMemSize;
int sMemSize_MCC;
dim3 FIELD_GRID,FIELD_BLOCK;
dim3 FIELD_GRID2,FIELD_BLOCK2;
dim3 DEPOSIT_GRID,DEPOSIT_BLOCK;
dim3 EFIELD_GRID,EFIELD_BLOCK;
dim3 MOVE_GRID, MOVE_BLOCK;
dim3 SORT_GRID, SORT_BLOCK;
dim3 MCC_GRID, MCC_BLOCK;
dim3 DIAG_G_GRID, DIAG_G_BLOCK;
dim3 DIAG_NSPG_GRID, DIAG_NSPG_BLOCK;

float *dev_MCC_rate,*dev_ave_MCC_rate;
float *dev_RCstack,*dev_stack;
CollF *dev_Coll_Flag;
ArCollD *dev_ArCX;
O2CollD *dev_O2CX;
ArO2CollD *dev_ArO2CX;
MCC_sigmav *Host_SigmaV;
MCC_sigmav *dev_SigmaV;

long int seed; // related to Random number
double t;  // real time
int tstep; // number of time step
int cstep; // Number of Cycle step
int DumpFlag; // Dump File ON,OFF
char InputFile[80]; // INPUT FILE NAME
char DumpFile[80];  // DUMP FILE NAME
char *ConstBFile;  // DUMP FILE NAME
Fluid *dev_FG;
GFG *dev_FG_Den, *dev_FG_Src;

int Total_maxnp;
GCP *Host_sp, *dev_sp;
GPG *dev_G_sp;
Species *dev_info_sp;
cudaDeviceProp prop;
float dt;   // timestsep for PIC
float dtc; // time step for continuity equation
float dt_dx, dt_dy;
int DT_MCCn; // mcc count for each step
float dt_mcc; // timestsep for MCC Module
int CYCLE_NUM; // Minimum frequency number of cycle
int *Efield_Flag, *Cond_Source_num, *Cond_count,**Cond_Power_ID;
float *Cond_Power;
float **CC_a; //Circuit_Const a
int Hcount;
int hist_count, hist_ave_count;
int nave_count;
float *t_array;
float *iter_array;    
Hist *HistPt,*HistFG;
float *t_ave_array;
Hist *Hist_ave_Pt,*Hist_ave_Pt_stack;
float **Current_hist,**Surf_charge_hist,**Volt_hist,**Volt_cond_hist;
float ***SP_current_hist;
float *Current_Now;
float *vec_Potential, *sum_Potential, *ave_Potential; // [Gsize] potential
float *vec_Source, *sum_Source, *ave_Source; // [Gsize] Charge density 
float *vec_Sigma, *sum_Sigma, *ave_Sigma; // [Gsize] [Dielectric] surface charge density, [Conductor] Surface current
float *sum_Ex, *ave_Ex, *sum_Ey, *ave_Ey;
int *Stack_ave_np, *ave_np;
float *new_ave_np, *old_ave_np;

int PT_Movie_S_count;