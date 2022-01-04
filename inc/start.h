#include "xypic.h"

extern double t;  // real time
extern float dt;   // timestsep for PIC
extern float dtc; // time step for continuity equationex
extern int tstep; // number of time step
extern int cstep; // Number of Cycle step
extern int DumpFlag; // Dump File ON,OFF
extern char InputFile[80]; // INPUT FILE NAME
extern char DumpFile[80];  // DUMP FILE NAME
extern char *ConstBFile;  // DUMP FILE NAME
extern void Argon_CrossSectionSET(CollF *CF);
extern void Oygen_CrossSectionSET(CollF *CF);
extern void ArO2_CrossSectionSET(CollF *CF);
extern void CG_Matrix_Setting(float *A, int *Ai, int *Aj, float **b, float *M, float *Atemp, float *btemp);
extern void SetParticleLoad(int isp, float Ninit, int load_type, float x_left,float x_right, float y_top, float y_bottom, float vti);
#ifndef __START_H__
#define __START_H__
// C Variable declaration
extern int PRINT_Flag;
extern int device_num;
extern float xlength,ylength,zlength;
extern int ngx,ngy,Gsize;
extern int ncx,ncy,Csize;
extern float dx,dy;
extern float idxy,idx,idy,dx2,dy2,dxdy2,hdx,hdy,r_eps0;
extern float dt_dx, dt_dy;
extern float fncx,fncy,fngx,fngy;
extern float *x_Garray,*y_Garray;
extern float *x_Carray,*y_Carray;
extern int BoundaryNUM;
extern int *BoundaryX0,*BoundaryY0,*BoundaryX1,*BoundaryY1,*BoundaryBC;
extern float *BoundaryTEMP;
extern int CondNUM,CondNUMR;  // CondNUMR is Real Conductor number
extern int*CondM_ID, *CondX0,*CondX1,*CondY0,*CondY1;
extern float *CondTEMP,*CondR,*CondL,*CondC;
extern int SrcNUM;
extern int *SrcM_ID;
extern float *SrcDC, *SrcPOWER, *SrcAC, *SrcFREQ, *Src2piFREQ, *SrcPHASE, *SrcRPHASE;
extern int External_Flag; // 0 : Voltage driven, 1: Power driven
extern float Min_FREQ, Max_FREQ;
extern int DielNUM,DielNUMR;
extern int *DielM_ID, *DielX0, *DielX1, *DielY0, *DielY1;
extern float *DielEPS;
extern GGA *vec_G;
extern GCA *vec_C;
extern int **StructureIndex;
extern int *vec_StructureIndex;
//
extern int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen
extern int nsp, nfsp, nBG;
extern Species *SP;// particle species
extern Fluid *FG;	// fluid species
extern BackG *BG;	// background species
extern HCP *PtD;
extern GPG *Host_G_sp;
extern float Total_Pressure;
//
extern int DT_PIC;  // Number of 1 cycle step
extern int DT_CONTI; // How many times PIC dt?
extern int Lap_Field_Solver_Test,Lap_Field_Solver_Flag,Lap_Field_Solver_Save;
extern float PCGtol;
extern float PCGtol2;
extern int HISTMAX;
extern int dHIST;
extern int NP_LIMIT; //Each of particle limit
extern int N_ave;
extern int N_smt;  // Number of smoothing every timestep
extern int ConstB_Flag; // Magnetic field 
extern int PD_intv;
extern float PD_Ratio;
extern int CSS_Flag;
//
extern float EndTime;
extern float Margin_ave_np;
extern int Flag_ave_np, Same_ave_np;
extern int Basic_Flag; // 0 : Basic, 1: OTHERS
//
extern int nRct_cx,nRct_rc; // Number of reaction _ cross section or Reaction rate
extern int TnRct; // Total Number of reaction 
extern int mMnum;
extern CollF *Coll_Flag;
extern int Msize;
extern float *MCC_rate,*ave_MCC_rate;
//
extern int FieldIter;
extern int A_size;
extern float *A_val,*TA_val;
extern int *Ai,*Aj;
extern int **A_idx;
extern float *MatM,**cond_b,*temp_b;
extern float *phi_cond;
extern float **AM,*V_t,*b_t,*extq,*extq_1,*extq_2,*extq_3;
extern float *CondCharge;
extern float *Surf_charge,*Old_Surf_charge,*Old2_Surf_charge;
//
extern int init_dump_num;
extern int OVER_dump_order;
extern int dump_order;
extern int dump_num;
extern float *dump_cycle;
extern int TecplotS_CX_Flag;
extern int TecplotS_Gsize_Flag;
extern int TecplotS_Particle_Flag;
extern int TecplotS_Particle_Num;
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
//
void InputRead();
void Geometry_setting();
float Face_To_Area(int Face);
void FieldSolverSetting();
void GasSetting();
int ***TIMalloc(int sizeX,int sizeY,int sizeZ);
float ***TFMalloc(int sizeX,int sizeY,int sizeZ);
int **MIMalloc(int sizeX,int sizeY);
float **MFMalloc(int sizeX,int sizeY);
int *VIMalloc(int size);
float *VFMalloc(int size);
void TFFree(float ***T,int sizeX,int sizeY);
void TIFree(int ***T,int sizeX,int sizeY);
void MFFree(float **M,int sizeX);
void MIFree(int **M,int sizeX);
void TFInit(float ***T,float C,int sizeX,int sizeY,int sizeZ);
void TIInit(int ***T,int C,int sizeX,int sizeY,int sizeZ);
void MFInit(float **M,float C,int sizeX,int sizeY);
void MIInit(int **M,int C,int sizeX,int sizeY);
void VFInit(float *V,float C,int size);
void VIInit(int *V,int C,int size);
void MFCopy(float **M,float **C,int sizeX,int sizeY);
void MICopy(int **M,int **C,int sizeX,int sizeY);
void VFCopy(float *V,float *C,int size);
void VICopy(int *V,int *C,int size);
int IVnZC(char A[50],int Value); // Int Value non Zero CHECK
float FIVnZC(char A[50],float Value);// float Value non Zero CHECK
void MFDigonal(float **M,float C,float D,int sizeX,int sizeY);
#endif