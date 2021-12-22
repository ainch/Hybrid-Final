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
int PRINT_Flag;
int device_num;
float xlength,ylength,zlength;
int ngx,ngy,Gsize;
int ncx,ncy,Csize;
float dx,dy;
float idxy,idx,idy,dx2,dy2,dxdy2,hdx,hdy,r_eps0;
float dt_dx, dt_dy;
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
//
int MainGas; // Gas type 0:argon, 1:oxygen, 2:argon/oxygen
int nsp, nfsp, nBG;
Species *SP;// particle species
Fluid *FG;	// fluid species
BackG *BG;	// background species
HCP *PtD;
GPG *Host_G_sp;
float Total_Pressure;
//
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
//
float EndTime;
float Margin_ave_np;
int Flag_ave_np, Same_ave_np;
int Basic_Flag; // 0 : Basic, 1: OTHERS
//
int nRct_cx,nRct_rc; // Number of reaction _ cross section or Reaction rate
int TnRct; // Total Number of reaction 
int mMnum;
CollF *Coll_Flag;
int Msize;
float *MCC_rate,*ave_MCC_rate;
//
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
//
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