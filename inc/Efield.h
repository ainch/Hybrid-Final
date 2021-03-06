#include "xypic.h"
extern float dx,dy;
extern int CondNUM,CondNUMR;  // CondNUMR is Real Conductor number
extern int*CondM_ID, *CondX0,*CondX1,*CondY0,*CondY1;
extern float *CondTEMP,*CondR,*CondL,*CondC;
extern int SrcNUM;
extern int *SrcM_ID;
extern float *SrcDC, *SrcPOWER, *SrcAC, *SrcFREQ, *Src2piFREQ, *SrcPHASE, *SrcRPHASE;
extern int External_Flag; // 0 : no External , 1: External
extern float Min_FREQ, Max_FREQ;
extern int DT_PIC;  // Number of 1 cycle step
extern int DT_CONTI; // How many times PIC dt?
extern int DT_MCCn; // mcc count for each step
#ifndef __EFIELD_H__
#define __EFIELD_H__
extern float dt;   // timestsep for PIC
extern float dtc; // time step for continuity equation
extern float dt_dx, dt_dy;
extern int DT_MCCn; // mcc count for each step
extern float dt_mcc; // timestsep for MCC Module
extern int CYCLE_NUM; // Minimum frequency number of cycle
extern int *Efield_Flag, *Cond_Source_num, *Cond_count,**Cond_Power_ID;
extern float *Cond_Power;
extern float **CC_a; //Circuit_Const a
void Source_setting();
#endif