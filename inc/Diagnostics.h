#include "xypic.h"
extern int nsp,nfsp;
extern int HISTMAX;
extern int dHIST;
#ifndef __DIAGNOSTICS_H__
#define __DIAGNOSTICS_H__
int Hcount;
int hist_count;
int nave_count;
float *t_array;
float *iter_array;    
Hist *HistPt,*HistFG;
void Diagnostic_Setting();
#endif

