#include "Diagnostics.h"

void Diagnostic_Setting(){
	int isp;
	nave_count = 0;	
	hist_count=0;
   	Hcount = 1;
   	t_array = (float *) malloc(HISTMAX * sizeof(float));
   	iter_array = (float *) malloc(HISTMAX * sizeof(float));
   	HistPt = (Hist *) malloc(nsp * sizeof(Hist));
   	HistFG = (Hist *) malloc(nfsp * sizeof(Hist));
   	for(isp=0;isp<nsp;isp++){
      	HistPt[isp].np = (float *) malloc(HISTMAX * sizeof(float));
      	VFInit(HistPt[isp].np,0.0,HISTMAX);
   	}
   	for(isp=0;isp<nfsp;isp++){
      	HistFG[isp].np = (float *) malloc(HISTMAX * sizeof(float));
      	VFInit(HistFG[isp].np,0.0,HISTMAX);
   	}
}
