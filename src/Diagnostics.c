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
	// Potential
	// Source : Total Charge density 
   	// Potential : current potential
	// Sigma : [Dielectric] surface charge density, [Conductor] Surface current
   	vec_Potential = VFMalloc(Gsize); VFInit(vec_Potential,0.0,Gsize); // current potential
   	sum_Potential = VFMalloc(Gsize); VFInit(sum_Potential,0.0,Gsize); // accumulation potential
   	ave_Potential = VFMalloc(Gsize); VFInit(ave_Potential,0.0,Gsize); // averaged potential
   	vec_Source = VFMalloc(Gsize); VFInit(vec_Source,0.0,Gsize);
   	sum_Source = VFMalloc(Gsize); VFInit(sum_Source,0.0,Gsize);
   	ave_Source = VFMalloc(Gsize); VFInit(ave_Source,0.0,Gsize);
	vec_Sigma = VFMalloc(Gsize); VFInit(vec_Sigma,0.0,Gsize);
   	sum_Sigma = VFMalloc(Gsize); VFInit(sum_Sigma,0.0,Gsize);
   	ave_Sigma = VFMalloc(Gsize); VFInit(ave_Sigma,0.0,Gsize);
	sum_Ex = VFMalloc(Gsize); VFInit(sum_Ex,0.0,Gsize);
   	ave_Ex = VFMalloc(Gsize); VFInit(ave_Ex,0.0,Gsize);
	sum_Ey = VFMalloc(Gsize); VFInit(sum_Ey,0.0,Gsize);
   	ave_Ey = VFMalloc(Gsize); VFInit(ave_Ey,0.0,Gsize);
}
