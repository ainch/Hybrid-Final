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
      	VFInit(HistPt[isp].np,0.0f,HISTMAX);
   	}
   	for(isp=0;isp<nfsp;isp++){
      	HistFG[isp].np = (float *) malloc(HISTMAX * sizeof(float));
      	VFInit(HistFG[isp].np,0.0f,HISTMAX);
   	}
	hist_ave_count=0;
	t_ave_array = VFMalloc(HISTMAX); VFInit(t_ave_array,0.0f,HISTMAX);
	Hist_ave_Pt = (Hist *) malloc(nsp * sizeof(Hist));
	Hist_ave_Pt_stack = (Hist *) malloc(nsp * sizeof(Hist));
	for(isp=0;isp<nsp;isp++){
      	Hist_ave_Pt[isp].np = (float *) malloc(HISTMAX * sizeof(float));
      	VFInit(Hist_ave_Pt[isp].np,0.0f,HISTMAX);
		Hist_ave_Pt_stack[isp].np = (float *) malloc(HISTMAX * sizeof(float));
      	VFInit(Hist_ave_Pt_stack[isp].np,0.0f,HISTMAX);
   	}
	// Source
	Current_hist = MFMalloc(CondNUMR,HISTMAX); MFInit(Current_hist,0.0f,CondNUMR,HISTMAX);
	Surf_charge_hist = MFMalloc(CondNUMR,HISTMAX); MFInit(Surf_charge_hist,0.0f,CondNUMR,HISTMAX);
	Volt_hist = MFMalloc(CondNUMR,HISTMAX); MFInit(Volt_hist,0.0F,CondNUMR,HISTMAX);  // Voltage
	Volt_cond_hist = MFMalloc(CondNUMR,HISTMAX); MFInit(Volt_cond_hist,0.0f,CondNUMR,HISTMAX); // Circuit Total Voltage
	SP_current_hist = TFMalloc(nsp,CondNUMR,HISTMAX); TFInit(SP_current_hist,0.0f,nsp,CondNUMR,HISTMAX);
	Current_Now = VFMalloc(CondNUMR);	 VFInit(Current_Now, 0.0f, CondNUMR);	// Current now each of conductor
	
	// Potential
	// Source : Total Charge density 
   	// Potential : current potential
	// Sigma : [Dielectric] surface charge density, [Conductor] Surface current
   	vec_Potential = VFMalloc(Gsize); VFInit(vec_Potential,0.0f,Gsize); // current potential
   	ave_Potential = VFMalloc(Gsize); VFInit(ave_Potential,0.0f,Gsize); // averaged potential
   	vec_Source = VFMalloc(Gsize); VFInit(vec_Source,0.0f,Gsize);
   	ave_Source = VFMalloc(Gsize); VFInit(ave_Source,0.0f,Gsize);
	vec_Sigma = VFMalloc(Gsize); VFInit(vec_Sigma,0.0f,Gsize);
   	ave_Sigma = VFMalloc(Gsize); VFInit(ave_Sigma,0.0f,Gsize);
   	ave_Ex = VFMalloc(Gsize); VFInit(ave_Ex,0.0f,Gsize);
   	ave_Ey = VFMalloc(Gsize); VFInit(ave_Ey,0.0f,Gsize);
	
	// average np check
	ave_np = VIMalloc(nsp); VIInit(ave_np,0,nsp);
	Stack_ave_np = VIMalloc(nsp); VIInit(Stack_ave_np,0,nsp);
    new_ave_np = VFMalloc(nsp); VFInit(new_ave_np,0.0f,nsp);
    old_ave_np = VFMalloc(nsp); VFInit(old_ave_np,0.0f,nsp);
	for(isp = 0;isp<nsp;isp++){
        old_ave_np[isp] = SP[isp].np;
		Hist_ave_Pt[isp].np[hist_ave_count] = SP[isp].np;
    }
	hist_ave_count++;
}
