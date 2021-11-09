#include "Viewer.h"

void Make_Value(){
	int isp,i,j,k,index;
	float sum;
	// Particle
	np = VIMalloc(nsp);
	sp_x = MFMalloc(nsp,(MAX_VIEW_NP-1));
	sp_y = MFMalloc(nsp,(MAX_VIEW_NP-1));
	sp_vx = MFMalloc(nsp,(MAX_VIEW_NP-1));
	sp_vy = MFMalloc(nsp,(MAX_VIEW_NP-1));
	sp_vz = MFMalloc(nsp,(MAX_VIEW_NP-1));
	VIInit(np,0,nsp);
	MFInit(sp_x, 0.0, nsp, (MAX_VIEW_NP-1));
	MFInit(sp_y, 0.0, nsp, (MAX_VIEW_NP-1));
	MFInit(sp_vx, 0.0, nsp, (MAX_VIEW_NP-1));
	MFInit(sp_vy, 0.0, nsp, (MAX_VIEW_NP-1));
	MFInit(sp_vz, 0.0, nsp, (MAX_VIEW_NP-1));
	int k1,k2;
	k1=0;k2=0;
	for(isp=0;isp<nsp;isp++){
        np[isp] = SP[isp].np;
		printf("Number of Particle %s : %d\n",SP[isp].name,SP[isp].np);
        if(np[isp]<(MAX_VIEW_NP-1)){
			for (i = 0; i <np[isp]; i++) {
				sp_x[isp][i] = dx*PtD[isp].x[i];
				sp_y[isp][i] = dy*PtD[isp].y[i];
				sp_vx[isp][i] = PtD[isp].vx[i];
				sp_vy[isp][i] = PtD[isp].vy[i];
				sp_vz[isp][i] = PtD[isp].vz[i];
        	}
		}else{
			np[isp] = (MAX_VIEW_NP-1);
			for (i = 0; i < np[isp]; i++) {
				index=frand()*((float)SP[isp].np-1.0f);
				sp_x[isp][i]=dx*PtD[isp].x[index];
				sp_y[isp][i]=dy*PtD[isp].y[index];
				sp_vx[isp][i]=PtD[isp].vx[index];
				sp_vy[isp][i]=PtD[isp].vy[index];
				sp_vz[isp][i]=PtD[isp].vz[index];
        	}
		}
    }
	sp_den = TFMalloc(nsp,ngx,ngy); 	TFInit(sp_den,0.0,nsp,ngx,ngy);  
	sp_ave_den = TFMalloc(nsp,ngx,ngy);	TFInit(sp_ave_den,0.0,nsp,ngx,ngy);
	sp_sigma = TFMalloc(nsp,ngx,ngy);	TFInit(sp_sigma,0.0,nsp,ngx,ngy);
	for (isp = 0; isp < nsp; isp++) {
		for (i = 0; i < ngx; i++) {
			for (j = 0; j < ngy; j++) {
				index = isp * Gsize + i * ngy + j;
				sp_den[isp][i][j] = Host_G_sp[index].den;
				sp_ave_den[isp][i][j] = Host_G_sp[index].ave_den;
				sp_sigma[isp][i][j] = Host_G_sp[index].sigma;
			}
		}
	}
	FG_den = TFMalloc(nfsp,ncx,ncy); 	TFInit(FG_den,0.0,nfsp,ncx,ncy);  
	FG_source = TFMalloc(nfsp,ncx,ncy);	TFInit(FG_source,0.0,nfsp,ncx,ncy);
	FG_flux_x = TFMalloc(nfsp,ngx,ngy); TFInit(FG_flux_x,0.0,nfsp,ngx,ngy);  
	FG_flux_y = TFMalloc(nfsp,ngx,ngy); TFInit(FG_flux_y,0.0,nfsp,ngx,ngy);  
	for (isp = 0; isp < nfsp; isp++) {
		for (i = 0; i < ncx; i++) {
			for (j = 0; j < ncy; j++) {
				index = isp * Csize + i * ncy + j;
				FG_den[isp][i][j] = 0.0f;
				//FG_source[isp][i][j] = 0.0f;
			}
		}
		for (i = 0; i < ngx; i++) {
			for (j = 0; j < ngy; j++) {
				//index = isp * Gsize + i * ngy + j;
				//FG_flux_x[isp][i][j] = 0.0f;
				//FG_flux_y[isp][i][j] = 0.0f;
			}
		}
	}
	GasTemp = MFMalloc(ngx,ngy); MFInit(GasTemp,0.0,ngx,ngy);
	GasDens = MFMalloc(ngx,ngy); MFInit(GasDens,0.0,ngx,ngy);
	GasVel = MFMalloc(ngx,ngy); MFInit(GasVel,0.0,ngx,ngy);
	GasDens2 = MFMalloc(ngx,ngy); MFInit(GasDens2,0.0,ngx,ngy);
	GasVel2 = MFMalloc(ngx,ngy); MFInit(GasVel2,0.0,ngx,ngy);
	LAP_pot = MFMalloc(ngx,ngy); MFInit(LAP_pot,0.0,ngx,ngy);
	POIS_pot = MFMalloc(ngx,ngy); MFInit(POIS_pot,0.0,ngx,ngy);
	Ex = MFMalloc(ngx,ngy); MFInit(Ex,0.0,ngx,ngy);
	Ey = MFMalloc(ngx,ngy); MFInit(Ey,0.0,ngx,ngy);
	Potential = MFMalloc(ngx,ngy); MFInit(Potential,0.0,ngx,ngy);
	AvePotential = MFMalloc(ngx,ngy); MFInit(AvePotential,0.0,ngx,ngy);
	AveSource = MFMalloc(ngx,ngy); MFInit(AveSource,0.0,ngx,ngy);
	AveSigma = MFMalloc(ngx,ngy); MFInit(AveSigma,0.0,ngx,ngy);
	AveEx = MFMalloc(ngx,ngy); MFInit(AveEx,0.0,ngx,ngy);
	AveEy = MFMalloc(ngx,ngy); MFInit(AveEy,0.0,ngx,ngy);
	for (i = 0; i < ngx; i++) {
		for (j = 0; j < ngy; j++) {
			index = i * ngy + j;
			GasTemp[i][j] = vec_G[index].Temp;
			GasDens[i][j] = vec_G[index].BackDen1;
			GasVel[i][j] = vec_G[index].BackVel1;
			Potential[i][j] = vec_Potential[index];
			AvePotential[i][j] = ave_Potential[index];
			AveSource[i][j] = ave_Source[index];
			AveSigma[i][j] = ave_Sigma[index];
			AveEx[i][j] = ave_Ex[index];
			AveEy[i][j] = ave_Ey[index];
			if(MainGas == 2){
				GasDens2[i][j] = vec_G[index].BackDen2;
				GasVel2[i][j] = vec_G[index].BackVel2;
			}
			LAP_pot[i][j] = vec_G[index].Lap_Pot;
			POIS_pot[i][j] = vec_G[index].Pois_Pot;
			Ex[i][j] = vec_G[index].Ex;
			Ey[i][j] = vec_G[index].Ey;
		}
	}
	MCC_Diag = TFMalloc(TnRct,ngx,ngy); TFInit(MCC_Diag,0.0,TnRct,ngx,ngy); 
	MCC_Value = VFMalloc(TnRct); VFInit(MCC_Value,0.0f,TnRct);
	for (isp = 0; isp < TnRct; isp++) {
		k=0;
		sum = 0.0f;
		for (i = 0; i < ngx; i++) {
			for (j = 0; j < ngy; j++) {
				index = i * ngy + j;
				MCC_Diag[isp][i][j] = ave_MCC_rate[index*TnRct + isp];
				sum += ave_MCC_rate[index*TnRct + isp];
				if(vec_G[i*ngy+j].DensRegion){
					k++;
				}
			}
		}
		MCC_Value[isp] = sum/k;
	}
}
void Set_Time_History(){
	char buffer[50];
    int isp,i,k;
	XGSet2D("linlin", "t (s)", "A Number of ptcl(t)", "closed", 100, 500, 1.0,
			1.0, TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
	for (isp = 0; isp < nsp; isp++)
		XGCurve(t_array, HistPt[isp].np, &hist_count, isp);

    //XGSet2D("linlog", "t (s)", "A Number of Fluid(t)", "closed", 100, 500, 1.0,
	//		1.0, TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
	//for (isp = 0; isp < nfsp; isp++)
	//	XGCurve(t_array, HistFG[isp].np, &hist_count, isp);

	XGSet2D("linlin", "t (s)", "History Averaged Particle (t)", "closed", 100, 500, 1.0,
			1.0, TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
	for (isp = 0; isp < nsp; isp++)
		XGCurve(t_ave_array, Hist_ave_Pt[isp].np, &hist_ave_count, isp);
	
	XGSet2D("linlin", "t (s)", "History Averaged Particle stack (t)", "closed", 100, 500, 1.0,
			1.0, TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
	for (isp = 0; isp < nsp; isp++)
		XGCurve(t_ave_array, Hist_ave_Pt_stack[isp].np, &hist_ave_count, isp);

	XGSet2D("linlin", "t (s)", "History PCG Iter (t)", "closed", 100, 500, 1.0, 1.0,
	TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
	XGScat2D(t_array, iter_array, &hist_count, 0);
	
	for (k = 0; k < CondNUMR; k++) {
		sprintf(buffer, "History Cond %d, Species Current", k, SP[isp].name);
		XGSet2D("linlin", "t (s)", buffer, "closed", 100, 500, 1.0, 1.0, TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
		for (isp = 0; isp < nsp; isp++) XGCurve(t_array, SP_current_hist[isp][k], &hist_count, isp);
		
		sprintf(buffer, "History Cond %d, Total Current", k);
		XGSet2D("linlin", "t (s)", buffer, "closed", 100, 500, 1.0, 1.0, TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
		XGCurve(t_array, Current_hist[k], &hist_count, 0);
		
		sprintf(buffer, "History Cond %d, Surf_charge_hist", k);
		XGSet2D("linlin", "t (s)", buffer, "closed", 100, 500, 1.0, 1.0, TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
		XGCurve(t_array, Surf_charge_hist[k], &hist_count, 0);
		
		sprintf(buffer, "History Cond %d, Constant Voltage", k);
		XGSet2D("linlin", "t (s)", buffer, "closed", 100, 500, 1.0, 1.0, TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
		XGCurve(t_array, Volt_hist[k], &hist_count, 0);
		
		sprintf(buffer, "History Cond %d, Circuit Voltage", k);
		XGSet2D("linlin", "t (s)", buffer, "closed", 100, 500, 1.0, 1.0, TRUE, TRUE, 0.0, 0.0, 0.0, 0.0);
		XGCurve(t_array, Volt_cond_hist[k], &hist_count, 0);
	}
}
void Set_PhaseSpace(){
	char buffer[50];
	int isp,i,k;
	for (isp = 0; isp < nsp; isp++) {
		sprintf(buffer, "A x-y Phase Space %s", SP[isp].name);
		XGSet2D("linlin", "x (m)", buffer, "closed", ncx, ncy, 1.0, 1.0, FALSE,
		FALSE, 0.0, xlength, 0.0, ylength);
		for (k = 0; k < CondNUM; k++)
			XGStructure(5, FILLED, 4, 4, CondX0[k] * dx, CondY0[k] * dy,
					CondX0[k] * dx, CondY1[k] * dy, CondX1[k] * dx,
					CondY1[k] * dy, CondX1[k] * dx, CondY0[k] * dy,
					CondX0[k] * dx, CondY0[k] * dy);
		for (k = 0; k < DielNUM; k++)
			XGStructure(5, FILLED, 5, 6, DielX0[k] * dx, DielY0[k] * dy,
					DielX0[k] * dx, DielY1[k] * dy, DielX1[k] * dx,
					DielY1[k] * dy, DielX1[k] * dx, DielY0[k] * dy,
					DielX0[k] * dx, DielY0[k] * dy);
		XGScat2D(sp_x[isp], sp_y[isp], &np[isp], isp);

		sprintf(buffer, "Vx-X Phase Space %s", SP[isp].name);
		XGSet2D("linlin", "x (m)", buffer, "closed", 200, 200, 1.0, 1.0, FALSE,
		TRUE, 0.0, xlength, 0.0, 0.0);
		XGScat2D(sp_x[isp], sp_vx[isp], &np[isp], isp);

		sprintf(buffer, "Vx-Y Phase Space %s", SP[isp].name);
		XGSet2D("linlin", "y (m)", buffer, "closed", 200, 200, 1.0, 1.0, FALSE,
		TRUE, 0.0, ylength, 0.0, 0.0);
		XGScat2D(sp_y[isp], sp_vx[isp], &np[isp], isp);

		sprintf(buffer, "Vy-X Phase Space %s", SP[isp].name);
		XGSet2D("linlin", "x (m)", buffer, "closed", 200, 200, 1.0, 1.0, FALSE,
		TRUE, 0.0, xlength, 0.0, 0.0);
		XGScat2D(sp_x[isp], sp_vy[isp], &np[isp], isp);

		sprintf(buffer, "Vy-Y Phase Space %s", SP[isp].name);
		XGSet2D("linlin", "y (m)", buffer, "closed", 200, 200, 1.0, 1.0, FALSE,
		TRUE, 0.0, ylength, 0.0, 0.0);
		XGScat2D(sp_y[isp], sp_vy[isp], &np[isp], isp);

		sprintf(buffer, "Vz-X Phase Space %s", SP[isp].name);
		XGSet2D("linlin", "x (m)", buffer, "closed", 200, 200, 1.0, 1.0, FALSE,
		TRUE, 0.0, xlength, 0.0, 0.0);
		XGScat2D(sp_x[isp], sp_vz[isp], &np[isp], isp);

		sprintf(buffer, "Vz-Y Phase Space %s", SP[isp].name);
		XGSet2D("linlin", "y (m)", buffer, "closed", 200, 200, 1.0, 1.0, FALSE,
		TRUE, 0.0, ylength, 0.0, 0.0);
		XGScat2D(sp_y[isp], sp_vz[isp], &np[isp], isp);

		sprintf(buffer, "Vx-Vy Phase Space %s", SP[isp].name);
		XGSet2D("linlin", "Vx (m/s)", buffer, "closed", 200, 200, 1.0, 1.0,
		TRUE,TRUE, 0.0, 0.0, 0.0, 0.0);
		XGScat2D(sp_vx[isp], sp_vy[isp], &np[isp], isp);
	}
}
void Set_Particle_2D(){
	char buffer[50];
	int isp;
	for (isp = 0; isp < nsp; isp++) {
		sprintf(buffer, "Density %s", SP[isp].name);
		XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 0.0, 60.0, "closed", 100,
				150, 1.0, 1.0, SP[isp].np2c / dx/dy, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
		XGSurf(x_Garray, y_Garray, sp_den[isp], &ngx, &ngy, 3);
		sprintf(buffer, "A Averaged Density %s", SP[isp].name);
		XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
		XGSurf(x_Garray, y_Garray, sp_ave_den[isp], &ngx, &ngy, 3);
		sprintf(buffer, "Sigma %s", SP[isp].name);
		XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
		XGSurf(x_Garray, y_Garray, sp_sigma[isp], &ngx, &ngy, 3);
	}
}
void Set_Fluid_2D(){
	char buffer[50];
	int isp;
	for (isp = 0; isp < nfsp; isp++) {
		sprintf(buffer, "Density %s", FG[isp].name);
		XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
		XGSurf(x_Carray, y_Carray, FG_den[isp], &ncx, &ncy, 3);
		/*
		sprintf(buffer, "Source %s", FG[isp].name);
		XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
		XGSurf(x_Carray, y_Carray, FG_source[isp], &ncx, &ncy, 3);
		sprintf(buffer, "Flux_x %s", FG[isp].name);
		XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
		XGSurf(x_Garray, y_Garray, FG_flux_x[isp], &ngx, &ngy, 3);
		sprintf(buffer, "Flux_y %s", FG[isp].name);
		XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
		XGSurf(x_Garray, y_Garray, FG_flux_y[isp], &ngx, &ngy, 3);
		*/
	}
}
void Set_Background_2D(){
	char buffer[50];
	int isp;
	sprintf(buffer, "Gas Temperature");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, GasTemp, &ngx, &ngy, 3);
	sprintf(buffer, "Gas Density %s", BG[0].name);
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, GasDens, &ngx, &ngy, 3);
	sprintf(buffer, "Gas Velocity %s", BG[0].name);
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, GasVel, &ngx, &ngy, 3);
	if(MainGas==2){
		sprintf(buffer, "Gas Density %s", BG[1].name);
		XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
		XGSurf(x_Garray, y_Garray, GasDens2, &ngx, &ngy, 3);
		sprintf(buffer, "Gas Velocity %s", BG[1].name);
		XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
		XGSurf(x_Garray, y_Garray, GasVel2, &ngx, &ngy, 3);
	}
	sprintf(buffer, "Ex");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, Ex, &ngx, &ngy, 3);
	sprintf(buffer, "Ey");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, Ey, &ngx, &ngy, 3);
	sprintf(buffer, "Potential Laplace (V)");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, LAP_pot, &ngx, &ngy, 3);
	sprintf(buffer, "Potential Poisson (V)");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, POIS_pot, &ngx, &ngy, 3);
	sprintf(buffer, "Potential Total (V)");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, Potential, &ngx, &ngy, 3);
}
void Set_FIELD_2D(){
	char buffer[50];
	int isp;
	sprintf(buffer, "Averaged Potential");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, AvePotential, &ngx, &ngy, 3);
	sprintf(buffer, "Averaged Source");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, AveSource, &ngx, &ngy, 3);
	sprintf(buffer, "Averaged Sigma");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, AveSigma, &ngx, &ngy, 3);
	sprintf(buffer, "Averaged Ex");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, AveEx, &ngx, &ngy, 3);
	sprintf(buffer, "Averaged Ey");
	XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 45.0, -60.0, "closed", 100,
			150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
			xlength, 0.0, ylength, 0.0, 0.0);
	XGSurf(x_Garray, y_Garray, AveEy, &ngx, &ngy, 3);
	
}
void Set_MCC_RATE(){
	char buffer[50];
	int i;
	if(MainGas == ARGON){
		for(i=0;i<5;i++){
			sprintf(buffer, "Z MCC E Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=5;i<7;i++){
			sprintf(buffer, "Z MCC Arion Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=7;i<TnRct;i++){
			sprintf(buffer, "Z MCC Constant Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
	}else if(MainGas == OXYGEN){
		for(i=0;i<41;i++){
			sprintf(buffer, "Z MCC E Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=41;i<47;i++){
			sprintf(buffer, "Z MCC Onion Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=47;i<53;i++){
			sprintf(buffer, "Z MCC O2ion Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=53;i<58;i++){
			sprintf(buffer, "Z MCC Oion Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=58;i<TnRct;i++){
			sprintf(buffer, "Z MCC Constant Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
	}else if(MainGas == ARO2){
		for(i=0;i<46;i++){
			sprintf(buffer, "Z MCC E Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=46;i<52;i++){
			sprintf(buffer, "Z MCC Onion Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=52;i<60;i++){
			sprintf(buffer, "Z MCC O2ion Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=60;i<65;i++){
			sprintf(buffer, "Z MCC Oion Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=65;i<68;i++){
			sprintf(buffer, "Z MCC Arion Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	
		for(i=68;i<TnRct;i++){
			sprintf(buffer, "Z MCC Constant Reaction %d ave=%1.1g",i,MCC_Value[i]);
			XGSet3D("linlinlin", "X (m)", "Y (m)", buffer, 20.0, 60.0, "closed", 100,
				150, 1.0, 1.0, 1.0, FALSE, FALSE, TRUE, 0.0,
				xlength, 0.0, ylength, 0.0, 0.0);
			XGSurf(x_Garray, y_Garray, MCC_Diag[i], &ngx, &ngy, 3);
		}	

	}
	
}
void Make_Plot(){
	// Time History
    Set_Time_History();
	// Particle
	Set_PhaseSpace();
	// DENSITY
	Set_Particle_2D();
	Set_Fluid_2D();
	Set_Background_2D();
	Set_FIELD_2D();
	Set_MCC_RATE();
}
void Dump(){
}
void Quit(){
}
