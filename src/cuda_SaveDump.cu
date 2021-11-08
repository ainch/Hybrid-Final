#include "cuda_SaveDump.cuh"
void V001_DUMP(FILE *SF){

}
void V000_DUMP(FILE *SF){
    int isp,i;
    // time
    fwrite(&t, 8, 1, SF);
    fwrite(&tstep, 4, 1, SF);
    fwrite(&cstep, 4, 1, SF);
    // Particle Data gpu >> cpu
    cudaMemcpy(SP, dev_info_sp,  nsp * sizeof(Species), cudaMemcpyDeviceToHost);
    cudaMemcpy(Host_sp, dev_sp, Total_maxnp * sizeof(GCP), cudaMemcpyDeviceToHost);
    cudaMemcpy(Host_G_sp, dev_G_sp, nsp * Gsize * sizeof(GPG), cudaMemcpyDeviceToHost);
    Copy_GCPtoHCP(SP, Host_sp, PtD, Host_G_sp);
    // Gas info
    fwrite(&MainGas, 4, 1, SF);
    for(isp=0;isp<nsp;isp++){
        // NP2C
        fwrite(&SP[isp].np2c, 4, 1, SF);
        fwrite(&SP[isp].np, 4, 1, SF);
        for (i = 0; i < SP[isp].np; i++) {
            fwrite(&PtD[isp].CellID[i], 4, 1, SF);
            fwrite(&PtD[isp].x[i], 4, 1, SF);
            fwrite(&PtD[isp].y[i], 4, 1, SF);
            fwrite(&PtD[isp].vx[i], 4, 1, SF);
            fwrite(&PtD[isp].vy[i], 4, 1, SF);
            fwrite(&PtD[isp].vz[i], 4, 1, SF);
        }
    }
    // History data
    fwrite(&hist_count, 4, 1, SF);
    fwrite(&dHIST, 4, 1, SF);
	fwrite(t_array, 4, hist_count, SF);
	for (isp = 0; isp < nsp; isp++) {
		fwrite(HistPt[isp].np, 4, hist_count, SF);
	}
	fwrite(iter_array, 4, hist_count, SF);
    for (i = 0; i < CondNUMR; i++) {
		fwrite(Current_hist[i], 4, hist_count, SF);
		for (isp = 0; isp < nsp; isp++) {
			fwrite(SP_current_hist[isp][i], 4, hist_count, SF);
		}
		fwrite(Volt_hist[i], 4, hist_count, SF);
		fwrite(Volt_cond_hist[i], 4, hist_count, SF);
		fwrite(Surf_charge_hist[i], 4, hist_count, SF);
	}
    // 2D data
    for (i = 0; i < nsp*Gsize; i++) {
        fwrite(&Host_G_sp[i].den, 4, 1, SF);    // pt density : GPG
        fwrite(&Host_G_sp[i].ave_den, 4, 1, SF); // pt Ave_density : GPG
        fwrite(&Host_G_sp[i].sigma, 4, 1, SF); // pt sigma
    }
    cudaMemcpy(vec_G, dev_GvecSet, Gsize * sizeof(GGA), cudaMemcpyDeviceToHost);
    for (i = 0; i < Gsize; i++) {
        fwrite(&vec_G[i].Temp, 4, 1, SF); // BG density,Tem,vel : GGA
        fwrite(&vec_G[i].BackVel1, 4, 1, SF); // BG density,Tem,vel : GGA
        fwrite(&vec_G[i].BackDen1, 4, 1, SF); // BG density,Tem,vel : GGA
        fwrite(&vec_G[i].BackVel2, 4, 1, SF); // BG density,Tem,vel : GGA
        fwrite(&vec_G[i].BackDen2, 4, 1, SF); // BG density,Tem,vel : GGA
        fwrite(&vec_G[i].Lap_Pot, 4, 1, SF); // BG Lap_Pot
        fwrite(&vec_G[i].Pois_Pot, 4, 1, SF); // BG Pois_Pot
        fwrite(&vec_G[i].Ex, 4, 1, SF); // BG Ex
        fwrite(&vec_G[i].Ey, 4, 1, SF); // BG Ey
    } 
    // Field
    cudaMemcpy(vec_Potential, TotPotential, Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    fwrite(vec_Potential, 4, Gsize, SF); 
    cudaMemcpy(ave_Potential, dev_ave_Potential, Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    fwrite(ave_Potential, 4, Gsize, SF); 
    cudaMemcpy(ave_Source, dev_ave_Source, Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    fwrite(ave_Source, 4, Gsize, SF); 
    cudaMemcpy(ave_Sigma, dev_ave_Sigma, Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    fwrite(ave_Sigma, 4, Gsize, SF); 
    cudaMemcpy(ave_Ex, dev_ave_Ex, Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    fwrite(ave_Ex, 4, Gsize, SF); 
    cudaMemcpy(ave_Ey, dev_ave_Ey, Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    fwrite(ave_Ey, 4, Gsize, SF); 
    // MCC
    cudaMemcpy(ave_MCC_rate, dev_ave_MCC_rate, TnRct * Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    fwrite(ave_MCC_rate, 4, TnRct * Gsize, SF); // MCC AVE RATE : mccrate
}
void SaveDumpFile(int KEY2,int KEY1,int KEY0){
    FILE *SaveFile;
	char FileName[512];
    int isp;
    float time_sum;
    /////// Function access management ///////
    if (init_dump_num){
		while (1){ // find dump_order
			if (tstep <= dump_cycle[dump_order])
				break;
			else
				dump_order++;
		}
		init_dump_num--;
	}
	if (tstep < dump_cycle[dump_order]){
		return;
    }
    if(dump_order >= dump_num){
        OVER_dump_order++;
        if(OVER_dump_order == 100 * DT_PIC){
            OVER_dump_order = 0;
        }else{
            return;
        }
    }
	dump_order++;
    //////////////////////////////////////////
    // Dump load version SAVE.  Ver.[KEY2].[KEY1].[KEY0]
    //  KEY0 : 0~9 If you just add and remove the storage variable
    //  KEY1 : 0~9 If you just add and remove the storage variable
    //  KEY2 : When there is a change in a significant calculation module
    // Version History
    //  Ver.0.0.0 : Time, Particle information, np2c, Number of particle,
    //
    /// Open Dump File
	fprintf(stderr,"\n-------------------------Dumping File Ver.[%d][%d][%d]---------------------------\n",KEY2,KEY1,KEY0);
    sprintf(FileName, "%s.dmp%d", InputFile, dump_order);
	if ((SaveFile = fopen(FileName, "w")) == NULL) {
		puts("Dump: open failed");
		exit(-1);
	}
    //Save Start
    fwrite(&KEY2, 4, 1, SaveFile);
    fwrite(&KEY1, 4, 1, SaveFile);
    fwrite(&KEY0, 4, 1, SaveFile);
    if(KEY2>=0 && KEY1>=0 && KEY0>=0)  V000_DUMP(SaveFile);
    if(KEY2>=0 && KEY1>=0 && KEY0>=1)  V001_DUMP(SaveFile);
    //Save End
    fclose(SaveFile);

    // time calculate
	while(totaltime > 1000){
		TotalT_S++;
		totaltime = totaltime - 1000;
	}
	while(TotalT_S >= 60){
			TotalT_M++;
			TotalT_S -= 60;
	}
	while(TotalT_M >= 60){
			TotalT_H++;
			TotalT_M -= 60;
	}
	while(TotalT_H >= 24){
			TotalT_D++;
			TotalT_H -= 24;
	}
    fprintf(stderr, "Dump at t=%1.5e(s),Step[%d]Cycle[%d] %s\n", t,tstep,cstep,FileName);
	for (isp = 0; isp < nsp; isp++){
		fprintf(stderr, "%s : %d,  ", SP[isp].name, SP[isp].np);
	}fprintf(stderr, "\n");
	fprintf(stderr, "Domain size : %d X %d =%d,  ", ngx, ngy, ngx * ngy);
	fprintf(stderr, "Time: %d(d), %d(h), %d(m), %d(s)\n",TotalT_D,TotalT_H,TotalT_M,TotalT_S);
	time_sum = gputime_field+gputime_efield+gputime_move+gputime_sort+gputime_mcc+gputime_continue+gputime_deposit+gputime_diag+gputime_Tec+gputime_dump;
	fprintf(stderr, "Total : time = %2.8f	(s)\n", time_sum * 0.001);
	fprintf(stderr, "Field	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_field * 0.001, gputime_field * 100 / time_sum);
	fprintf(stderr, "Efield	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_efield * 0.001, gputime_efield * 100 / time_sum);
	fprintf(stderr, "Move	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_move * 0.001, gputime_move * 100 / time_sum);
	fprintf(stderr, "Sort	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_sort * 0.001, gputime_sort * 100 / time_sum);
	fprintf(stderr, "Mcc	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_mcc * 0.001, gputime_mcc * 100 / time_sum);
	fprintf(stderr, "CONTI	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_continue * 0.001, gputime_continue * 100 / time_sum);
	fprintf(stderr, "Depo	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_deposit * 0.001, gputime_deposit * 100 / time_sum);
	fprintf(stderr, "Diag	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_diag * 0.001, gputime_diag * 100 / time_sum);
	fprintf(stderr, "Tecplot: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_Tec * 0.001, gputime_Tec * 100 / time_sum);
	fprintf(stderr, "Dump	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_dump * 0.001, gputime_dump * 100 / time_sum);
	fprintf(stderr, "------------------------------------------------------------------------------\n");
}