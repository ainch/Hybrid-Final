#include "cuda_Diagnostic.cuh"

void Diagnostic_Basic(){
    int i, j, k, isp, index;
    static int power_init = 0;
    float cond_current, dis_current, phi_now, power_total;
    int buf1=0, now_np;
    float buf;
    if(nave_count==N_ave){  // average calculate
        Conti_Flag = 1;  
        Average_Particle_Density<<<DIAG_NSPG_GRID, DIAG_NSPG_BLOCK>>>(nsp, Gsize, N_ave, dev_info_sp, dev_G_sp);
        Average_Field_Data<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, N_ave, TotPotential,dev_Source,dev_Sigma,dev_GvecSet
                        ,dev_sum_Potential,dev_sum_Source,dev_sum_Sigma,dev_sum_Ex,dev_sum_Ey
                        ,dev_ave_Potential,dev_ave_Source,dev_ave_Sigma,dev_ave_Ex,dev_ave_Ey);
        switch(MainGas){ 
		case ARGON:
            Average_Argon_MCC_rate<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, N_ave, dt, dev_MCC_rate, dev_ave_MCC_rate, dev_info_sp);
            break;
		case OXYGEN: 
            Average_Oxygen_MCC_rate<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, N_ave, dt, dev_MCC_rate, dev_ave_MCC_rate, dev_info_sp);
			break;
		case ARO2: 
            Average_ArO2_MCC_rate<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, N_ave, dt, dev_MCC_rate, dev_ave_MCC_rate, dev_info_sp);
			break;
        }
        nave_count = 0;
    }else{ // accomulation data
        Accomulate_Particle_Density<<<DIAG_NSPG_GRID, DIAG_NSPG_BLOCK>>>(nsp, Gsize, dev_G_sp);
        Accomulate_Field_Data<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize,TotPotential,dev_Source,dev_Sigma,dev_GvecSet
                        ,dev_sum_Potential,dev_sum_Source,dev_sum_Sigma,dev_sum_Ex,dev_sum_Ey);
        nave_count++;
    }  
    // Steady-state check
    if(Flag_ave_np){
        for(isp = 0;isp<nsp;isp++){
            cudaMemcpy(&now_np, &dev_info_sp[isp].np, sizeof(int),cudaMemcpyDeviceToHost);
            ave_np[isp] += (int)((float)now_np/(float)DT_PIC);
        }
        if(tstep%DT_PIC == 0){
            buf1 = 0;
            for(isp = 0;isp<nsp;isp++){
                new_ave_np[isp] = (float)ave_np[isp];
                ave_np[isp] = 0;
                buf = 100.0f * fabs((new_ave_np[isp]-old_ave_np[isp])/new_ave_np[isp]);
                old_ave_np[isp] = new_ave_np[isp];
                if(buf < Margin_ave_np){
                    Stack_ave_np[isp]++;
                }else{
                    Stack_ave_np[isp] = 0;
                }
                if(Stack_ave_np[isp] >= Same_ave_np){
                    buf1++;
                }
                t_ave_array[hist_ave_count] = (float) t;
                Hist_ave_Pt[isp].np[hist_ave_count] = new_ave_np[isp];
                Hist_ave_Pt_stack[isp].np[hist_ave_count] = (float)Stack_ave_np[isp];
            }
            hist_ave_count++;
            if(buf1>=nsp){
                Flag_ave_np = 0;
                Basic_Flag = 0;
            }
        }
    }
    if(Basic_Flag<-1) if(cstep > abs(Basic_Flag)) Basic_Flag = 0;
    // Calculate Current for Power driven or External circuit 
    cudaMemcpy(Host_G_buf, dev_Sigma, Gsize * sizeof(float),	cudaMemcpyDeviceToHost);
	for (i = 0; i < CondNUMR; i++) {
		Old2_Surf_charge[i] = Old_Surf_charge[i];
		Old_Surf_charge[i] = Surf_charge[i];
		Surf_charge[i] = 0.0f;
	}
	for (i = 0; i < Gsize; i++) {
		if (vec_G[i].CondID) {
			index = vec_G[i].CondID - 1;
			Surf_charge[index] += Host_G_buf[i] * vec_G[i].Area;
		}
	}
    // Power driven and Dual frequency;
    cudaMemcpy(CondCharge, dev_CondCharge, nsp * CondNUMR * sizeof(float),cudaMemcpyDeviceToHost);
    power_total = 0;
    for (i = 0; i < CondNUMR; i++) {
        // Current calculator
        Current_Now[i] = 0.0f;
        cond_current = 0.0f;
        for (isp = 0; isp < nsp; isp++) {
			cond_current += SP[isp].q_density * CondCharge[isp*CondNUMR + i] / dt;
		}
        dis_current = 0.5 * (3 * Surf_charge[i] - 4 * Old_Surf_charge[i] + Old2_Surf_charge[i]) / dt;
		Current_Now[i] = dis_current - cond_current;
        // Power calculator
        if(Cond_Source_num[i] == 1){ // Single Power Calculator
            //Accumulation
            phi_now = V_t[i];
			Cond_Power[i] += phi_now * Current_Now[i];
            Cond_count[i]++;
            //Averaged power
            if(Cond_count[i] == PD_intv){
                Cond_Power[i] /= (float)Cond_count[i];
                power_total+=Cond_Power[i];
                Cond_count[i] = 0;
                index = Cond_Power_ID[i][0];
                if(Cond_Power[i]>0) phi_now = SrcAC[index]*pow(fabs(SrcPOWER[index]/Cond_Power[i]),0.2);
                if(power_init) {
					if(phi_now>(1.0f + PD_Ratio)*SrcAC[index]){
						SrcAC[index]*=1.05;
					}
					else{
						SrcAC[index]=phi_now;
					}
				}
				printf("\nConductor %d : Fix Power %g, Current Power %g, Voltage %g\n",i,SrcPOWER[index],Cond_Power[i],SrcAC[index]);
				Cond_Power[i] = 0.0f;
                Cond_count[i]=0;
				power_init++;
            }
        }else if(Cond_Source_num[i] == 2){ // Dual Power Calculator

        }
    }
    if (hist_count >= HISTMAX) {
        for (k = 1, i = 4; k < HISTMAX / 4; k++, i += 4) {
			t_array[k] = t_array[i];
			for (isp = 0; isp < nsp; isp++) {
                HistPt[isp].np[k] =  HistPt[isp].np[i];
			}
			iter_array[k] = iter_array[i];
            for (j = 0; j < CondNUMR; j++) {
				for (isp = 0; isp < nsp; isp++) {
					SP_current_hist[isp][j][k] = SP_current_hist[isp][j][i];
				}
				Current_hist[j][k] = Current_hist[j][i];
				Surf_charge_hist[j][k] = Surf_charge_hist[j][i];
				Volt_hist[j][k] = Volt_hist[j][i];
				Volt_cond_hist[j][k] = Volt_cond_hist[j][i];
			}
		}
		hist_count = k;
		dHIST *= 4;
    }
    if((--Hcount)==0){
        t_array[hist_count] = (float) t;
        iter_array[hist_count] = (float) *FIter;
        for (isp = 0; isp < nsp; isp++) {
            cudaMemcpy(&now_np, &dev_info_sp[isp].np, sizeof(int),cudaMemcpyDeviceToHost);
            HistPt[isp].np[hist_count] = now_np;
        }
        for (i = 0; i < CondNUMR; i++) {
            Current_hist[i][hist_count] = Current_Now[i];
            Surf_charge_hist[i][hist_count] = Surf_charge[i];
            Volt_hist[i][hist_count] = V_t[i];
		    Volt_cond_hist[i][hist_count] = phi_cond[i];
            for (isp = 0; isp < nsp; isp++) {
			    SP_current_hist[isp][i][hist_count] = SP[isp].q_density * CondCharge[isp*CondNUMR + i] / dt;
		    }
        }
        hist_count++;
        Hcount = dHIST;
    }
}
__global__ void Accomulate_Field_Data(int Gsize, float *TotPot, float *Source, float *Sigma, GGA *vecG
                        ,float *sum_Potential, float *sum_Source, float *sum_Sigma, float *sum_Ex, float *sum_Ey){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize) return;
    sum_Potential[TID] += TotPot[TID];
    //if(TotPot[TID] !=0) printf("TotPot[%d]= %g\n",TID,TotPot[TID]);
    sum_Source[TID] += Source[TID];
    sum_Sigma[TID] += Sigma[TID];
    sum_Ex[TID] += vecG[TID].Ex;
    sum_Ey[TID] += vecG[TID].Ey;           
}
__global__ void Average_Field_Data(int Gsize, int N_ave, float *TotPot, float *Source, float *Sigma, GGA *vecG
                        ,float *sum_Potential, float *sum_Source, float *sum_Sigma, float *sum_Ex, float *sum_Ey
                        ,float *ave_Potential, float *ave_Source, float *ave_Sigma, float *ave_Ex, float *ave_Ey){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize) return;
    float oneofN = 1/(float)N_ave;
    ave_Potential[TID] = oneofN * sum_Potential[TID];// average
    //printf("ave_Potential[%d]= %g\n",TID,oneofN);
    sum_Potential[TID] = TotPot[TID];// Init && start
    ave_Source[TID] = oneofN * sum_Source[TID];
    sum_Source[TID] = Source[TID];
    ave_Sigma[TID] = oneofN * sum_Sigma[TID];
    sum_Sigma[TID] = Sigma[TID];
    ave_Ex[TID] = oneofN * sum_Ex[TID];
    sum_Ex[TID] = vecG[TID].Ex;
    ave_Ey[TID] = oneofN * sum_Ey[TID];    
    sum_Ey[TID] = vecG[TID].Ey;   
}
__global__ void Accomulate_Particle_Density(int nsp, int Gsize, GPG *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=nsp*Gsize) return;
    data[TID].sum_den += data[TID].den;
}
__global__ void Average_Particle_Density(int nsp, int Gsize, int N_ave, Species *info, GPG *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=nsp*Gsize) return;
    int isp = TID/Gsize;
    float temp = info[isp].Denscale/(float)N_ave;
    data[TID].sum_den += data[TID].den;
    data[TID].ave_den = temp * data[TID].sum_den;
    data[TID].sum_den = 0.0f;
}
__global__ void Average_Argon_MCC_rate(int Gsize, int TnRct, int N_ave, float dt, float *MCCR, float *ave_MCCR, Species *info){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize) return;
    int i;
    float temp = info[0].Denscale/(float)N_ave/dt;

    for(i=0;i<7;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    ave_MCCR[TID*TnRct+8] = temp * MCCR[TID*TnRct+8];
    MCCR[TID*TnRct+8] = 0.0f;
}
__global__ void Average_Oxygen_MCC_rate(int Gsize, int TnRct, int N_ave, float dt, float *MCCR, float *ave_MCCR, Species *info){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize) return;
    int i;
    float temp;

    temp = info[0].Denscale/(float)N_ave/dt;
    for(i=0;i<41;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    temp = info[3].Denscale/(float)N_ave/dt;
    for(i=41;i<47;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    temp = info[1].Denscale/(float)N_ave/dt;
    for(i=47;i<53;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    temp = info[2].Denscale/(float)N_ave/dt;
    for(i=53;i<58;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    temp = 1/(float)N_ave/dt;
    ave_MCCR[TID*TnRct+58] = info[3].Denscale * temp * MCCR[TID*TnRct+58];
    MCCR[TID*TnRct+58] = 0.0f;
}
__global__ void Average_ArO2_MCC_rate(int Gsize, int TnRct, int N_ave, float dt, float *MCCR, float *ave_MCCR, Species *info){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize) return;
    int i;
    float temp;

    temp = info[0].Denscale/(float)N_ave/dt;
    for(i=0;i<46;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    temp = info[4].Denscale/(float)N_ave/dt;
    for(i=46;i<52;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    temp = info[2].Denscale/(float)N_ave/dt;
    for(i=52;i<60;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    temp = info[3].Denscale/(float)N_ave/dt;
    for(i=60;i<65;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    temp = info[1].Denscale/(float)N_ave/dt;
    for(i=65;i<68;i++){
        ave_MCCR[TID*TnRct+i] = temp * MCCR[TID*TnRct+i];
        MCCR[TID*TnRct+i] = 0.0f;
    }
    temp = 1/(float)N_ave/dt;
    ave_MCCR[TID*TnRct+69] = info[3].Denscale * temp * MCCR[TID*TnRct+69];
    MCCR[TID*TnRct+69] = 0.0f;
    ave_MCCR[TID*TnRct+70] = info[3].Denscale * temp * MCCR[TID*TnRct+70];
    MCCR[TID*TnRct+70] = 0.0f;
    ave_MCCR[TID*TnRct+79] = info[1].Denscale * temp * MCCR[TID*TnRct+79];
    MCCR[TID*TnRct+79] = 0.0f;
    ave_MCCR[TID*TnRct+80] = info[1].Denscale * temp * MCCR[TID*TnRct+80];
    MCCR[TID*TnRct+80] = 0.0f;
    ave_MCCR[TID*TnRct+81] = info[0].Denscale * temp * MCCR[TID*TnRct+81];
    MCCR[TID*TnRct+81] = 0.0f;
}
void Set_Diagnostic_cuda(){
    // Host BUF VECTOR
    Host_G_buf = VFMalloc(Gsize);
    Host_C_buf = VFMalloc(Csize);
    VFInit(Host_G_buf,0.0,Gsize);
    VFInit(Host_C_buf,0.0,Csize);

    checkCudaErrors(cudaMalloc((void**)&dev_sum_Potential, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_sum_Potential, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&dev_ave_Potential, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_ave_Potential, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemcpy(dev_ave_Potential, ave_Potential, Gsize * sizeof(float), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&dev_sum_Source, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_sum_Source, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&dev_ave_Source, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_ave_Source, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemcpy(dev_ave_Source, ave_Source, Gsize * sizeof(float), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&dev_sum_Sigma, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_sum_Sigma, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&dev_ave_Sigma, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_ave_Sigma, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemcpy(dev_ave_Sigma, ave_Sigma, Gsize * sizeof(float), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&dev_sum_Ex, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_sum_Ex, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&dev_ave_Ex, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_ave_Ex, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemcpy(dev_ave_Ex, ave_Ex, Gsize * sizeof(float), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&dev_sum_Ey, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_sum_Ey, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**)&dev_ave_Ey, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_ave_Ey, 0.0f, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemcpy(dev_ave_Ey, ave_Ey, Gsize * sizeof(float), cudaMemcpyHostToDevice));
}
	