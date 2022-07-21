#include "cuda_Diagnostic.cuh"

#include <cooperative_groups.h>
#include <cuda_runtime.h>

namespace cg = cooperative_groups;

static float * Host_G_buf_tmp;
static float * Surf_charge_tmp;

static void __global__ update
(
	float * Surf_charge_tmp,
	GGA * const Field, 
	float * const dev_phi_buf,
	int Gsize,
	int CondNUMR
)
{
	cg::thread_block group_block  = cg::this_thread_block();
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	__shared__ float a_shd[128];
	
	int id_cond = 0;
	float area;
	float potential;
	if(tid < Gsize)
	{
		id_cond = Field[tid].CondID;
		area = Field[tid].Area;
		potential = dev_phi_buf[tid];
        if(id_cond != 0)
           atomicAdd(Surf_charge_tmp + id_cond - 1, area * potential);
	}
}

void Diagnostic(){
    int i, j, k, isp, index;
    static int power_init = 0;
    float cond_current, dis_current, phi_now, power_total;
    int buf1=0, now_np[nsp];
    float oldDen[nfsp],newDen[nfsp];
    float buf;
    // Cal NP
    for (isp = 0; isp < nsp; isp++) {
		cudaMemset((void *) &dev_info_sp[isp].np, 0, sizeof(int));
		SumReductionINT1024All<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, isp, dev_G_sp, &dev_info_sp[isp].np);
		cudaMemcpy(&SP[isp].np, &dev_info_sp[isp].np, sizeof(int),cudaMemcpyDeviceToHost);
        now_np[isp] = SP[isp].np;
	}
    if(nave_count==N_ave){  // average calculate
        Conti_Flag = 1;  
        Average_Particle_Density<<<DIAG_NSPG_GRID, DIAG_NSPG_BLOCK>>>(nsp, Gsize, N_ave, dev_info_sp, dev_G_sp);
        Average_Field_Data<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, N_ave, TotPotential,dev_Source,dev_Sigma,dev_GvecSet
                        ,dev_sum_Potential,dev_sum_Source,dev_sum_Sigma,dev_sum_Ex,dev_sum_Ey
                        ,dev_ave_Potential,dev_ave_Source,dev_ave_Sigma,dev_ave_Ex,dev_ave_Ey);
        switch(MainGas){ 
		case ARGON:
            Average_Argon_MCC_rate<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, N_ave, dt, dev_MCC_rate, dev_ave_MCC_rate, dev_info_sp);
            Argon_Update_Source<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, ngy, TnRct, dev_ave_MCC_rate, dev_Coll_Flag, dev_G_sp, dev_FG_Den, dev_FG_Src, dev_GvecSet);
            break;
		case OXYGEN: 
            Average_Oxygen_MCC_rate<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, N_ave, dt, dev_MCC_rate, dev_ave_MCC_rate, dev_info_sp);
			Oxygen_Update_Source<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, dev_ave_MCC_rate, dev_Coll_Flag, dev_G_sp, dev_FG_Den, dev_FG_Src, dev_GvecSet);
            break;
		case ARO2: 
            Average_ArO2_MCC_rate<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, N_ave, dt, dev_MCC_rate, dev_ave_MCC_rate, dev_info_sp);
			ArO2_Update_Source<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, ngy, TnRct, dev_ave_MCC_rate, dev_Coll_Flag, dev_G_sp, dev_FG_Den, dev_FG_Src, dev_GvecSet);
            break;
        }
        cudaMemcpy(Fluid_Src, dev_FG_Src, nfsp * Gsize * sizeof(GFG), cudaMemcpyDeviceToHost);
        if(CSS_Flag){
            for(isp = 0;isp<nfsp;isp++){ // Continuity on
                FG[isp].CSS_Flag = 1;
            }
        }
        Sync_Fluid_GFGtoGFC_forSource(Fluid_Src, Fluid_sp);
        nave_count = 0;
    }else{ // accomulation data
        Accomulate_Particle_Density<<<DIAG_NSPG_GRID, DIAG_NSPG_BLOCK>>>(nsp, Gsize, dev_G_sp);
        Accomulate_Field_Data<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize,TotPotential,dev_Source,dev_Sigma,dev_GvecSet
                        ,dev_sum_Potential,dev_sum_Source,dev_sum_Sigma,dev_sum_Ex,dev_sum_Ey);
        nave_count++;
    }  
    if(Conti_Flag){//Fluid density copy
        if(CSS_Flag){
            for(isp = 0;isp<nfsp;isp++){
                oldDen[isp] = FG[isp].ave_Den;
            }
        }
        Sync_Fluid_GFCtoGFG_forDen(Fluid_sp, Fluid_Den); 
        cudaMemcpy(dev_FG_Den, Fluid_Den, nfsp * Gsize * sizeof(GFG), cudaMemcpyHostToDevice);
        if(CSS_Flag){
            for(isp = 0;isp<nfsp;isp++){
                newDen[isp] = FG[isp].ave_Den;
                if(nave_count > FG[isp].CSS_Check){
                    if(abs(oldDen[isp]-newDen[isp])/newDen[isp] <= FG[isp].CSS_Conver)
                        FG[isp].CSS_Flag = 1;
                }
            }
        }   
    }
    // Steady-state check
    if(Flag_ave_np){
        for(isp = 0;isp<nsp;isp++){
            ave_np[isp] += (int)((float)now_np[isp]/(float)DT_PIC);
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
    
    if(true)
    {
        checkCudaErrors(cudaMemset((void *) Surf_charge_tmp, 0.0, CondNUMR * sizeof(float)));
            
        dim3 dim_num_block = dim3(Gsize / 128 + 1);
        dim3 dim_size_block = dim3(128);
    
        cudaMemcpy
        (
            Host_G_buf_tmp, Host_G_buf, Gsize * sizeof(float), cudaMemcpyHostToDevice
        );
    
        update<<<dim_num_block, dim_size_block>>>
        (
            Surf_charge_tmp, dev_GvecSet,
            Host_G_buf_tmp, Gsize, CondNUMR
        );
    
        checkCudaErrors
        (
            cudaMemcpy
            (
                Surf_charge, Surf_charge_tmp, CondNUMR * sizeof(float), cudaMemcpyDeviceToHost
            )
        );
    }
    else
    {
        for (i = 0; i < Gsize; i++) {
            if (vec_G[i].CondID) {
                index = vec_G[i].CondID - 1;
                Surf_charge[index] += Host_G_buf[i] * vec_G[i].Area;
            }
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
            for (isp = 0; isp < nfsp; isp++) {
                HistFG[isp].np[k] = HistFG[isp].np[i];
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
            HistPt[isp].np[hist_count] = now_np[isp];
        }
        for (isp = 0; isp < nfsp; isp++) {
            HistFG[isp].np[hist_count] = FG[isp].ave_Den;
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
void Diagnostic_Basic(){
    int i, j, k, isp, index;
    static int power_init = 0;
    float cond_current, dis_current, phi_now, power_total;
    int buf1=0, now_np[nsp];
    float oldDen[nfsp],newDen[nfsp];
    float buf;
    // Cal NP
    for (isp = 0; isp < nsp; isp++) {
		cudaMemset((void *) &dev_info_sp[isp].np, 0, sizeof(int));
		SumReductionINT1024All<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, isp, dev_G_sp, &dev_info_sp[isp].np);
		cudaMemcpy(&SP[isp].np, &dev_info_sp[isp].np, sizeof(int),cudaMemcpyDeviceToHost);
        now_np[isp] = SP[isp].np;
	}
    if(nave_count==N_ave){  // average calculate
        Conti_Flag = 1;  
        Average_Particle_Density<<<DIAG_NSPG_GRID, DIAG_NSPG_BLOCK>>>(nsp, Gsize, N_ave, dev_info_sp, dev_G_sp);
        Average_Field_Data<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, N_ave, TotPotential,dev_Source,dev_Sigma,dev_GvecSet
                        ,dev_sum_Potential,dev_sum_Source,dev_sum_Sigma,dev_sum_Ex,dev_sum_Ey
                        ,dev_ave_Potential,dev_ave_Source,dev_ave_Sigma,dev_ave_Ex,dev_ave_Ey);
        switch(MainGas){ 
		case ARGON:
            Average_Argon_MCC_rate<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, N_ave, dt, dev_MCC_rate, dev_ave_MCC_rate, dev_info_sp);
            Argon_Update_Source<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, ngy, TnRct, dev_ave_MCC_rate, dev_Coll_Flag, dev_G_sp, dev_FG_Den, dev_FG_Src, dev_GvecSet);
            break;
		case OXYGEN: 
            Average_Oxygen_MCC_rate<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, N_ave, dt, dev_MCC_rate, dev_ave_MCC_rate, dev_info_sp);
			Oxygen_Update_Source<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, dev_ave_MCC_rate, dev_Coll_Flag, dev_G_sp, dev_FG_Den, dev_FG_Src, dev_GvecSet);
            break;
		case ARO2: 
            Average_ArO2_MCC_rate<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, TnRct, N_ave, dt, dev_MCC_rate, dev_ave_MCC_rate, dev_info_sp);
			ArO2_Update_Source<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize, ngy, TnRct, dev_ave_MCC_rate, dev_Coll_Flag, dev_G_sp, dev_FG_Den, dev_FG_Src, dev_GvecSet);
            break;
        }
        cudaMemcpy(Fluid_Src, dev_FG_Src, nfsp * Gsize * sizeof(GFG), cudaMemcpyDeviceToHost);
        Sync_Fluid_GFGtoGFC_forSource(Fluid_Src, Fluid_sp);
        if(CSS_Flag){
            for(isp = 0;isp<nfsp;isp++){ // Continuity on
                FG[isp].CSS_Flag = 1;
            }
        }
        nave_count = 0;
    }else{ // accomulation data
        Accomulate_Particle_Density<<<DIAG_NSPG_GRID, DIAG_NSPG_BLOCK>>>(nsp, Gsize, dev_G_sp);
        Accomulate_Field_Data<<<DIAG_G_GRID, DIAG_G_BLOCK>>>(Gsize,TotPotential,dev_Source,dev_Sigma,dev_GvecSet
                        ,dev_sum_Potential,dev_sum_Source,dev_sum_Sigma,dev_sum_Ex,dev_sum_Ey);
        nave_count++;
    }  
    if(Conti_Flag){//Fluid density copy
        if(CSS_Flag){
            for(isp = 0;isp<nfsp;isp++){
                oldDen[isp] = FG[isp].ave_Den;
            }
        }
        Sync_Fluid_GFCtoGFG_forDen(Fluid_sp, Fluid_Den); 
        cudaMemcpy(dev_FG_Den, Fluid_Den, nfsp * Gsize * sizeof(GFG), cudaMemcpyHostToDevice);
        if(CSS_Flag){
            for(isp = 0;isp<nfsp;isp++){
                newDen[isp] = FG[isp].ave_Den;
                if(nave_count > FG[isp].CSS_Check){
                    if(abs(oldDen[isp]-newDen[isp])/newDen[isp] <= FG[isp].CSS_Conver)
                        FG[isp].CSS_Flag = 1;
                }
            }
        }   
    }
    // Steady-state check
    if(Flag_ave_np){
        for(isp = 0;isp<nsp;isp++){
            ave_np[isp] += (int)((float)now_np[isp]/(float)DT_PIC);
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
            for (isp = 0; isp < nfsp; isp++) {
                HistFG[isp].np[k] = HistFG[isp].np[i]; 
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
            HistPt[isp].np[hist_count] = now_np[isp];
        }
        for (isp = 0; isp < nfsp; isp++) {
            HistFG[isp].np[hist_count] = FG[isp].ave_Den;
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
__global__ void SumReductionINT1024All(int n, int isp, GPG *g_data, int *g_max){
	__shared__ int sdata[1024];
	unsigned int TID=threadIdx.x;
	unsigned int i=blockDim.x*blockIdx.x+TID ; // global thread index
    
	//if(i>=n) g_data[isp*n + i].PtNumInCell = 0;
    if(i>=n) return;
    
	sdata[TID] = g_data[isp*n + i].PtNumInCell;
	__syncthreads();

	if(TID<512) sdata[TID]=sdata[TID]+sdata[TID+512];
		__syncthreads();

	if(TID<256) sdata[TID]=sdata[TID]+sdata[TID+256];
		__syncthreads();

	if(TID<128) sdata[TID]=sdata[TID]+sdata[TID+128];
		__syncthreads();

	if(TID<64) sdata[TID]=sdata[TID]+sdata[TID+64];
		__syncthreads();

	if(TID<32) warpSumReduceINT(sdata,TID);

	if(TID==0) {
		atomicAdd(g_max,sdata[0]);
	}
}
__device__ void warpSumReduceINT(volatile int* sdata,int TID)
{
	sdata[TID]+=sdata[TID+32];
	sdata[TID]+=sdata[TID+16];
	sdata[TID]+=sdata[TID+8];
	sdata[TID]+=sdata[TID+4];
	sdata[TID]+=sdata[TID+2];
	sdata[TID]+=sdata[TID+1];
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
__global__ void Argon_Update_Source(int Gsize, int ngy, int TnRct, float*MCCR, CollF *CollP, GPG *SP, GFG *FG, GFG *FG_S, GGA *BG){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize) return;
    //float Eden = 0.25*(SP[TID].ave_den+SP[TID+ngy].ave_den+SP[TID+1].ave_den+SP[TID+ngy+1].ave_den);
    float Eden = 0.25*SP[TID].ave_den;
    MCCR[TID*TnRct + 7] = CollP[7].RR*Eden*FG[TID].n;
    MCCR[TID*TnRct + 9] = CollP[9].RR*BG[TID].BackDen1*FG[TID].n;      
    FG_S[TID].n = MCCR[TID*TnRct + 2] - MCCR[TID*TnRct + 4] - MCCR[TID*TnRct + 7] - 2 * MCCR[TID*TnRct + 8] - MCCR[TID*TnRct + 9];
} 
__global__ void Oxygen_Update_Source(int Gsize, int TnRct, float*MCCR, CollF *CollP, GPG *SP, GFG *FG, GFG *FG_S, GGA *BG){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize) return;
    int RID;
    float Sum_O2A, Diff_O2A, Sum_O2B, Diff_O2B,Sum_OP, Diff_OP,Sum_OD, Diff_OD;
    RID = TnRct * TID;

    MCCR[RID + 59] = CollP[59].RR*FG[2*Gsize + TID].n*FG[3*Gsize + TID].n;
    MCCR[RID + 60] = CollP[60].RR*FG[3*Gsize + TID].n*BG[TID].BackDen1;
    MCCR[RID + 61] = CollP[61].RR*FG[3*Gsize + TID].n*BG[TID].BackDen1;
    MCCR[RID + 62] = CollP[62].RR*FG[3*Gsize + TID].n*BG[TID].BackDen1;
    MCCR[RID + 63] = CollP[63].RR*FG[TID].n*FG[2*Gsize + TID].n;
    MCCR[RID + 64] = CollP[64].RR*FG[TID].n*BG[TID].BackDen1;   
    MCCR[RID + 65] = CollP[65].RR*FG[TID].n*FG[TID].n;
    MCCR[RID + 66] = CollP[66].RR*FG[Gsize + TID].n*BG[TID].BackDen1;  
    //
    Sum_O2A  = MCCR[RID + 3] + MCCR[RID + 61];
    Diff_O2A = MCCR[RID + 14] + MCCR[RID + 15] + MCCR[RID + 16] + MCCR[RID + 17] + MCCR[RID + 18]
                + MCCR[RID + 19] + MCCR[RID + 20] + MCCR[RID + 21] + MCCR[RID + 46] + MCCR[RID + 51]
                 + MCCR[RID + 56] + MCCR[RID + 63] + MCCR[RID + 64] + 2*MCCR[RID + 65];
    Sum_O2B = MCCR[RID + 4] + MCCR[RID + 62];
    Diff_O2B = MCCR[RID + 22] + MCCR[RID + 23] + MCCR[RID + 24] + MCCR[RID + 25] + MCCR[RID + 26]
                 + MCCR[RID + 27] + MCCR[RID + 28] + MCCR[RID + 29] + MCCR[RID + 52] + MCCR[RID + 57]
                 + MCCR[RID + 58] + MCCR[RID + 66];
    Sum_OP = MCCR[RID + 6] + 2 * MCCR[RID + 7] + MCCR[RID + 8] + MCCR[RID + 11] + MCCR[RID + 13]
                + MCCR[RID + 15] + 2*MCCR[RID + 18] + MCCR[RID + 19] + MCCR[RID + 21] + MCCR[RID + 23]
                 + 2*MCCR[RID + 26] + MCCR[RID + 27] + MCCR[RID + 29] + MCCR[RID + 30] + MCCR[RID + 31]
                  + MCCR[RID + 40] + MCCR[RID + 42] + MCCR[RID + 44] + 2*MCCR[RID + 45] + MCCR[RID + 46]
                   + MCCR[RID + 50] + MCCR[RID + 53] + MCCR[RID + 56] + MCCR[RID + 57] + MCCR[RID + 58]
                    + MCCR[RID + 59]  + MCCR[RID + 60]  + MCCR[RID + 61] + MCCR[RID + 62];
    Diff_OP =  MCCR[RID + 33] + MCCR[RID + 34] + MCCR[RID + 35] + MCCR[RID + 36] + MCCR[RID + 37]
                 + MCCR[RID + 38] + MCCR[RID + 43] + MCCR[RID + 47];
    Sum_OD = MCCR[RID + 8] + 2*MCCR[RID + 9] + MCCR[RID + 19] + 2*MCCR[RID + 20] + MCCR[RID + 27]
                 + 2*MCCR[RID + 28] + MCCR[RID + 31] + MCCR[RID + 33];
    Diff_OD = MCCR[RID + 39] + MCCR[RID + 40] + MCCR[RID + 59] + MCCR[RID + 60] + MCCR[RID + 61] + MCCR[RID + 62];

    FG_S[TID].n = Sum_O2A - Diff_O2A;
    FG_S[Gsize + TID].n = Sum_O2B - Diff_O2B;
    FG_S[2*Gsize + TID].n = Sum_OP - Diff_OP;
    FG_S[3*Gsize + TID].n = Sum_OD - Diff_OD;   
}
__global__ void ArO2_Update_Source(int Gsize, int ngy, int TnRct, float*MCCR, CollF *CollP, GPG *SP, GFG *FG, GFG *FG_S, GGA *BG){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize) return;
    int RID;
    float Eden;
    float Sum_ARM, Diff_ARM,Sum_O2A, Diff_O2A, Sum_O2B, Diff_O2B,Sum_OP, Diff_OP,Sum_OD, Diff_OD;
    float N1,N2,N4,N5;
    //Eden = 0.25*(SP[TID].ave_den+SP[TID+ngy].ave_den+SP[TID+1].ave_den+SP[TID+ngy+1].ave_den);
    Eden = 0.25*SP[TID].ave_den;
    RID = TnRct * TID;
    N1 = FG[TID].n;
    N2 = FG[Gsize + TID].n;
    N4 = FG[3*Gsize + TID].n;
    N5 = FG[4*Gsize + TID].n;
    MCCR[RID + 68] = CollP[68].RR*Eden*N1;
    MCCR[RID + 71] = CollP[71].RR*N4*N5;
    MCCR[RID + 72] = CollP[72].RR*N5*BG[TID].BackDen2;
    MCCR[RID + 73] = CollP[73].RR*N5*BG[TID].BackDen2;
    MCCR[RID + 74] = CollP[74].RR*N5*BG[TID].BackDen2;
    MCCR[RID + 75] = CollP[75].RR*N2*N4;
    MCCR[RID + 76] = CollP[76].RR*N2*BG[TID].BackDen2;
    MCCR[RID + 77] = CollP[77].RR*N2*N2;
    MCCR[RID + 78] = CollP[78].RR*FG[2*Gsize + TID].n*BG[TID].BackDen2;
    MCCR[RID + 82] = CollP[82].RR*N1*BG[TID].BackDen1;
    MCCR[RID + 83] = CollP[83].RR*N1*N4;
    MCCR[RID + 84] = CollP[84].RR*N1*N4;
    MCCR[RID + 85] = CollP[85].RR*N1*BG[TID].BackDen2;
    MCCR[RID + 86] = CollP[86].RR*N1*BG[TID].BackDen2;
    MCCR[RID + 87] = CollP[87].RR*N1*BG[TID].BackDen2;

    Sum_ARM = MCCR[RID + 2];
    Diff_ARM = MCCR[RID + 4] + MCCR[RID + 68] + 2*MCCR[RID + 81] + MCCR[RID + 82]
                + MCCR[RID + 83]+ MCCR[RID + 84]+ MCCR[RID + 85]
                + MCCR[RID + 86] + MCCR[RID + 87];
    Sum_O2A  = MCCR[RID + 8] + MCCR[RID + 73];
    Diff_O2A = MCCR[RID + 19] + MCCR[RID + 20] + MCCR[RID + 21] + MCCR[RID + 22] + MCCR[RID + 23]
                + MCCR[RID + 24] + MCCR[RID + 25] + MCCR[RID + 26] + MCCR[RID + 51] + MCCR[RID + 56]
                 + MCCR[RID + 30] + MCCR[RID + 75] + MCCR[RID + 76] + 2*MCCR[RID + 77];
    Sum_O2B = MCCR[RID + 9] + MCCR[RID + 74];
    Diff_O2B = MCCR[RID + 27] + MCCR[RID + 28] + MCCR[RID + 29] + MCCR[RID + 30] + MCCR[RID + 31]
                 + MCCR[RID + 32] + MCCR[RID + 33] + MCCR[RID + 34] + MCCR[RID + 57] + MCCR[RID + 64]
                 + MCCR[RID + 70] + MCCR[RID + 78];
    Sum_OP = MCCR[RID + 11] + 2 * MCCR[RID + 12] + MCCR[RID + 13] + MCCR[RID + 16] + MCCR[RID + 18]
                + MCCR[RID + 20] + 2*MCCR[RID + 23] + MCCR[RID + 24] + MCCR[RID + 26] + MCCR[RID + 28]
                 + 2*MCCR[RID + 31] + MCCR[RID + 32] + MCCR[RID + 34] + MCCR[RID + 35] + MCCR[RID + 36]
                  + MCCR[RID + 45] + MCCR[RID + 47] + MCCR[RID + 49] + 2*MCCR[RID + 50] + MCCR[RID + 51]
                   + MCCR[RID + 55] + MCCR[RID + 60] + MCCR[RID + 63] + MCCR[RID + 64] + MCCR[RID + 70]
                   + MCCR[RID + 71] + MCCR[RID + 72] + MCCR[RID + 73] + MCCR[RID + 74] + MCCR[RID + 69]
                   +2*MCCR[RID + 85] + MCCR[RID + 86];
    Diff_OP =  MCCR[RID + 38] + MCCR[RID + 39] + MCCR[RID + 40] + MCCR[RID + 41] + MCCR[RID + 42]
                 + MCCR[RID + 43] + MCCR[RID + 48] + MCCR[RID + 52] + MCCR[RID + 79] + MCCR[RID + 83];
    Sum_OD = MCCR[RID + 13] + 2*MCCR[RID + 14] + MCCR[RID + 24] + 2*MCCR[RID + 25] + MCCR[RID + 32]
                 + 2*MCCR[RID + 33] + MCCR[RID + 36] + MCCR[RID + 38] + MCCR[RID + 83] + MCCR[RID + 86];
    Diff_OD = MCCR[RID + 44] + MCCR[RID + 45] + MCCR[RID + 71]+ MCCR[RID + 72]+ MCCR[RID + 73]
                    + MCCR[RID + 74];
    FG_S[TID].n = Sum_ARM - Diff_ARM;
    FG_S[Gsize + TID].n = Sum_O2A - Diff_O2A;
    FG_S[2*Gsize + TID].n = Sum_O2B - Diff_O2B;
    FG_S[3*Gsize + TID].n = Sum_OP - Diff_OP;
    FG_S[4*Gsize + TID].n = Sum_OD - Diff_OD;
}
void Set_Diagnostic_cuda(){
    // Host BUF VECTOR
    Host_G_buf = VFMalloc(Gsize);
    Host_C_buf = VFMalloc(Csize);
    VFInit(Host_G_buf,0.0,Gsize);
    VFInit(Host_C_buf,0.0,Csize);

    checkCudaErrors(cudaMalloc((void**)&Host_G_buf_tmp, sizeof(float) * Gsize));
    checkCudaErrors(cudaMalloc((void**) &Surf_charge_tmp, CondNUMR * sizeof(float)));

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
	
