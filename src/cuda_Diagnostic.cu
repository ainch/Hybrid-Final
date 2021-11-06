#include "cuda_Diagnostic.cuh"

void Diagnostic(){
    int i, k, isp;
    // Data Copy
    cudaMemcpy(vec_Potential, TotPotential, Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vec_Source, dev_Source, Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vec_Sigma, dev_Sigma, Gsize * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(vec_G, dev_GvecSet, Gsize * sizeof(GGA), cudaMemcpyDeviceToHost);
    if(nave_count==N_ave){  // average calculate
        //Conti_Flag = 1;  
        Average_Particle_Density<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(nsp, Gsize, N_ave, dev_info_sp, dev_G_sp);
        // field data
        for(i=0;i<Gsize;i++){
            ave_Potential[i] = sum_Potential[i]/nave_count;
            ave_Source[i] += sum_Source[i]/nave_count;
            ave_Sigma[i] += sum_Sigma[i]/nave_count;
            ave_Ex[i] += sum_Ex[i]/nave_count;
            ave_Ey[i] += sum_Ey[i]/nave_count;
            sum_Potential[i] = 0.0f;
            sum_Source[i] = 0.0f;
            sum_Sigma[i] = 0.0f;
            sum_Ex[i] = 0.0f;
            sum_Ey[i] = 0.0f;
            sum_Potential[i] += vec_Potential[i];
            sum_Source[i] += vec_Source[i];
            sum_Sigma[i] += vec_Sigma[i];
            sum_Ex[i] += vec_G[i].Ex;
            sum_Ey[i] += vec_G[i].Ey;
        }
        nave_count = 0;
    }else{ // accomulation data
        Accomulate_Particle_Density<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(nsp, Gsize, dev_G_sp);

        // field data
        for(i=0;i<Gsize;i++){
            sum_Potential[i] += vec_Potential[i];
            sum_Source[i] += vec_Source[i];
            sum_Sigma[i] += vec_Sigma[i];
            sum_Ex[i] += vec_G[i].Ex;
            sum_Ey[i] += vec_G[i].Ey;
        }
        nave_count++;
    }
    if (hist_count >= HISTMAX) {
        for (k = 1, i = 4; k < HISTMAX / 4; k++, i += 4) {
			t_array[k] = t_array[i];
			for (isp = 0; isp < nsp; isp++) {
                HistPt[isp].np[k] =  HistPt[isp].np[i];
			}
			iter_array[k] = iter_array[i];
		}
		hist_count = k;
		dHIST *= 4;
    }
    
    
    if((--Hcount)==0){
        t_array[hist_count] = (float) t;
        iter_array[hist_count] = (float) *FIter;
        cudaMemcpy(SP, dev_info_sp, nsp * sizeof(Species), cudaMemcpyDeviceToHost);
        for (isp = 0; isp < nsp; isp++) {
            HistPt[isp].np[hist_count] = SP[isp].np;
        }
        hist_count++;
        Hcount = dHIST;
    }
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
    float temp = info[isp].Denscale/N_ave;
    data[TID].sum_den += data[TID].den;
    data[TID].ave_den = temp * data[TID].sum_den;
    data[TID].sum_den = 0.0f;
}
void Set_Diagnostic_cuda(){
    // Host BUF VECTOR
    Host_G_buf = VFMalloc(Gsize);
    Host_C_buf = VFMalloc(Csize);
    VFInit(Host_G_buf,0.0,Gsize);
    VFInit(Host_C_buf,0.0,Csize);

}