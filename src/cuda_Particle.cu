#include "cuda_Particle.cuh"

void Set_Particle_cuda(){  
    // This function calculates the following variables :
    // 1. dev_info_sp [Species] - [nsp]
    // 2. dev_sp [GCP] - [nsp * Sum of MAXNP]
    // 3. dev_G_sp [GPG] - [nsp * Gsize]
    int isp,i;
    int Total_maxnp = 0;
    printf(" Particle copy CPU to GPU. --> ");
    for(isp=0;isp<nsp;isp++){
        Total_maxnp += SP[isp].MAXNP;
    }
    //printf(" Total Maxnp = %d\n",Total_maxnp);
    // CPU >> CPU
    Host_sp = (GCP *)malloc(Total_maxnp * sizeof(GCP));
    GCPInit(Total_maxnp,Host_sp);
    Copy_HCPtoGCP(SP, PtD, Host_sp, Host_G_sp);
    // CPU >> GPU 
    checkCudaErrors(cudaMalloc((void**)&dev_info_sp, nsp * sizeof(Species)));
    checkCudaErrors(cudaMemcpy(dev_info_sp, SP, nsp * sizeof(Species), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc((void**)&dev_sp, Total_maxnp * sizeof(GCP)));
    checkCudaErrors(cudaMemcpy(dev_sp, Host_sp, Total_maxnp * sizeof(GCP), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc((void**)&dev_G_sp, Gsize * nsp * sizeof(GPG)));
    checkCudaErrors(cudaMemcpy(dev_G_sp, Host_G_sp, Gsize * nsp * sizeof(GPG), cudaMemcpyHostToDevice));
    printf("Complete!\n");
}
void Copy_HCPtoGCP(Species *info, HCP *A, GCP *B, GPG *C){
    // OUTPUT : B[MAXNP*nsp], C[Gsize*nsp]
    int isp,k,XID,YID,index;
    // particle data >> buf data
    for (isp = 0; isp < nsp; isp++){
	    for (k = 0; k < info[isp].np; k++) {
            XID = A[isp].x[k];
		    YID = A[isp].y[k];
            // index = nsp? + cellID? + ptNumInCell?
            index = info[isp].St_num + A[isp].CellID[k] + 
                        C[isp*Gsize + A[isp].CellID[k]].PtNumInCell * Gsize;
            //if(k == 0) printf("x=%d y=%d id=%d\n",XID,YID,index);
            B[index].CellID = A[isp].CellID[k];
		    B[index].x = A[isp].x[k] - XID;
		    B[index].y = A[isp].y[k] - YID;
		    B[index].vx = A[isp].vx[k];
		    B[index].vy = A[isp].vy[k];
		    B[index].vz = A[isp].vz[k];
            C[isp*Gsize + A[isp].CellID[k]].PtNumInCell++;
            //if(k<5){
                //printf("ISP[%d][%d]:x[%g]y[%g]>ID[%d]:x[%g]\n",isp,k,A[isp].x[k],A[isp].y[k],A[isp].CellID[k],B[index].x);
                //printf("np[%d] = %d\n",isp*Gsize + A[isp].CellID[k],C[isp*Gsize + A[isp].CellID[k]].PtNumInCell);
            //} 
        }
        for(k=0;k<Gsize;k++){
            C[isp*Gsize + A[isp].CellID[k]].den = C[isp*Gsize + A[isp].CellID[k]].PtNumInCell;
            //printf("ISP[%d],ID[%d] = NP[%d], ",isp,k,C[isp*Gsize + k].PtNumInCell);
            //printf("MaxNP[%d]\n",C[isp*Gsize + k].MaxPtNumInCell);
        }
    }
}
void Copy_GCPtoHCP(Species *info, GCP *A,HCP *B, GPG *C){
    // OUTPUT : B[isp][NP_LIM], C[Gsize*nsp]
    int isp,k,XID,YID,index;
    exit(1);
    for (isp = 0; isp < nsp; isp++){
    }
}
void GCPInit(int size, GCP *A){
    int i;
    for (i = 0; i < size; i++){
        A[i].CellID = -1;
        A[i].x = 0.0;
        A[i].y = 0.0;
        A[i].vx = 0.0;
        A[i].vy = 0.0;
        A[i].vz = 0.0;
    }
}