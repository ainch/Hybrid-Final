#include "cuda_Particle.cuh"

void Set_Particle_cuda(){  
    // This function calculates the following variables :
    // 1. dev_info_sp [Species] - [nsp]
    // 2. dev_sp [GCP] - [nsp * Sum of MAXNP]
    // 3. dev_G_sp [GPG] - [nsp * Gsize]
    int isp;
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
    int isp,k,PID,GID,XID,YID,index;
    int ID,SPID,SGID,PNC;
    float sum = 0.0;
    int count = 0;
    printf("Gsize=%d GID = %d\n",Gsize,Gsize*(nsp));
    // particle data >> buf data
    for (isp = 0; isp < nsp; isp++){
        SPID = info[isp].St_num;
        SGID = isp * Gsize;
        printf(" - %s : Start address = %d\n",info[isp].name,SPID);
        printf(" - %s : Number of PT = %d\n",info[isp].name,info[isp].np);
        for (k = 0; k < info[isp].np; k++) {
            //printf("[%d].np = %d, k = %d\n",isp,info[isp].np,k);
            XID = A[isp].x[k];
		    YID = A[isp].y[k];
            ID = XID * ngy + YID;
            PID = SPID + ID;
            GID = SGID + ID;
            PNC = C[GID].PtNumInCell;
            index = PID + PNC * Gsize;
            B[index].CellID = GID;
		    B[index].x = A[isp].x[k] - XID;
		    B[index].y = A[isp].y[k] - YID;
		    B[index].vx = A[isp].vx[k];
		    B[index].vy = A[isp].vy[k];
		    B[index].vz = A[isp].vz[k];
            C[GID].PtNumInCell++;
            //if(isp ==1 && XID <8) printf("np[%d] = XID[%d]\n",k,XID);
        }
        sum = 0.0;
        count = 0;
        for(k=0;k<Gsize;k++){
            GID = isp*Gsize + k;
            //printf("[%d]:GID=[%d]\n",k,GID);
            C[GID].den = info[isp].Denscale * C[GID].PtNumInCell;
            sum += C[GID].den;
            if( C[GID].den != 0)  count++;
        }
        printf(" - %s : Average Density = %g\n",info[isp].name,sum/count);
    }
}
void Copy_GCPtoHCP(Species *info, GCP *A, HCP *B, GPG *C){
    // INPUT : A[nsp * MAXNP]
    // OUTPUT : B[isp][NP_LIM], C[Gsize*nsp]
    int isp,i,j,k;
    int CID,AID;
    int XID,YID;
    for(isp=0;isp<nsp;isp++){
        k=0;
        for (i = 0; i < Gsize; i++){
            CID = isp * Gsize + i;
            for(j=0;j<C[CID].PtNumInCell;j++){
                AID = info[isp].St_num + i + j * Gsize;
                B[isp].CellID[k] = A[AID].CellID - isp * Gsize;
                //printf("[%d][%d]=[%d]\n",isp,k,B[isp].CellID[k]);
                XID = (int)(B[isp].CellID[k]/ngy);
                YID = (int)(B[isp].CellID[k]%ngy);
                B[isp].x[k] = XID + A[AID].x;
                B[isp].y[k] = YID + A[AID].y;
                B[isp].vx[k] = A[AID].vx;
                B[isp].vy[k] = A[AID].vy;
                B[isp].vz[k] = A[AID].vz;
                //printf("[%d][%d]=[%g][%g][%g][%g][%g]\n",isp,k,B[isp].x[k],B[isp].y[k],B[isp].vx[k],B[isp].vy[k],B[isp].vz[k]);
                k++;
            }
        }
        info[isp].np = k;
    }
    //exit(1);
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