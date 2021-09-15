#include "cuda_Deposit.cuh"

void Deposit_cuda(){
    int i,j;
    // Function list
    // 1. DepositAtom - Particle > density
    // 2. boundary density up
    // 3. density smoothing
    // 4. charge density cal --> dev_Source
    // 5. dev_Source --> dev_b  (PCG setting)
    // <<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>
    printf("Deposit start.\n");
    DepositAtom<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(Gsize,ngy,dev_info_sp,dev_sp,dev_G_sp);
    DepositBoundary<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(Gsize,ngy,nsp,dev_GvecSet,dev_G_sp);
    // Smoothing

    // 4. charge density cal --> dev_Source

    // 5. dev_Source --> dev_b  (PCG setting)

    // TEST
    cudaMemcpy(Host_G_sp, dev_G_sp, Gsize * nsp * sizeof(GPG), cudaMemcpyDeviceToHost);
    for(j=ngy-1;j>=0;j--){
        for(i=0;i<ngx;i++){
            printf("%8.2g ",Host_G_sp[i*ngy+j].den);
        } printf("\n");
    }

}
__global__ void DepositBoundary(int Gsize, int ngy, int nsp, GGA *vecSet, GPG *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>Gsize*nsp) return;
    int ID;
    ID = (int)TID%Gsize;

	if(vecSet[ID].Face==LEFT || vecSet[ID].Face==RIGHT || vecSet[ID].Face==UP || vecSet[ID].Face==DOWN) data[TID].den *= 2;
	else if(vecSet[ID].Face==UL_CORN || vecSet[ID].Face==LL_CORN || vecSet[ID].Face==UR_CORN || vecSet[ID].Face==LR_CORN) data[TID].den *= 4;
	if(vecSet[ID].Boundary==NEUMANN) data[TID].den *= 2;
	else if(vecSet[ID].Boundary==DIRICHLET && vecSet[ID].Face!=NO_FACE && vecSet[ID].CondID!=0) data[TID].den *= 2;
}
__global__ void DepositAtom(int Gsize, int ngy, Species *info, GCP *sp, GPG *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
	int i,k,PNC,isp,ID;
    float lx,ly;
	float WS,WN,ES,EN;

    if(TID>Gsize*info[0].spnum) return;

    data[TID].den = 0.0;

    PNC = data[TID].PtNumInCell;
    if(PNC==0) return;

    isp = (int)TID/Gsize;
    ID = (int)TID%Gsize;

    //printf("TID[%d] = [%d]\n",TID,isp);
    WS=0.0; WN=0.0; ES=0.0; EN=0.0;
	i=info[isp].St_num + ID;

	for(k=0;k<PNC;k++){
		lx=sp[i].x; ly=sp[i].y;
		WS+=(1-lx)*(1-ly);
		WN+=(1-lx)*ly;
		ES+=lx*(1-ly);
		EN+=lx*ly;
		i+=Gsize;
	}

	atomicAdd(&data[TID].den,WS);
	atomicAdd(&data[TID+1].den,WN);
	atomicAdd(&data[TID+ngy].den,ES);
	atomicAdd(&data[TID+ngy+1].den,EN); 
}