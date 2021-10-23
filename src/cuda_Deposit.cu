#include "cuda_Deposit.cuh"

void Deposit_cuda(){
    int i;
    // DepositAtom - Particle > density
	DepositInitDensity<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(Gsize,dev_info_sp,dev_G_sp);
	//start plot Gsize values
	cudaDeviceSynchronize();
    DepositAtom<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(Gsize,ngy,dev_info_sp,dev_sp,dev_G_sp);
	cudaDeviceSynchronize();
    DepositBoundary<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(Gsize,ngy,nsp,dev_GvecSet,dev_G_sp);
    // Smoothing
    for(i=0;i<N_smt;i++){
        Smooth_121_A<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(ngx, ngy, nsp, dev_GvecSet, dev_G_sp);
        Smooth_121_B<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(ngx, ngy, nsp, dev_GvecSet, dev_G_sp);
    }
    // charge density cal --> dev_Source
    SumSource<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(ngx, ngy, dev_info_sp, dev_GvecSet, dev_G_sp, dev_Sigma, dev_Source);
    // dev_Source --> dev_b  (PCG setting)
    PCG_Set<<<DEPOSIT_GRID,DEPOSIT_BLOCK>>>(Gsize,dev_A_idx, dev_GvecSet, dev_Source, dev_R);
	cudaDeviceSynchronize();
}
__global__ void PCG_Set(int Gsize, int *IDX, GGA *vecSet, float *Source, float *B){
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
	if(TID>=Gsize) return;
	if(IDX[TID]) B[IDX[TID]-1]=Source[TID];
}
__global__ void SumSource(int ngx, int ngy, Species *info, GGA *vecSet, GPG *data, float *Sigma, float *Source){
     int TID = threadIdx.x + blockIdx.x * blockDim.x;
     if(TID>=ngx*ngy) return;
     int i;
     int x,y;
     x=TID/ngy; y=TID%ngy;
	float sum, SD;
    float SumSig = 0.0;
     
	if(x==0) {
		if(vecSet[TID+ngy].Boundary==DIELECTRIC) {
            for(i=0;i<info[0].spnum;i++) SumSig += data[TID + i*ngx*ngy].sigma * info[i].q_density;
			Sigma[TID] = 2 * SumSig/vecSet[TID].Area;
		    SD = SumSig;
		}
		else if(vecSet[TID+ngy].Boundary==CONDUCTOR) {
               for(i=0;i<info[0].spnum;i++) data[TID + i*ngx*ngy].sigma = 0;
			SD=0;
		}
	}
	else if(x==ngx-1) {
		if(vecSet[TID-ngy].Boundary==DIELECTRIC) {
			for(i=0;i<info[0].spnum;i++) SumSig += data[TID + i*ngx*ngy].sigma * info[i].q_density;
			Sigma[TID] = 2 * SumSig/vecSet[TID].Area;
		    SD = SumSig;
		}
		else if(vecSet[TID-ngy].Boundary==CONDUCTOR) {
			for(i=0;i<info[0].spnum;i++) data[TID + i*ngx*ngy].sigma = 0;
			SD=0;
		}
	}
	else if(y==0) {
		if(vecSet[TID+1].Boundary==DIELECTRIC) {
			for(i=0;i<info[0].spnum;i++) SumSig += data[TID + i*ngx*ngy].sigma * info[i].q_density;
			Sigma[TID] = 2 * SumSig/vecSet[TID].Area;
		    SD = SumSig;
		}
		else if(vecSet[TID+1].Boundary==CONDUCTOR) {
			for(i=0;i<info[0].spnum;i++) data[TID + i*ngx*ngy].sigma = 0;
			SD=0;
		}
	}
	else if(y==ngy-1) {
		if(vecSet[TID-1].Boundary==DIELECTRIC) {
			for(i=0;i<info[0].spnum;i++) SumSig += data[TID + i*ngx*ngy].sigma * info[i].q_density;
			Sigma[TID] = 2 * SumSig/vecSet[TID].Area;
		    SD = SumSig;
		}
		else if(vecSet[TID-1].Boundary==CONDUCTOR) {
			for(i=0;i<info[0].spnum;i++) data[TID + i*ngx*ngy].sigma = 0;
			SD=0;
		}
	}
	else if(vecSet[TID].Boundary==DIELECTRIC) {
		for(i=0;i<info[0].spnum;i++) SumSig += data[TID + i*ngx*ngy].sigma * info[i].q_density;
		Sigma[TID] = SumSig/vecSet[TID].Area;
		SD = SumSig;
	}
	else if(vecSet[TID].Boundary==CONDUCTOR) {
		for(i=0;i<info[0].spnum;i++) data[TID + i*ngx*ngy].sigma = 0;
		SD=0;
	}
	else {
		SD=0;
	}
    SumSig = 0.0;
    for(i=0;i<info[0].spnum;i++) SumSig += data[TID + i*ngx*ngy].den * info[i].q_density;
	sum=(SumSig + SD)/EPS0;
	Source[TID]=sum;
}
__global__ void Smooth_121_A(int ngx, int ngy, int nsp, GGA *vecSet, GPG *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=ngx*ngy*nsp) return;
    if(data[TID].den == 0) return;
    int ID;
    int x,y;
    ID = (int)TID%(ngx*ngy);
    x=ID/ngy; y=ID%ngy;
    // den >> smt_den at x direction
     if(vecSet[ID].Boundary==0){
        data[TID].smt_den = 0.25 * (data[TID-ngy].den + 2  *data[TID].den + data[TID+ngy].den);
	}else if(vecSet[ID].Boundary==CONDUCTOR || vecSet[ID].Boundary==DIELECTRIC) {
		if(vecSet[ID].Face==LEFT || vecSet[ID].Face==UL_CORN || vecSet[ID].Face==LL_CORN) 
               data[TID].smt_den = 0.5 * (data[TID].den + data[TID-ngy].den);
		else if(vecSet[ID].Face==RIGHT || vecSet[ID].Face==UR_CORN || vecSet[ID].Face==LR_CORN) 
               data[TID].smt_den = 0.5 * (data[TID].den + data[TID+ngy].den);
		else if(vecSet[ID].Face==UP || vecSet[ID].Face==DOWN) 
               data[TID].smt_den = 0.25 * (data[TID-ngy].den + 2*data[TID].den + data[TID+ngy].den);
	}else if(vecSet[ID].Boundary==DIRICHLET) {
		if(vecSet[ID].Face==LEFT || vecSet[ID].Face==UL_CORN || vecSet[ID].Face==LL_CORN) 
               data[TID].smt_den = 0.5 * (data[TID].den + data[TID-ngy].den);
		else if(vecSet[ID].Face==RIGHT || vecSet[ID].Face==UR_CORN || vecSet[ID].Face==LR_CORN) 
               data[TID].smt_den = 0.5 * (data[TID].den + data[TID+ngy].den);
		else if(vecSet[ID].Face==UP || vecSet[ID].Face==DOWN) {
			if(x==0) 
                    data[TID].smt_den = 0.5 * (data[TID].den+data[TID+ngy].den);
			else if(x==ngx-1) 
                    data[TID].smt_den = 0.5 * (data[TID].den+data[TID-ngy].den);
			else if(y==0 || y==ngy-1) 
                    data[TID].smt_den = 0.25 * (data[TID-ngy].den+2*data[TID].den+data[TID+ngy].den);
		}
	}else if(vecSet[ID].Boundary==NEUMANN) {
		if(vecSet[ID].Face==NO_FACE && vecSet[ID].CondID==0) {
			if(x==0) 
                    data[TID].smt_den=0.5*(data[TID].den+data[TID+ngy].den);
			else if(x==ngx-1) 
                    data[TID].smt_den=0.5*(data[TID].den+data[TID-ngy].den);
			else if(y==0 || y==ngy-1) 
                    data[TID].smt_den=0.25*(data[TID-ngy].den+2*data[TID].den+data[TID+ngy].den);
		}else if(vecSet[ID].Face==LEFT || vecSet[ID].Face==UL_CORN || vecSet[ID].Face==LL_CORN) 
               data[TID].smt_den=0.5*(data[TID].den+data[TID-ngy].den);
		else if(vecSet[ID].Face==RIGHT || vecSet[ID].Face==UR_CORN || vecSet[ID].Face==LR_CORN) 
               data[TID].smt_den=0.5*(data[TID].den+data[TID+ngy].den);
		else if(vecSet[ID].Face==UP || vecSet[ID].Face==DOWN) {
			if(x==0) 
                    data[TID].smt_den=0.5*(data[TID].den+data[TID+ngy].den);
			else if(x==ngx-1) 
                    data[TID].smt_den=0.5*(data[TID].den+data[TID-ngy].den);
			else if(y==0 || y==ngy-1) 
                    data[TID].smt_den=0.25*(data[TID-ngy].den+2*data[TID].den+data[TID+ngy].den);
		}
	}
}
__global__ void Smooth_121_B(int ngx, int ngy, int nsp, GGA *vecSet, GPG *data){
     int TID = threadIdx.x + blockIdx.x * blockDim.x;
     if(TID>=ngx*ngy*nsp) return;
     if(data[TID].smt_den == 0) return;
     int ID;
     int x,y;
     ID = (int)TID%(ngx*ngy);
     x=ID/ngy; y=ID%ngy;
     // smt_den >> den at y direction
     if(vecSet[ID].Boundary==0) {
          data[TID].den=0.25*(data[TID-1].smt_den+2*data[TID].smt_den+data[TID+1].smt_den);
	}else if(vecSet[ID].Boundary==CONDUCTOR || vecSet[ID].Boundary==DIELECTRIC) {
		if(vecSet[ID].Face==UP || vecSet[ID].Face==UL_CORN || vecSet[ID].Face==UR_CORN) 
               data[TID].den=0.5*(data[TID].smt_den+data[TID+1].smt_den);
		else if(vecSet[ID].Face==DOWN || vecSet[ID].Face==LR_CORN || vecSet[ID].Face==LL_CORN) 
               data[TID].den=0.5*(data[TID].smt_den+data[TID-1].smt_den);
		else if(vecSet[ID].Face==LEFT || vecSet[ID].Face==RIGHT) 
               data[TID].den=0.25*(data[TID-1].smt_den+2*data[TID].smt_den+data[TID+1].smt_den);
	}else if(vecSet[ID].Boundary==DIRICHLET) {
		if(vecSet[ID].Face==UP || vecSet[ID].Face==UL_CORN || vecSet[ID].Face==UR_CORN) 
               data[TID].den=0.5*(data[TID].smt_den+data[TID+1].smt_den);
		else if(vecSet[ID].Face==DOWN || vecSet[ID].Face==LR_CORN || vecSet[ID].Face==LL_CORN) 
               data[TID].den=0.5*(data[TID].smt_den+data[TID-1].smt_den);
		else if(vecSet[ID].Face==LEFT || vecSet[ID].Face==RIGHT) {
			if(y==0) 
                    data[TID].den=0.5*(data[TID].smt_den+data[TID+1].smt_den);
			else if(y==ngy-1) 
                    data[TID].den=0.5*(data[TID].smt_den+data[TID-1].smt_den);
			else if(x==0 || x==ngx-1) 
                    data[TID].den=0.25*(data[TID-1].smt_den+2*data[TID].smt_den+data[TID+1].smt_den);
		}
     }else if(vecSet[ID].Boundary==NEUMANN) {
		if(vecSet[ID].Face==NO_FACE && vecSet[ID].CondID==0) {
			if(y==0) 
                    data[TID].den=0.5*(data[TID].smt_den+data[TID+1].smt_den);
			else if(y==ngy-1) 
                    data[TID].den=0.5*(data[TID].smt_den+data[TID-1].smt_den);
			else if(x==0 || x==ngx-1) 
                    data[TID].den=0.25*(data[TID-1].smt_den+2*data[TID].smt_den+data[TID+1].smt_den);
		}else if(vecSet[ID].Face==UP || vecSet[ID].Face==UL_CORN || vecSet[ID].Face==UR_CORN) 
               data[TID].den=0.5*(data[TID].smt_den+data[TID+1].smt_den);
		else if(vecSet[ID].Face==DOWN || vecSet[ID].Face==LR_CORN || vecSet[ID].Face==LL_CORN) 
               data[TID].den=0.5*(data[TID].smt_den+data[TID-1].smt_den);
		else if(vecSet[ID].Face==LEFT || vecSet[ID].Face==RIGHT) {
			if(y==0) 
                    data[TID].den=0.5*(data[TID].smt_den+data[TID+1].smt_den);
			else if(y==ngy-1) 
                    data[TID].den=0.5*(data[TID].smt_den+data[TID-1].smt_den);
			else if(x==0 || x==ngx-1) 
                    data[TID].den=0.25*(data[TID-1].smt_den+2*data[TID].smt_den+data[TID+1].smt_den);
		}
	}
}
__global__ void DepositBoundary(int Gsize, int ngy, int nsp, GGA *vecSet, GPG *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize*nsp) return;
    int ID;
    ID = (int)TID%Gsize;

	if(vecSet[ID].Face==LEFT || vecSet[ID].Face==RIGHT || vecSet[ID].Face==UP || vecSet[ID].Face==DOWN) data[TID].den *= 2;
	else if(vecSet[ID].Face==UL_CORN || vecSet[ID].Face==LL_CORN || vecSet[ID].Face==UR_CORN || vecSet[ID].Face==LR_CORN) data[TID].den *= 4;
	if(vecSet[ID].Boundary==NEUMANN) data[TID].den *= 2;
	else if(vecSet[ID].Boundary==DIRICHLET && vecSet[ID].Face!=NO_FACE && vecSet[ID].CondID!=0) data[TID].den *= 2;

}
__global__ void DepositAtom(int Gsize, int ngy, Species *info, GCP *sp, GPG *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
	if(TID>=Gsize*info[0].spnum) return;
	int ID,isp;
	isp = (int)TID/Gsize;
    ID = (int)TID%Gsize;

	int i,k,PNC;
    float lx,ly;
	float WS,WN,ES,EN;
   
    PNC = data[TID].PtNumInCell;
    if(PNC==0) return;
    WS=0.0; WN=0.0; ES=0.0; EN=0.0;
	i = info[isp].St_num + ID;
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
	atomicAdd(&info[isp].np,PNC); 
}
__global__ void DepositInitDensity(int Gsize, Species *info, GPG *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize*info[0].spnum) return;
	int isp = (int)TID/Gsize;
    data[TID].den = 0.0;
	info[isp].np = 0;
}