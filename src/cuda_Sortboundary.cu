#include "cuda_Sortboundary.cuh"


void SortBounndary_cuda(){
	SortBoundary_Basic<<<SORT_GRID, SORT_BLOCK>>>(Gsize,ngy,CondNUMR,dt_dx,dt_dy, dev_StructureIndex, dev_info_sp, dev_sp, dev_G_sp, dev_GvecSet, dev_CondVec, dev_ReArgFlag);
	cudaDeviceSynchronize();
}
void Set_SortBoundary_cuda(){
	int size;
	size = (ncx + 2) * (ncy + 2);
	checkCudaErrors(cudaMalloc((void**) &dev_StructureIndex, size * sizeof(int)));
    checkCudaErrors(cudaMemset((void *) dev_StructureIndex, 0.0, size * sizeof(int)));
	checkCudaErrors(cudaMemcpy(dev_StructureIndex, vec_StructureIndex , size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMalloc((void**) &dev_ReArgFlag, sizeof(int)));
	checkCudaErrors(cudaMemset((void *) dev_ReArgFlag, 0, sizeof(int)));
}	
__global__ void SortBoundary_Basic(int Gsize,int ngy, int CondNum, float dt_dx,float dt_dy,int *StructureIndex, Species *info, GCP *sp, GPG *data, GGA *Field, GCondA *Cond, int *ReArgFlag){
	int TID = threadIdx.x + blockIdx.x * blockDim.x;
	int isp,ID;
    if(TID>=Gsize*info[0].spnum) return;
    isp = (int)TID/Gsize; //species number [< nsp]
    ID = (int)TID%Gsize; // Grid ID [< Gsize]
	
	int MPNC,PNMIC,oldPNC;
	int CID,MCID,SMCID,DID;
	int PNC;
	int SIndex;
	int index,i,x,k;
	int del_pa,del_pb;
	float delta,del_a,del_b;
	int Flag_type,Flag_x,Flag_y,Flag_vx,Flag_vy;
	int Sum_CID;

	PNMIC = data[TID].PtNumMoveInterCell;
	if(PNMIC==0) return;

	MPNC = data[TID].MaxPtNumInCell;
	PNC = data[TID].PtNumInCell;
	if(PNC+PNMIC>0.9*MPNC) {
		*ReArgFlag=1;
	}
	x = ID/ngy;
	i = info[isp].St_num + ID + (MPNC-1)*Gsize;
	for(k=0;k<PNMIC;k++){
		// Where do the particles go?
		MCID = sp[i].CellID;
		if(MCID>1) // ngy, ngy+1, ngy-1
			SMCID = (ID + ngy + 2 + x) + (MCID+1);	
		else if(MCID<-1) // -ngy, -ngy+1, -ngy-1
			SMCID = (ID + ngy + 2 + x) + (MCID-1);
		else // +1, -1
			SMCID = (ID + ngy + 2 + x) + (MCID);	
		
		CID = TID + MCID; // Destination
		SIndex=StructureIndex[SMCID]; // Destination information
		// 0: PLASMA
		// 1 ~ : Dielectric >> del_a, del_b
		// 100 ~ : Conductor >> charge
		// -2, -4, -6, -8 : Neumann >>  position velocity update
		// -12, -16, -20, -24 : Neumann edge >>  position velocity update

		del_pa = 0; del_pb = 0;
		del_a = 0.0; del_b = 0.0;

		Flag_type = 1;  
		if(SIndex==0){ // Plasma
			Sum_CID = 0;
			Flag_x = 1;	Flag_y = 1;	Flag_vx = 1; Flag_vy = 1;
		}else if(SIndex>=1 && SIndex<100){ // Dielectric
			Flag_type = 0; 
			if(MCID == 1){ // top
				delta = sp[i].x - sp[i].vx/sp[i].vy * sp[i].y;
 				del_pa = CID;		del_a = 1-delta;
				del_pb = CID+ngy;	del_b = delta;
			}else if(MCID == -1){ //bottom
				delta = sp[i].x+sp[i].vx/sp[i].vy*(1-sp[i].y);
				del_pa = TID;		del_a = 1 - delta;
				del_pb = TID+ngy;	del_b = delta;
			}else if(MCID==ngy){		 	// right
				delta = sp[i].y-sp[i].vy/sp[i].vx*sp[i].x;
				del_pa = CID;		del_a = 1 - delta;
				del_pb = CID+1;		del_b = delta;
			}else if(MCID==-ngy){	    // left
				delta = sp[i].y+sp[i].vy/sp[i].vx*(1-sp[i].x);
				del_pa = TID;		del_a = 1 - delta;
				del_pb = TID+1;		del_b = delta;
			}else if(MCID==ngy+1){
				delta = sp[i].x-sp[i].vx/sp[i].vy*sp[i].y;
				if(delta>=0) {
					if(StructureIndex[SMCID-1]) { // right surf --> right cell  --> ngy+1
						delta = sp[i].y-sp[i].vy/sp[i].vx*sp[i].x;
						del_pa = CID - 1;	del_a = -delta;
						del_pb = CID;		del_b = 1 + delta;						
					}else {
						del_pa = CID;		del_a = 1 - delta;
						del_pb = CID+ngy;	del_b = delta;		
					}
				}else { 							// top surf --> top cell --> ngy+1
					if(StructureIndex[SMCID - ngy - 1]) {
						del_pa = TID + 1;	del_a = -delta;
						del_pb = CID;		del_b = 1 + delta;	
					}else {
						delta = sp[i].y-sp[i].vy/sp[i].vx*sp[i].x;
						del_pa = CID;		del_a = 1 - delta;
						del_pb = CID+1;		del_b = delta;	
					}
				}
			}
			else if(MCID==ngy-1){
				delta = sp[i].x+sp[i].vx/sp[i].vy*(1-sp[i].y);
				if(delta>=0) { // right surf --> right cell  --> ngy-1
					if(StructureIndex[SMCID+1]) {
						delta = sp[i].y-sp[i].vy/sp[i].vx*sp[i].x;
						del_pa = CID + 1;	del_a = 2 - delta;
						del_pb = CID + 2;	del_b = delta - 1;	
					}
					else {
						del_pa = CID + 1;	del_a = 1 - delta;
						del_pb = CID+1+ngy;	del_b = delta;	
					}
				}
				else {		 // Bottom surf --> Bottom cell  --> ngy-1
					if(StructureIndex[SMCID-ngy-1]) {
						del_pa = TID;		del_a = -delta;
						del_pb = CID + 1;	del_b = 1 + delta;	
					}
					else {
						delta = sp[i].y-sp[i].vy/sp[i].vx*sp[i].x;
						del_pa = CID;		del_a = 1 - delta;
						del_pb = CID+1;		del_b = delta;	
					}
				}
			}
			else if(MCID==-ngy+1){
				delta = sp[i].x-sp[i].vx/sp[i].vy*sp[i].y;
				if(delta<1) { // Left surf --> Left cell  --> -ngy+1
					if(StructureIndex[SMCID-1]) {
						delta = sp[i].y+sp[i].vy/sp[i].vx*(1-sp[i].x);
						del_pa = CID+ngy-1;	del_a = -delta;
						del_pb = CID+ngy;	del_b = 1 + delta;	
					}
					else {
						del_pa = CID;		del_a = 1-delta;
						del_pb = CID+ngy;	del_b = delta;	
					}
				}
				else {  // top surf --> top cell  --> -ngy+1
					if(StructureIndex[SMCID+ngy+1]) {
						del_pa = TID + 1;	del_a = 2-delta;
						del_pb = TID+1+ngy;	del_b = delta-1;	
					}
					else {
						delta = sp[i].y+sp[i].vy/sp[i].vx*(1-sp[i].x);
						del_pa = TID+1;		del_a = 1-delta;
						del_pb = TID+2;		del_b = delta;
					}
				}
			}
			else if(MCID==-ngy-1){
				delta = sp[i].x+sp[i].vx/sp[i].vy*(1-sp[i].y);
				if(delta<1) { // Left surf --> Left cell  --> -ngy-11
					if(StructureIndex[SMCID+1]) {
						delta = sp[i].y+sp[i].vy/sp[i].vx*(1-sp[i].x);
						del_pa = TID;		del_a = 2-delta;
						del_pb = TID+1;		del_b = delta-1;
					}
					else {
						del_pa = CID;		del_a = 1-delta;
						del_pb = TID;		del_b = delta;	
					}
				}
				else { // Bottom surf --> Bottom cell  --> -ngy-1
					if(StructureIndex[SMCID+ngy+1]) {
						del_pa = TID;		del_a = 2 - delta;
						del_pb = TID+ngy;	del_b = delta - 1;	
					}
					else {
						delta = sp[i].y+sp[i].vy/sp[i].vx*(1-sp[i].x);
						del_pa = TID-1;		del_a = 1 - delta;
						del_pb = TID;		del_b = delta;	
					}
				}
			}
		}else if(SIndex>=100){ // Conductor
			Flag_type = 2;
			index = SIndex - 101;
		}else if(SIndex==-1){ // Dirichlet B.C
			Flag_type = 3; 
		}else if(SIndex==-2){ // Neumann B.C 
			Sum_CID = ngy;
			Flag_x = 0;	Flag_y = 1;	Flag_vx = -1; Flag_vy = 1;
		}else if(SIndex==-6){ // Neumann B.C 
			Sum_CID = -ngy;
			Flag_x = 0;	Flag_y = 1;	Flag_vx = -1; Flag_vy = 1;
		}else if(SIndex==-4){ // Neumann B.C 
			Sum_CID = 1;
			Flag_x = 1;	Flag_y = 0;	Flag_vx = 1; Flag_vy = -1;
		}else if(SIndex==-8){ // Neumann B.C 
			Sum_CID = -1;
			Flag_x = 1;	Flag_y = 0;	Flag_vx = 1; Flag_vy = -1;
		}else if(SIndex==-12){ // Neumann B.C 
			Sum_CID = ngy+1;
			Flag_x = 0;	Flag_y = 0;	Flag_vx = -1; Flag_vy = -1;
		}else if(SIndex==-16){ // Neumann B.C 
			Sum_CID = 1-ngy;
			Flag_x = 0;	Flag_y = 0;	Flag_vx = -1; Flag_vy = -1;
		}else if(SIndex==-20){ // Neumann B.C 
			Sum_CID = ngy-1;
			Flag_x = 0;	Flag_y = 0;	Flag_vx = -1; Flag_vy = -1;
		}else if(SIndex==-24){ // Neumann B.C 
			Sum_CID = -ngy-1;
			Flag_x = 0;	Flag_y = 0;	Flag_vx = -1; Flag_vy = -1;
		}
		if(Flag_type==0){
			atomicAdd(&data[del_pa].sigma,del_a);
			atomicAdd(&data[del_pb].sigma,del_b);
		}else if(Flag_type==1){
			CID+=Sum_CID;
			oldPNC = atomicAdd(&data[CID].PtNumInCell,1);
			DID = (int)CID%Gsize;
			index = info[isp].St_num + DID + oldPNC * Gsize; // ??
			sp[index].CellID = CID;
			if(Flag_x==1) sp[index].x = sp[i].x;
			else sp[index].x = 1-sp[i].x;
			if(Flag_y==1) sp[index].y = sp[i].y;
			else sp[index].y = 1-sp[i].y;
			sp[index].vx = Flag_vx * sp[i].vx;
			sp[index].vy = Flag_vy * sp[i].vy;
			sp[index].vz = sp[i].vz;
		}else if(Flag_type==2){
			atomicAdd(&Cond[isp*CondNum + index].Charge,1.0);
		}
		i-=Gsize;
	}		
}