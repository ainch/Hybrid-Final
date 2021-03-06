#include "cuda_Move.cuh"

void Move_cuda() {
    MoveE_Basic<<<MOVE_GRID, MOVE_BLOCK>>>(Gsize, ngy, dt_dx, dt_dy, dev_info_sp, dev_sp, dev_G_sp, dev_GvecSet);
	//MoveE_Gas_Basic<<<MOVE_GRID, MOVE_BLOCK>>>(Gsize, ngy, dt_dx, dt_dy, dev_info_sp, dev_sp, dev_G_sp, dev_GvecSet);
	cudaDeviceSynchronize();
}
__global__ void MoveE_Gas_Basic(int Gsize,int ngy, float dt_dx,float dt_dy, Species *info, GCP *sp, GPG *data, GGA *Field){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize) return;
	Move_Electron(TID, Gsize,ngy, dt_dx,dt_dy, info, sp, data, Field);
	Move_ion(TID, Gsize,ngy, dt_dx,dt_dy, info, sp, data, Field);
}
__device__ void Move_Electron(int TID, int Gsize,int ngy, float dt_dx,float dt_dy, Species *info, GCP *sp, GPG *data, GGA *Field){
	int PNC;
    PNC = data[TID].PtNumInCell;
    if(PNC==0){
        data[TID].PtNumMoveInterCell=0;
        return;
    } 
	float ex_ws,ex_wn,ex_es,ex_en;
	float ey_ws,ey_wn,ey_es,ey_en;
	float lx,ly;
	int MPNC;
	int PNMC,index;
	int i,k,xp;
	float mvX,mvY,mvZ;
	float del_vx,del_vy;
    float id_cell;
	PNMC=0;
	ex_ws=Field[TID].Ex; 
	ex_wn=Field[TID+1].Ex; 
	ex_es=Field[TID+ngy].Ex; 
	ex_en=Field[TID+ngy+1].Ex;
	ey_ws=Field[TID].Ey; 
	ey_wn=Field[TID+1].Ey; 
	ey_es=Field[TID+ngy].Ey; 
	ey_en=Field[TID+ngy+1].Ey;

	MPNC = data[TID].MaxPtNumInCell;
	i = TID;
	for(k=0;k<PNC;k++){
		lx=sp[i].x; 
		ly=sp[i].y;
		del_vx=ex_ws*(1-lx)*(1-ly)+ex_wn*(1-lx)*ly+ex_es*lx*(1-ly)+ex_en*lx*ly;
		del_vy=ey_ws*(1-lx)*(1-ly)+ey_wn*(1-lx)*ly+ey_es*lx*(1-ly)+ey_en*lx*ly;

		mvX=sp[i].vx+del_vx*info[0].Ascale;
		mvY=sp[i].vy+del_vy*info[0].Ascale;
		mvZ=sp[i].vz;

		lx+=mvX*dt_dx;
		ly+=mvY*dt_dy;

        if(ly>=1 || ly<0 || lx>=1 || lx<0){ // out of cell
            PNMC++;
            index = TID + (MPNC-PNMC)*Gsize;
			id_cell = 0.0f;
            if(ly>=1){								//top
				id_cell+=1;
				ly-=1.0;
			}
			else if(ly<0){							//bottom
				id_cell-=1;
				ly+=1.0;
			}
			if(lx>=1){								//right
				id_cell+=ngy;
				lx-=1.0;
			}
			else if(lx<0){							//left
				id_cell-=ngy;
				lx+=1.0;
			}
			//while(ly>=1 || ly<0 || lx>=1 || lx<0){
			//	if(ly>=1) ly-=1.0;
			//	else if(ly<0) ly+=1.0;
			//	if(lx>=1) lx-=1.0;
			//	else if(lx<0) lx+=1.0;
			//}
            sp[index].CellID = id_cell;
        }else{
            index = i-PNMC*Gsize;
        }
        sp[index].vx=mvX;
		sp[index].vy=mvY;
		sp[index].vz=mvZ;
        sp[index].x=lx;
		sp[index].y=ly;
		//if(sp[index].vx==0) printf("vx[%d]: B[%g]->A[%g]\n",isp,sp[i].x,sp[index].x);
		//if(sp[index].x>1.0f) printf("[%d]=[%d,%d]isp[%d][%g]E[%g,%g]: B[%g]->A[%g] = %g\n",sp[index].CellID,(int)(pp/ngy),(int)(pp%ngy),isp,info[isp].Ascale,del_vx,del_vy,sp[i].vx,sp[index].vx, mvX*dt_dx);
       	i+=Gsize;
    }
	data[TID].PtNumMoveInterCell=PNMC;
	data[TID].PtNumInCell-=PNMC;
}
__device__ void Move_ion(int TID, int Gsize,int ngy, float dt_dx,float dt_dy, Species *info, GCP *sp, GPG *data, GGA *Field){
	int ID = Gsize + TID;
	int PNC;
    PNC = data[ID].PtNumInCell;
    if(PNC==0){
        data[ID].PtNumMoveInterCell=0;
        return;
    } 
	float ex_ws,ex_wn,ex_es,ex_en;
	float ey_ws,ey_wn,ey_es,ey_en;
	float lx,ly;
	int MPNC;
	int PNMC,index;
	int i,k,xp;
	float mvX,mvY,mvZ;
	float del_vx,del_vy;
    float id_cell;
	PNMC=0;
	ex_ws=Field[TID].Ex; 
	ex_wn=Field[TID+1].Ex; 
	ex_es=Field[TID+ngy].Ex; 
	ex_en=Field[TID+ngy+1].Ex;
	ey_ws=Field[TID].Ey; 
	ey_wn=Field[TID+1].Ey; 
	ey_es=Field[TID+ngy].Ey; 
	ey_en=Field[TID+ngy+1].Ey;

	MPNC = data[ID].MaxPtNumInCell;
	i = info[1].St_num + TID;
	for(k=0;k<PNC;k++){
		lx=sp[i].x; 
		ly=sp[i].y;
		del_vx=ex_ws*(1-lx)*(1-ly)+ex_wn*(1-lx)*ly+ex_es*lx*(1-ly)+ex_en*lx*ly;
		del_vy=ey_ws*(1-lx)*(1-ly)+ey_wn*(1-lx)*ly+ey_es*lx*(1-ly)+ey_en*lx*ly;

		mvX=sp[i].vx+del_vx*info[1].Ascale;
		mvY=sp[i].vy+del_vy*info[1].Ascale;
		mvZ=sp[i].vz;

		lx+=mvX*dt_dx;
		ly+=mvY*dt_dy;

        if(ly>=1 || ly<0 || lx>=1 || lx<0){ // out of cell
            PNMC++;
            index = info[1].St_num + TID + (MPNC-PNMC)*Gsize;
			id_cell = 0.0f;
            if(ly>=1){								//top
				id_cell+=1;
				ly-=1.0;
			}
			else if(ly<0){							//bottom
				id_cell-=1;
				ly+=1.0;
			}
			if(lx>=1){								//right
				id_cell+=ngy;
				lx-=1.0;
			}
			else if(lx<0){							//left
				id_cell-=ngy;
				lx+=1.0;
			}
			while(ly>=1 || ly<0 || lx>=1 || lx<0){
				if(ly>=1) ly-=1.0;
				else if(ly<0) ly+=1.0;
				if(lx>=1) lx-=1.0;
				else if(lx<0) lx+=1.0;
			}
            sp[index].CellID = id_cell;
        }else{
            index = i-PNMC*Gsize;
        }
        sp[index].vx=mvX;
		sp[index].vy=mvY;
		sp[index].vz=mvZ;
        sp[index].x=lx;
		sp[index].y=ly;
		//if(sp[index].vx==0) printf("vx[%d]: B[%g]->A[%g]\n",isp,sp[i].x,sp[index].x);
		//if(sp[index].x>1.0f) printf("[%d]=[%d,%d]isp[%d][%g]E[%g,%g]: B[%g]->A[%g] = %g\n",sp[index].CellID,(int)(pp/ngy),(int)(pp%ngy),isp,info[isp].Ascale,del_vx,del_vy,sp[i].vx,sp[index].vx, mvX*dt_dx);
       	i+=Gsize;
    }
	data[ID].PtNumMoveInterCell=PNMC;
	data[ID].PtNumInCell-=PNMC;
}
__global__ void MoveE_Basic(int Gsize,int ngy, float dt_dx,float dt_dy, Species *info, GCP *sp, GPG *data, GGA *Field){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=Gsize*info[0].spnum) return;
	int PNC,isp,ID;
	isp = (int)TID/Gsize; //species number [< nsp]
    ID = (int)TID%Gsize; // Grid ID [< Gsize]
    PNC = data[TID].PtNumInCell;
    if(PNC==0){
        data[TID].PtNumMoveInterCell=0;
        return;
    } 
	float ex_ws,ex_wn,ex_es,ex_en;
	float ey_ws,ey_wn,ey_es,ey_en;
	float lx,ly;
	int MPNC;
	int PNMC,index;
	int i,k,xp;
	float mvX,mvY,mvZ;
	float del_vx,del_vy;
    float id_cell;
	
	PNMC=0;
	
	ex_ws=Field[ID].Ex; 
	ex_wn=Field[ID+1].Ex; 
	ex_es=Field[ID+ngy].Ex; 
	ex_en=Field[ID+ngy+1].Ex;
	ey_ws=Field[ID].Ey; 
	ey_wn=Field[ID+1].Ey; 
	ey_es=Field[ID+ngy].Ey; 
	ey_en=Field[ID+ngy+1].Ey;

	MPNC = data[TID].MaxPtNumInCell;
	i = info[isp].St_num + ID;
	for(k=0;k<PNC;k++){
		lx=sp[i].x; 
		ly=sp[i].y;
		del_vx=ex_ws*(1-lx)*(1-ly)+ex_wn*(1-lx)*ly+ex_es*lx*(1-ly)+ex_en*lx*ly;
		del_vy=ey_ws*(1-lx)*(1-ly)+ey_wn*(1-lx)*ly+ey_es*lx*(1-ly)+ey_en*lx*ly;
		mvX=sp[i].vx+del_vx*info[isp].Ascale;
		mvY=sp[i].vy+del_vy*info[isp].Ascale;
		mvZ=sp[i].vz;

		lx+=mvX*dt_dx;
		ly+=mvY*dt_dy;

        if(ly>=1 || ly<0 || lx>=1 || lx<0){ // out of cell
            PNMC++;
            index = info[isp].St_num + ID + (MPNC-PNMC)*Gsize;
			id_cell = 0.0f;
            if(ly>=1){								//top
				id_cell+=1;
				ly-=1.0;
			}
			else if(ly<0){							//bottom
				id_cell-=1;
				ly+=1.0;
			}
			if(lx>=1){								//right
				id_cell+=ngy;
				lx-=1.0;
			}
			else if(lx<0){							//left
				id_cell-=ngy;
				lx+=1.0;
			}
			/*
			while(ly>=1 || ly<0 || lx>=1 || lx<0){
				if(ly>=1) ly-=1.0;
				else if(ly<0) ly+=1.0;
				if(lx>=1) lx-=1.0;
				else if(lx<0) lx+=1.0;
			}
			*/
            sp[index].CellID = id_cell;
        }else{
            index = i-PNMC*Gsize;
        }
		sp[index].x=lx;
		sp[index].y=ly;
        sp[index].vx=mvX;
		sp[index].vy=mvY;
		sp[index].vz=mvZ;
		//if(sp[index].vx==0) printf("vx[%d]: B[%g]->A[%g]\n",isp,sp[i].x,sp[index].x);
		//if(sp[index].x>1.0f) printf("[%d]=[%d,%d]isp[%d][%g]E[%g,%g]: B[%g]->A[%g] = %g\n",sp[index].CellID,(int)(pp/ngy),(int)(pp%ngy),isp,info[isp].Ascale,del_vx,del_vy,sp[i].vx,sp[index].vx, mvX*dt_dx);
       	i+=Gsize;
    }
	data[TID].PtNumMoveInterCell=PNMC;
	data[TID].PtNumInCell-=PNMC;
}
