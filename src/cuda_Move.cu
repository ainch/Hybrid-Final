#include "cuda_Move.cuh"

void Move_cuda() {
    MoveE_Basic<<<MOVE_GRID, MOVE_BLOCK>>>(Gsize, ngy, dt_dx, dt_dy, dev_info_sp, dev_sp, dev_G_sp, dev_GvecSet);
	cudaDeviceSynchronize();
}
__global__ void MoveE_Basic(int Gsize,int ngy, float dt_dx,float dt_dy, Species *info, GCP *sp, GPG *data, GGA *Field){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
	int PNC,isp,ID;
	isp = (int)TID/Gsize; //species number [< nsp]
    ID = (int)TID%Gsize; // Grid ID [< Gsize]
    if(TID>Gsize*info[isp].spnum) return;
	
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
            sp[index].CellID = id_cell;
        }else{
            index = i-PNMC*Gsize;
        }
        sp[index].vx=mvX;
		sp[index].vy=mvY;
		sp[index].vz=mvZ;
        sp[index].x=lx;
		sp[index].y=ly;
		if(sp[index].vx==0) printf("B[%g]->A[%g]\n",sp[i].x,sp[index].x);
        i+=Gsize;
    }
	data[TID].PtNumMoveInterCell=PNMC;
	data[TID].PtNumInCell-=PNMC;

}