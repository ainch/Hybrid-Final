#include "cuda_mccAr.cuh"
__device__ void Direct_Argon_Electron(int Gsize, int ngy, int ID, int MCCn, float dtm, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, GGA *BG, GFC *Fluid){
	int i,j,k,n,index;
	int CID,PNC,PNC2;
	int nx,ny,ngx;
	int Target,oldPNC;
	int Colltype;
	float Tprob,Prob1,Prob2;
    int Randn,AddPt1;
	float R1,R2;
	float VX,VY,VZ,VX_buf,VY_buf,VZ_buf;
	float dum,vel,vel2,engy,rengy;
	float SumSigma,SumEngyLoss;
	PNC = data[ID].PtNumInCell;
    PNC2 = data[ID+Gsize].PtNumInCell;
	curandState LocalStates = states[ID];
    nx = ID/ngy;
	ny = ID%ngy;
	ngx = Gsize/ngy;
	if(nx == ngx-1) nx--;
	if(ny == ngy-1) ny--;
	CID = ny + (ngy-1)*nx;
	Prob1 = 1.0f - exp(-1*dtm*sigv[0].val*BG[ID].BackDen1);
	Prob2 = 1.0f - exp(-1*dtm*sigv[1].val*Fluid[CID].ave_den);
	Tprob = Prob1 + Prob2;
    // Calculate total Collision probability.
	Randn = MCCn;
    AddPt1 = 0;
	i = info[0].St_num + ID;
	for(k=0;k<PNC;k++){
        Colltype = 0;
		for(j=0;j<Randn;j++){
			R1 = curand_uniform(&LocalStates);
			if(R1<Tprob){
                Colltype = 1;
                break;
            }
		}
        if(Colltype == 0){
            i+=Gsize;
			continue;
        }
		R1 = Tprob * R1;
        if(R1 <= Prob1)	Target = (int)0;
		else			Target = (int)1;
        // Calculate energy
		VX = sp[i].vx;
		VY = sp[i].vy;
		VZ = sp[i].vz;
		dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		VX/=vel; VY/=vel; VZ/=vel;
		engy = info[0].Escale * dum;
        Colltype = 2;
        switch(Target){
			case 0:{
                R2 = curand_uniform(&LocalStates) * sigv[0].val / vel;
				// 0. e + Ar > e + Ar 			Elastic Scattering
				SumSigma = Argon_CrossSection(0, engy, N_LOGX, idLOGX, CX);
				if(R2<=SumSigma){
				// 1. e + Ar > e + Ar* 			Excitation to Total Excited state
				}else if(engy > info_CX[1].Th_e && R2<=(SumSigma += Argon_CrossSection(1, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[1].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
				// 2. e + Ar > e + Ar* 			Excitation to AR4SM
				}else if(engy > info_CX[2].Th_e && R2<=(SumSigma += Argon_CrossSection(2, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[2].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
				// 3. e + Ar > e + e + Ar^		Direct ionization
				}else if(engy > info_CX[3].Th_e && R2<=(SumSigma += Argon_CrossSection(3, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3;
					engy-=info_CX[3].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					AddPt1++;
				}else{
					Colltype = 1;
				}
				break;
			}
			case 1:{
                R2 = curand_uniform(&LocalStates)*sigv[1].val / vel;
				// 4. e + Ar* > e + e + Ar^		step ionization
				SumSigma = Argon_CrossSection(4, engy, N_LOGX, idLOGX, CX);
				if(engy > info_CX[4].Th_e && R2<=SumSigma){
					Colltype = 3;
					engy-=info_CX[4].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					AddPt1++;
				}else{
					Colltype = 1;
				}
				break;
			}
		}
        if(Colltype == 2){ // Just energy loss
			dev_anewvel(engy,vel,&VX,&VY,&VZ,0,info_CX[4].mofM,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			sp[i].vx = VX;
			sp[i].vy = VY;
			sp[i].vz = VZ;
		}else if(Colltype == 3){ //ionization 
			//printf("Ionization 1 ! \n");
            ///// scatter the created electron
			index = info[0].St_num + ID + (PNC + AddPt1 - 1) * Gsize; 
			sp[index].CellID = sp[i].CellID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
            dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,info_CX[4].mofM,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
            ///// assign velocities to the created ion
            index = info[1].St_num + ID + (PNC2 + AddPt1 - 1) * Gsize; 
            sp[index].CellID = sp[i].CellID + Gsize;
            sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
            sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
            n = (nvel-1)*curand_uniform(&LocalStates);
			dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[ID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			///// scatter the incident electron
			dev_anewvel(engy,vel2,&VX,&VY,&VZ,0,info_CX[4].mofM,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
            sp[i].vx = VX;
			sp[i].vy = VY;
			sp[i].vz = VZ;
			//AddPt1--;
		}
		i+=Gsize;
	}
    data[ID].PtNumInCell = PNC + AddPt1;
    data[ID+Gsize].PtNumInCell = PNC2 + AddPt1;
	states[ID]=LocalStates;
}
__device__ void Direct_Argon_ArIon(int Gsize, int ngy, int ID, int MCCn, float dt, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, GGA *BG, GFC *Fluid){
	int i,j,k,n,index;
	int PNC;
	float Prob;
	float R1;
	float VX,VY,VZ,VX_buf,VY_buf,VZ_buf;
	float dum,vel,engy;
	float SumSigma,SumEngyLoss;
	float vneutx,vneuty,vneutz;

    PNC = data[ID+Gsize].PtNumInCell;
	curandState LocalStates = states[ID];
	Prob = 1.0f - exp(-1*dt*sigv[2].val*BG[ID].BackDen1);

    // Calculate total Collision probability.
	i = info[1].St_num + ID;
	for(k=0;k<PNC;k++){
		R1 = curand_uniform(&LocalStates);
		if(R1>Prob){
			i+=Gsize;
			continue;
		} 
		n = (nvel-1)*curand_uniform(&LocalStates);
		dev_maxwellv(&vneutx,&vneuty,&vneutz,vsave[n],BG[ID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
		VX=sp[i].vx-vneutx;
	  	VY=sp[i].vy-vneutz;
	  	VZ=sp[i].vz-vneuty;
		dum=VX*VX+VY*VY+VZ*VZ;
		engy=info[1].Escale*dum;
	  	vel=sqrt(dum);
		R1 = curand_uniform(&LocalStates) * sigv[2].val / vel;
		SumSigma = Argon_CrossSection(5, engy, N_LOGX, idLOGX, CX);
		if(R1<=SumSigma){
			// 5. Ar + Ar^ > Ar + Ar^		Charge Exchange
			VX_buf = vneutx;
			VY_buf = vneuty;
			VZ_buf = vneutz;
		}else if(R1<=(SumSigma += Argon_CrossSection(6, engy, N_LOGX, idLOGX, CX))){
			// 6. AR + AR^ > AR + AR^		ELASTIC SCATTERING
			dev_newvel_IONSC(&VX,&VY,&VZ,vel,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			VX_buf = VX+vneutx;
			VY_buf = VY+vneuty;
			VZ_buf = VZ+vneutz;
		}else{
			VX_buf = sp[i].vx;
			VY_buf = sp[i].vy;
			VZ_buf = sp[i].vz;
		}
		sp[i].vx = VX_buf;
		sp[i].vy = VY_buf;
		sp[i].vz = VZ_buf;	
		i+=Gsize;
	}
	states[ID]=LocalStates;
}
__device__ void Collision_Check(int Gsize, int Csize, int ngy, int ID, int isp, float dt, int MCCn, float dtm, curandState *states, Species *info, GPG *data, GCP *sp, MCC_sigmav *sigv, GGA *BG, GFC *Fluid){
	int i,j,k,index,Randn;
	int TID,CID,PNC,PNMC,MPNC;
	int nx,ny,ngx;
	float Prob = 0.0f;
	float R1;
	TID = isp*Gsize + ID;
	nx = ID/ngy;
	ny = ID%ngy;
	ngx = Gsize/ngy;
	if(nx == ngx-1) nx--;
	if(ny == ngy-1) ny--;
	CID = ny + (ngy-1)*nx;
	curandState LocalStates = states[TID];
	PNC = data[TID].PtNumInCell;
	MPNC = data[TID].MaxPtNumInCell;
	PNMC = 0;
	// Calculate total Collision probability.
	if(isp == 0){ // Electron
		Prob = 2.0 - exp(-1*dtm*sigv[0].val*BG[ID].BackDen1) - exp(-1*dtm*sigv[1].val*Fluid[CID].ave_den);
		Randn = MCCn;
	}else{ // ion
		Prob = 1.0 - exp(-1*dt*sigv[2].val*BG[ID].BackDen1);
		Randn = 1;
	}
	i = info[isp].St_num + ID;
	for(k=0;k<PNC;k++){
		for(j=0;j<Randn;j++){
			R1 = curand_uniform(&LocalStates);
			if(R1<Prob) break;
		}
		if(R1 > Prob){ // no collision
			index = i - PNMC*Gsize;
		}else{ // collision
			PNMC++;
			index = info[isp].St_num + ID + (MPNC-PNMC)*Gsize;
		}
		sp[index].CellID = sp[i].CellID;
		sp[index].vx=sp[i].vx;
		sp[index].vy=sp[i].vy;
		sp[index].vz=sp[i].vz;
        sp[index].x=sp[i].x;
		sp[index].y=sp[i].y;
		i+=Gsize;
	}
	states[TID]=LocalStates;
	data[TID].PtNumMCCInCell=PNMC;
	data[TID].PtNumInCell-=PNMC;
}

__device__ void Argon_Ar_Collision(int Gsize, int TID, float dt, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, GGA *BG, GFC *Fluid){
	int i,k,n,index;
	int ID,PNMC,MPNC,PNC,AddPt;
	int Target,oldPNC;
	float prob;
	float cut1;
	float R1,R2;
	float VX,VY,VZ,VX_buf,VY_buf,VZ_buf;
	float dum,vel,vel2,engy,rengy;
	float SumSigma,SumEngyLoss;
	float vneutx,vneuty,vneutz;
	ID = Gsize + TID;
	PNC = data[ID].PtNumInCell;
	PNMC = data[ID].PtNumMCCInCell;
	MPNC = data[ID].MaxPtNumInCell;
	AddPt = 0;
	curandState LocalStates = states[ID];
	prob = 1.0f - exp(-1*dt*sigv[2].val*BG[TID].BackDen1);
	i = info[1].St_num + TID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
		// Calculate energy
		n = (nvel-1)*curand_uniform(&LocalStates);
		dev_maxwellv(&vneutx,&vneuty,&vneutz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
		VX = sp[i].vx - vneutx;
		VY = sp[i].vy - vneuty;
		VZ = sp[i].vz - vneutz;
		dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		engy = info[1].Escale * dum;
		R1 = curand_uniform(&LocalStates) * sigv[2].val / vel;
		SumSigma = Argon_CrossSection(5, engy, N_LOGX, idLOGX, CX);
		if(R1<=SumSigma){
			// 5. Ar + Ar^ > Ar + Ar^		Charge Exchange
			VX_buf = vneutx;
			VY_buf = vneuty;
			VZ_buf = vneutz;
		}else if(R1<=(SumSigma += Argon_CrossSection(6, engy, N_LOGX, idLOGX, CX))){
			// 6. AR + AR^ > AR + AR^		ELASTIC SCATTERING
			dev_newvel_IONSC(&VX,&VY,&VZ,vel,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			VX_buf = VX+vneutx;
			VY_buf = VY+vneuty;
			VZ_buf = VZ+vneutz;
		}else{
			VX_buf = sp[i].vx;
			VY_buf = sp[i].vy;
			VZ_buf = sp[i].vz;
		}
		index = info[1].St_num + TID + (PNC+AddPt) * Gsize; 
		sp[index].CellID = ID;
		sp[index].x = sp[i].x;
		sp[index].y = sp[i].y;
		sp[index].vx = VX_buf;
		sp[index].vy = VY_buf;
		sp[index].vz = VZ_buf;		
		i-=Gsize;
		AddPt++;
	}
	data[ID].PtNumInCell = PNC + AddPt;
	states[TID]=LocalStates;
}
__device__ void Argon_E_Collision(int Gsize, int ngy, int TID, int MCCn, float dtm, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, GGA *BG, GFC *Fluid){
	int i,k,n,index;
	int CID,PNMC,MPNC;
	int PNC,PNC_ion,AddPt,AddPt_ion;
	int nx,ny,ngx;
	int Target,oldPNC;
	int Colltype;
	float Tprob,Prob1,Prob2;
	float R1,R2;
	float VX,VY,VZ,VX_buf,VY_buf,VZ_buf;
	float dum,vel,vel2,engy,rengy;
	float SumSigma,SumEngyLoss;
	PNC = data[TID].PtNumInCell;
	PNC_ion = data[TID+Gsize].PtNumInCell;
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	AddPt = 0;
	AddPt_ion = 0;
	curandState LocalStates = states[TID];
	nx = TID/ngy;
	ny = TID%ngy;
	ngx = Gsize/ngy;
	if(nx == ngx-1) nx--;
	if(ny == ngy-1) ny--;
	CID = ny + (ngy-1)*nx;
	Prob1 = 1.0f - exp(-1*dtm*sigv[0].val*BG[TID].BackDen1);
	Prob2 = 1.0f - exp(-1*dtm*sigv[1].val*Fluid[CID].ave_den);
	Tprob = Prob1 + Prob2;
	i = info[0].St_num + TID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
		// Collision target setting
		R1 = Tprob * curand_uniform(&LocalStates);
		if(R1 <= Prob1)	Target = (int)0;
		else			Target = (int)1;
		// Calculate energy
		VX = sp[i].vx;
		VY = sp[i].vy;
		VZ = sp[i].vz;
		dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		VX/=vel; VY/=vel; VZ/=vel;
		engy = info[0].Escale * dum;
		Colltype = 2;
		switch(Target){
			case 0:{ // Target == Argon
				R2 = curand_uniform(&LocalStates) * sigv[0].val / vel;
				SumSigma = Argon_CrossSection(0, engy, N_LOGX, idLOGX, CX);
				if(R2<=SumSigma){ 
					// 0. e + Ar > e + Ar 	 	Elastic Scattering
				}else if(engy > info_CX[1].Th_e && R2<=(SumSigma += Argon_CrossSection(1, engy, N_LOGX, idLOGX, CX))){
					// 1. e + Ar > e + Ar* 		Excitation to Total Excited state
					engy-=info_CX[1].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
				}else if(engy > info_CX[2].Th_e && R2<=(SumSigma += Argon_CrossSection(2, engy, N_LOGX, idLOGX, CX))){
					// 2. e + Ar > e + Ar* 		Excitation to AR4SM
					engy-=info_CX[2].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
				}else if(engy > info_CX[3].Th_e && R2<=(SumSigma += Argon_CrossSection(3, engy, N_LOGX, idLOGX, CX))){
					// 3. e + Ar > e + e + Ar^	Direct ionization
					Colltype = 3;
					engy-=info_CX[3].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
				}else{
					Colltype = 1;
				}
				break;
			}
			case 1:{ // Target == Metastable argon
				R2 = curand_uniform(&LocalStates)*sigv[1].val / vel;
				SumSigma = Argon_CrossSection(4, engy, N_LOGX, idLOGX, CX);
				if(engy > info_CX[4].Th_e && R2<=SumSigma){
					// 4. e + Ar* > e + e + Ar^		step ionization
					Colltype = 3;
					engy-=info_CX[4].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
				}else{
					Colltype = 1;
				}
				break;
			}
		}
		//PNC,AddPt
		if(Colltype == 1){ // null Collision
			index = info[0].St_num + TID + (PNC+AddPt) * Gsize; 
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = sp[i].vx;
			sp[index].vy = sp[i].vy;
			sp[index].vz = sp[i].vz;
			AddPt++;
		}else if(Colltype == 2){ // Just energy loss
			dev_anewvel(engy,vel,&VX,&VY,&VZ,0,info_CX[4].mofM,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			index = info[0].St_num + TID + (PNC+AddPt) * Gsize; 
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			AddPt++;
		}else if(Colltype == 3){ //ionization 
			index = info[0].St_num + TID + (PNC+AddPt) * Gsize; 
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			dev_anewvel(rengy,vel2,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,info_CX[4].mofM,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			AddPt++;
			// New Electron create
			dev_anewvel(engy,vel2,&VX,&VY,&VZ,0,info_CX[4].mofM,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			index = info[0].St_num + TID + (PNC+AddPt) * Gsize; 
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			AddPt++;
			// ion create
			index = info[1].St_num + TID + (PNC_ion+AddPt_ion) * Gsize; 
			n = (nvel-1)*curand_uniform(&LocalStates);
			dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			sp[index].CellID = TID + Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			AddPt_ion++;
		}
		i-=Gsize;
	}
	data[TID].PtNumInCell = PNC + AddPt;
	data[TID+Gsize].PtNumInCell = PNC_ion + AddPt_ion;
	states[TID]=LocalStates;
}

__device__ float Argon_CrossSection(int R, float engy, int N_LOGX, float idLOGX, ArCollD *data){
	if(engy == 0) return 0.0;
	float lengy = log10(engy);
	float ee1, a1, a2;
	int ee2;
	lengy = lengy - data[0].xe;
	ee1 = idLOGX * lengy;
	ee2 = (int)ee1;
	a1 = ee1 - ee2;
	a2 = 1 - a1;
	switch (R) {
        case 0 : 
			if(lengy < data[0].xe){
				return data[0].cx_0;
			}else if(lengy > data[N_LOGX-1].xe){
				return data[N_LOGX-1].cx_0 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
			}
			return a2*data[ee2].cx_0+a1*data[ee2+1].cx_0;
			break;
        case 1 :
			if(lengy < data[0].xe){
				return data[0].cx_1;
			}else if(lengy > data[N_LOGX-1].xe){
				return data[N_LOGX-1].cx_1 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
			}
			return a2*data[ee2].cx_1+a1*data[ee2+1].cx_1;
			break;
        case 2 :
			if(lengy < data[0].xe){
				return data[0].cx_2;
			}else if(lengy > data[N_LOGX-1].xe){
				return data[N_LOGX-1].cx_2 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
			}
			return a2*data[ee2].cx_2+a1*data[ee2+1].cx_2;
			break;
        case 3 :
			if(lengy < data[0].xe){
				return data[0].cx_3;
			}else if(lengy > data[N_LOGX-1].xe){
				return data[N_LOGX-1].cx_3 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
			}
			return a2*data[ee2].cx_3+a1*data[ee2+1].cx_3;
			break;
        case 4 :
			if(lengy < data[0].xe){
				return data[0].cx_4;
			}else if(lengy > data[N_LOGX-1].xe){
				return data[N_LOGX-1].cx_4 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
			}
			return a2*data[ee2].cx_4+a1*data[ee2+1].cx_4;
			break;
        case 5 :
			if(lengy < data[0].xe){
				return data[0].cx_5;
			}else if(lengy > data[N_LOGX-1].xe){
				return data[N_LOGX-1].cx_5 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
			}
			return a2*data[ee2].cx_5+a1*data[ee2+1].cx_5;
			break;
        case 6 :
			if(lengy < data[0].xe){
				return data[0].cx_6;
			}else if(lengy > data[N_LOGX-1].xe){
				return data[N_LOGX-1].cx_6 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
			}
			return a2*data[ee2].cx_6+a1*data[ee2+1].cx_6;
			break;
        default :
            printf("\nError : Call about cross section data in ARMCC.\n\n");
            return 0.0;
    }
}