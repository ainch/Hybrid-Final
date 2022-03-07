#include "cuda_mccAr.cuh"
__device__ void MCC_Argon_RC(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
	int TnRct, float *MCCR, float *Stack, CollF *info_CX, Species *info, GPG *data, GCP *sp, GFG *FG){
	float PNC,M1;
	int oldPNC,index,index2,n;
	int xi = TID/ngy;
	int yi = TID%ngy;
	if(yi>ngy-2) return;
	if(xi>Gsize/ngy-2) return;
	curandState LocalStates = states[TID];
	M1 = 0;
	PNC = Stack[TID];
	PNC += info_CX[8].RR*FG[TID].n*FG[TID].n * info[0].MCCscale;
	while(PNC >= 1.0f){
		//printf("[%d]=[%d][%d]PNC = %g %g %g %g \n",TID,xi,yi,PNC,info_CX[8].RR,FG[CID].den,info[isp].MCCscale);
		oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
		index = info[0].St_num + TID + oldPNC*Gsize;
		sp[index].CellID = TID;
		sp[index].x = curand_uniform(&LocalStates);
		sp[index].y = curand_uniform(&LocalStates);
		n = (nvel-1) * curand_uniform(&LocalStates);
		dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],info[0].vti,curand_uniform(&LocalStates),curand_uniform(&LocalStates));	
		oldPNC = atomicAdd(&data[TID+Gsize].PtNumInCell,1);
		index2 = info[1].St_num + TID + oldPNC*Gsize;
		sp[index2].CellID = TID + Gsize;
		sp[index2].x = sp[index].x;
		sp[index2].y = sp[index].y;
		n = (nvel-1) * curand_uniform(&LocalStates);
		dev_maxwellv(&sp[index2].vx,&sp[index2].vy,&sp[index2].vz,vsave[n],info[1].vti,curand_uniform(&LocalStates),curand_uniform(&LocalStates));	
		PNC--;
		M1 += 0.25;
		break;
	}
	atomicAdd(&MCCR[TID*TnRct + 8] ,M1);
	atomicAdd(&MCCR[(TID+1)*TnRct + 8] ,M1);
	atomicAdd(&MCCR[(TID+ngy)*TnRct + 8] ,M1);
	atomicAdd(&MCCR[(TID+ngy+1)*TnRct + 8] ,M1);
	Stack[TID] = PNC;
    states[TID] = LocalStates;
}
__device__ void Direct_Argon_Electron(int Gsize, int ngy, int ID, int MCCn, float dtm, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, int TnRct, float*MCCR, GGA *BG, GFG *Fluid){
	int i,j,k,n,index;
	int PNC,PNC2;
	int Target;
	int Colltype;
	float Prob1,Prob2;
    int AddPt;
	float R1,R2;
	float VX,VY,VZ;
	float dum,vel,vel2,engy,rengy,massRatio;
	float SumSigma,SumEngyLoss;
	PNC = data[ID].PtNumInCell;
    PNC2 = data[ID+Gsize].PtNumInCell;
	curandState LocalStates = states[ID];
	Prob1 = 1.0f - exp(-1*dtm*sigv[0].val*BG[ID].BackDen1);
	Prob2 = Prob1 + 1.0f - exp(-1*dtm*sigv[1].val*Fluid[ID].n);
    // Calculate total Collision probability.
    AddPt = 0;
	massRatio = info_CX[4].mofM;
	i = info[0].St_num + ID;


	int x = PNC;
	int y = 0;
	
	if(true)
	{
		while(x != 0)
		{
			int mask0 = 0;
			int mask1 = 0;
			int mask;
			for(int j=0;j<30;j++)
			{
				Colltype = 0;
				for(int k=0;k<MCCn;k++){
					R1 = curand_uniform(&LocalStates);
					if(R1<Prob2){
						Colltype = 1;
						break;
					}
				}
				if(Colltype == 1)
				{				
					if(R1 <= Prob1)
					{
						//Target = (int)0;
						mask0 = mask0 | (1<<j);
					}	
					else
					{
						//Target = (int)1;
						mask1 = mask1 | (1<<j);
					}
				}
				--x;
				if(x == 0)break;
			}
			mask = mask0 | mask1;
			
			int st = 0;
			while(mask0 != 0)
			{
				Colltype = 1;
				int lb = __ffs(mask0);
				i = (y + st + lb - 1) * Gsize + info[0].St_num + ID;
				Target = (int)0;
				
				{
					// Calculate energy
					VX = sp[i].vx;
					VY = sp[i].vy;
					VZ = sp[i].vz;
					dum = VX*VX+VY*VY+VZ*VZ;
					vel = sqrt(dum);
					VX/=vel; VY/=vel; VZ/=vel;
					engy = info[0].Escale * dum;
					{
						R2 = curand_uniform(&LocalStates) * sigv[0].val / vel;
						// 0. e + Ar > e + Ar 			Elastic Scattering
						SumSigma = Argon_CrossSection(0, engy, N_LOGX, idLOGX, CX);
						if(R2<=SumSigma){
						// 1. e + Ar > e + Ar* 			Excitation to Total Excited state
							Colltype = 0;
							MCCR[ID*TnRct]++;
						}else if(engy > info_CX[1].Th_e && R2<=(SumSigma += Argon_CrossSection(1, engy, N_LOGX, idLOGX, CX))){
							Colltype = 0;
							engy-=info_CX[1].Th_e;
							vel=sqrt(fabs(engy)/info[0].Escale);
							MCCR[ID*TnRct+1]++;
						// 2. e + Ar > e + Ar* 			Excitation to AR4SM
						}else if(engy > info_CX[2].Th_e && R2<=(SumSigma += Argon_CrossSection(2, engy, N_LOGX, idLOGX, CX))){
							Colltype = 0;
							engy-=info_CX[2].Th_e;
							vel=sqrt(fabs(engy)/info[0].Escale);
							MCCR[ID*TnRct+2]++;
						// 3. e + Ar > e + e + Ar^		Direct ionization
						}else if(engy > info_CX[3].Th_e && R2<=(SumSigma += Argon_CrossSection(3, engy, N_LOGX, idLOGX, CX))){
							Colltype = 2;
							engy-=info_CX[3].Th_e;
							rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
							engy-=rengy;
							vel = sqrt(fabs(rengy)/info[0].Escale);
							vel2 = sqrt(fabs(engy)/info[0].Escale);
							MCCR[ID*TnRct+3]++;
						}
					}
					if(Colltype == 0){ // Just energy loss
						dev_anewvel(engy,vel,&VX,&VY,&VZ,0,massRatio,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
						sp[i].vx = VX;
						sp[i].vy = VY;
						sp[i].vz = VZ;
					}else if(Colltype == 2){ //ionization 
						//printf("Ionization 1 ! \n");
						///// scatter the created electron
						index = info[0].St_num + ID + (PNC + AddPt) * Gsize; 
						sp[index].CellID = ID;
						sp[index].x = sp[i].x;
						sp[index].y = sp[i].y;
						sp[index].vx = VX;
						sp[index].vy = VY;
						sp[index].vz = VZ;
						dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,massRatio,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
						///// assign velocities to the created ion
						index = info[1].St_num + ID + (PNC2 + AddPt) * Gsize; 
						sp[index].CellID = ID + Gsize;
						sp[index].x = sp[i].x;
						sp[index].y = sp[i].y;
						sp[index].vx = VX;
						sp[index].vy = VY;
						sp[index].vz = VZ;
						//printf("\n[%d][%d] ionization \n",ID,ID+Gsize);
						n = (nvel-1)*curand_uniform(&LocalStates);
						dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[ID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
						///// scatter the incident electron
						dev_anewvel(engy,vel2,&VX,&VY,&VZ,0,massRatio,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
						sp[i].vx = VX;
						sp[i].vy = VY;
						sp[i].vz = VZ;
						AddPt++;
					}
				}
	
				mask0 = mask0 >> lb;
				st += lb;
			}
			
			st = 0;
			while(mask1 != 0)
			{
				Colltype = 1;
				int lb = __ffs(mask);
				i = (y + st + lb - 1) * Gsize + info[0].St_num + ID;
				
				Target = (int)1;
				
				{
					// Calculate energy
					VX = sp[i].vx;
					VY = sp[i].vy;
					VZ = sp[i].vz;
					dum = VX*VX+VY*VY+VZ*VZ;
					vel = sqrt(dum);
					VX/=vel; VY/=vel; VZ/=vel;
					engy = info[0].Escale * dum;
					{
						R2 = curand_uniform(&LocalStates)*sigv[1].val / vel;
						// 4. e + Ar* > e + e + Ar^		step ionization
						SumSigma = Argon_CrossSection(4, engy, N_LOGX, idLOGX, CX);
						if(engy > info_CX[4].Th_e && R2<=SumSigma){
							Colltype = 2;
							engy-=info_CX[4].Th_e;
							rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
							engy-=rengy;
							vel = sqrt(fabs(rengy)/info[0].Escale);
							vel2 = sqrt(fabs(engy)/info[0].Escale);
							MCCR[ID*TnRct+4]++;
						}
					}
					if(Colltype == 0){ // Just energy loss
						dev_anewvel(engy,vel,&VX,&VY,&VZ,0,massRatio,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
						sp[i].vx = VX;
						sp[i].vy = VY;
						sp[i].vz = VZ;
					}else if(Colltype == 2){ //ionization 
						//printf("Ionization 1 ! \n");
						///// scatter the created electron
						index = info[0].St_num + ID + (PNC + AddPt) * Gsize; 
						sp[index].CellID = ID;
						sp[index].x = sp[i].x;
						sp[index].y = sp[i].y;
						sp[index].vx = VX;
						sp[index].vy = VY;
						sp[index].vz = VZ;
						dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,massRatio,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
						///// assign velocities to the created ion
						index = info[1].St_num + ID + (PNC2 + AddPt) * Gsize; 
						sp[index].CellID = ID + Gsize;
						sp[index].x = sp[i].x;
						sp[index].y = sp[i].y;
						sp[index].vx = VX;
						sp[index].vy = VY;
						sp[index].vz = VZ;
						//printf("\n[%d][%d] ionization \n",ID,ID+Gsize);
						n = (nvel-1)*curand_uniform(&LocalStates);
						dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[ID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
						///// scatter the incident electron
						dev_anewvel(engy,vel2,&VX,&VY,&VZ,0,massRatio,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
						sp[i].vx = VX;
						sp[i].vy = VY;
						sp[i].vz = VZ;
						AddPt++;
					}
				}
	
				mask1 = mask1 >> lb;
				st += lb;
			}

			y += 30;
		}
	}

	if(!true)
	{
		for(k=0;k<PNC;k++){
			Colltype = 0;
			for(j=0;j<MCCn;j++){
				R1 = curand_uniform(&LocalStates);
				if(R1<Prob2){
					Colltype = 1;
					break;
				}
			}
			
			if(Colltype == 1)
			{
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
				switch(Target){
					case 0:{
						R2 = curand_uniform(&LocalStates) * sigv[0].val / vel;
						// 0. e + Ar > e + Ar 			Elastic Scattering
						SumSigma = Argon_CrossSection(0, engy, N_LOGX, idLOGX, CX);
						if(R2<=SumSigma){
						// 1. e + Ar > e + Ar* 			Excitation to Total Excited state
							Colltype = 0;
							MCCR[ID*TnRct]++;
						}else if(engy > info_CX[1].Th_e && R2<=(SumSigma += Argon_CrossSection(1, engy, N_LOGX, idLOGX, CX))){
							Colltype = 0;
							engy-=info_CX[1].Th_e;
							vel=sqrt(fabs(engy)/info[0].Escale);
							MCCR[ID*TnRct+1]++;
						// 2. e + Ar > e + Ar* 			Excitation to AR4SM
						}else if(engy > info_CX[2].Th_e && R2<=(SumSigma += Argon_CrossSection(2, engy, N_LOGX, idLOGX, CX))){
							Colltype = 0;
							engy-=info_CX[2].Th_e;
							vel=sqrt(fabs(engy)/info[0].Escale);
							MCCR[ID*TnRct+2]++;
						// 3. e + Ar > e + e + Ar^		Direct ionization
						}else if(engy > info_CX[3].Th_e && R2<=(SumSigma += Argon_CrossSection(3, engy, N_LOGX, idLOGX, CX))){
							Colltype = 2;
							engy-=info_CX[3].Th_e;
							rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
							engy-=rengy;
							vel = sqrt(fabs(rengy)/info[0].Escale);
							vel2 = sqrt(fabs(engy)/info[0].Escale);
							MCCR[ID*TnRct+3]++;
						}
						break;
					}
					case 1:{
						R2 = curand_uniform(&LocalStates)*sigv[1].val / vel;
						// 4. e + Ar* > e + e + Ar^		step ionization
						SumSigma = Argon_CrossSection(4, engy, N_LOGX, idLOGX, CX);
						if(engy > info_CX[4].Th_e && R2<=SumSigma){
							Colltype = 2;
							engy-=info_CX[4].Th_e;
							rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
							engy-=rengy;
							vel = sqrt(fabs(rengy)/info[0].Escale);
							vel2 = sqrt(fabs(engy)/info[0].Escale);
							MCCR[ID*TnRct+4]++;
						}
						break;
					}
				}
				if(Colltype == 0){ // Just energy loss
					dev_anewvel(engy,vel,&VX,&VY,&VZ,0,massRatio,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
					sp[i].vx = VX;
					sp[i].vy = VY;
					sp[i].vz = VZ;
				}else if(Colltype == 2){ //ionization 
					//printf("Ionization 1 ! \n");
					///// scatter the created electron
					index = info[0].St_num + ID + (PNC + AddPt) * Gsize; 
					sp[index].CellID = ID;
					sp[index].x = sp[i].x;
					sp[index].y = sp[i].y;
					sp[index].vx = VX;
					sp[index].vy = VY;
					sp[index].vz = VZ;
					dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,massRatio,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
					///// assign velocities to the created ion
					index = info[1].St_num + ID + (PNC2 + AddPt) * Gsize; 
					sp[index].CellID = ID + Gsize;
					sp[index].x = sp[i].x;
					sp[index].y = sp[i].y;
					sp[index].vx = VX;
					sp[index].vy = VY;
					sp[index].vz = VZ;
					//printf("\n[%d][%d] ionization \n",ID,ID+Gsize);
					n = (nvel-1)*curand_uniform(&LocalStates);
					dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[ID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
					///// scatter the incident electron
					dev_anewvel(engy,vel2,&VX,&VY,&VZ,0,massRatio,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
					sp[i].vx = VX;
					sp[i].vy = VY;
					sp[i].vz = VZ;
					AddPt++;
				}
			}
			i+=Gsize;
		}
	}



	//if(AddPt != 0) printf("[%d] PNC[%d] addpt[%d]\n",ID,PNC,AddPt);
    data[ID].PtNumInCell = PNC + AddPt;
    data[ID+Gsize].PtNumInCell = PNC2+AddPt;
	states[ID]=LocalStates;
}
__device__ void Direct_Argon_ArIon(int Gsize, int ngy, int ID, int MCCn, float dt, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, int TnRct, float*MCCR, GGA *BG, GFG *Fluid){
	int i,k,n;
	int PNC,Null;
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

	int x = PNC;
	int y = 0;
	
	if(true)
	{
		while(x != 0)
		{
			int mask = 0;
			for(int j=0;j<30;j++)
			{
				R1 = curand_uniform(&LocalStates);
				if(R1 <= Prob)
				{
					mask = (1<<j) | mask;
				}
				--x;
				if(x == 0)break;
			}
			
			int st = 0;
			while(mask != 0)
			{
				int lb = __ffs(mask);
				i = (y + st + lb - 1) * Gsize + info[1].St_num + ID;
				
				{
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
						MCCR[ID*TnRct+5]++;
					}else if(R1<=(SumSigma += Argon_CrossSection(6, engy, N_LOGX, idLOGX, CX))){
						// 6. AR + AR^ > AR + AR^		ELASTIC SCATTERING
						dev_newvel_IONSC(&VX,&VY,&VZ,vel,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
						VX_buf = VX+vneutx;
						VY_buf = VY+vneuty;
						VZ_buf = VZ+vneutz;
						MCCR[ID*TnRct+6]++;
					}else{
						VX_buf = sp[i].vx;
						VY_buf = sp[i].vy;
						VZ_buf = sp[i].vz;
					}
					sp[i].vx = VX_buf;
					sp[i].vy = VY_buf;
					sp[i].vz = VZ_buf;
				}
	
				mask = mask >> lb;
				st += lb;
			}
			
			y += 30;
		}
	}

	if(!true)
	{
		for(k=0;k<PNC;k++)
		{
			R1 = curand_uniform(&LocalStates);
			if(R1<=Prob)
			{
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
					MCCR[ID*TnRct+5]++;
				}else if(R1<=(SumSigma += Argon_CrossSection(6, engy, N_LOGX, idLOGX, CX))){
					// 6. AR + AR^ > AR + AR^		ELASTIC SCATTERING
					dev_newvel_IONSC(&VX,&VY,&VZ,vel,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
					VX_buf = VX+vneutx;
					VY_buf = VY+vneuty;
					VZ_buf = VZ+vneutz;
					MCCR[ID*TnRct+6]++;
				}else{
					VX_buf = sp[i].vx;
					VY_buf = sp[i].vy;
					VZ_buf = sp[i].vz;
				}
				sp[i].vx = VX_buf;
				sp[i].vy = VY_buf;
				sp[i].vz = VZ_buf;
			}
			i+=Gsize;
		}
	}
	states[ID]=LocalStates;
}
__device__ void Ar_Collision_Check(int Gsize, int Csize, int ngy, int TID, float dt, int MCCn, float dtm, float dx, float dy,
                                        curandState *states, Species *info, GPG *data, GCP *sp, MCC_sigmav *sigv, GGA *BG, GFG *Fluid){
	int i,j,k,index,Randn;
	int ID,isp,PNMC,MPNC;
    int PNC,Flag;
	float Tprob,Prob1,Prob2;
	float R1;
	ID = TID%Gsize;
    isp = TID/Gsize;
	curandState LocalStates = states[TID];
	PNC = data[TID].PtNumInCell;
	MPNC = data[TID].MaxPtNumInCell;
	PNMC = 0;
	// Calculate total Collision probability.
    switch (isp){
    case 0: // Electron
		Prob1 = 1.0f - exp(-1*dtm*sigv[0].val*BG[ID].BackDen1);  // E + Ar
		Prob2 = Prob1 + 1.0f - exp(-1*dtm*sigv[1].val*Fluid[ID].n);  // E + Ar*
	    Tprob = Prob2; 
		Randn = MCCn;
        break;
	case 1: // Ar+
		Tprob = 1.0 - exp(-1*dt*sigv[2].val*BG[ID].BackDen1);
		Randn = 1;
		break;
    default:
        break;
    }
	i = info[isp].St_num + ID;
	for(k=0;k<PNC;k++){
        for(j=0;j<Randn;j++){
			R1 = curand_uniform(&LocalStates);
			if(R1<Tprob) break;
		}
		if(R1 >= Tprob){ // no collision
			index = i - PNMC*Gsize;
            Flag = sp[i].CellID;
		}else{ // collision
			PNMC++;
			index = info[isp].St_num + ID + (MPNC-PNMC)*Gsize;
            switch (isp){
            case 0:
                if(R1 <= Prob1)	        Flag = (int)0;
		        else			        Flag = (int)1;
                break;
            case 1:
				Flag = sp[i].CellID;
                break;
            default:
                break;
            }
		}
		sp[index].CellID = Flag;
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
__device__ void Ar_Collision_Check_v2(int Gsize, int Csize, int ngy, int TID, float dt, int MCCn, float dtm, float dx, float dy,
                                        curandState *states, Species *info, GPG *data, GCP *sp, MCC_sigmav *sigv, GGA *BG, GFG *Fluid){
	int i,j,k,index,Randn;
	int PNMC,MPNC;
    int PNC,Flag;
	float Prob1,Prob2;
	float R1;
	curandState LocalStates = states[TID];
	PNC = data[TID].PtNumInCell;
	MPNC = data[TID].MaxPtNumInCell;
	PNMC = 0;
	// Calculate total Collision probability.
	Prob1 = 1.0f - exp(-1*dtm*sigv[0].val*BG[TID].BackDen1);  // E + Ar
	Prob2 = Prob1 + 1.0f - exp(-1*dtm*sigv[1].val*Fluid[TID].n);  // E + Ar*
	Randn = MCCn;
	i = TID;
	for(k=0;k<PNC;k++){
        for(j=0;j<Randn;j++){
			R1 = curand_uniform(&LocalStates);
			if(R1<Prob2) break;
		}
		if(R1 >= Prob2){ // no collision
			index = i - PNMC*Gsize;
            Flag = sp[i].CellID;
		}else{ // collision
			PNMC++;
			index = TID + (MPNC-PNMC)*Gsize;
            if(R1 <= Prob1)	        Flag = (int)0;
		    else			        Flag = (int)1;            
		}
		sp[index].CellID = Flag;
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
__device__ void Ar_Electron_v2(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, int TnRct,float *MCCR, GGA *BG){
	int i,k,n,index;
	int PNMC,MPNC,Flag;
	int Colltype;
	int PNC,PNC1,Addpt,Addpt1;
	float mofm,R1;
	float VX,VY,VZ;
	float rengy,engy,dum,vel,vel2;		
	float SumSigma,SumEngyLoss;

    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	// Calculate total Collision probability
	i = TID + (MPNC-1)*Gsize;
	Addpt = 0; Addpt1 = 0;
	PNC = data[TID].PtNumInCell;
	PNC1 = data[Gsize + TID].PtNumInCell;
	//printf("PNMC = %d\n",PNMC);
	for(k=0;k<PNMC;k++){
        // Calculate energy
		VX = sp[i].vx;
		VY = sp[i].vy;
		VZ = sp[i].vz;
        Flag = sp[i].CellID;
        dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		VX/=vel; VY/=vel; VZ/=vel;
		engy = info[0].Escale * dum;
        Colltype = 0;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Energy loss
		// 2 : ionization 1 -> 3 charged
		R1 = curand_uniform(&LocalStates);
        switch(Flag){
			case 0:{ // E + Ar
				mofm = info_CX[0].mofM;
				R1 *= sigv[0].val / vel;
				if(engy > info_CX[0].Th_e && R1<=(SumSigma=Argon_CrossSection(0, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1;
					MCCR[TID*TnRct]++;
				}else if(engy > info_CX[1].Th_e && R1<=(SumSigma += Argon_CrossSection(1, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1;
					engy-=info_CX[1].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+1]++;
				}else if(engy > info_CX[2].Th_e && R1<=(SumSigma += Argon_CrossSection(2, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1;
					engy-=info_CX[2].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+2]++;
				}else if(engy > info_CX[3].Th_e && R1<=(SumSigma += Argon_CrossSection(3, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2;
					engy-=info_CX[3].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+3]++;
				}				
				break;
			}
			case 1:{ // E + Ar*
				mofm = info_CX[4].mofM;
				R1 *= sigv[1].val / vel;
				if(engy > info_CX[4].Th_e && R1<=(SumSigma=Argon_CrossSection(4, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2;
					engy-=info_CX[4].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+4]++;
				}
				break;
			}
			default:{
            	break;
        	}
		} 
        switch (Colltype){
        case 0: // 0 : Null collision
			index = TID + (PNC+Addpt)*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = sp[i].vx;
			sp[index].vy = sp[i].vy;
			sp[index].vz = sp[i].vz;
			Addpt++;
            break;
        case 1: // 1 : Energy loss
            dev_anewvel(engy,vel,&VX,&VY,&VZ,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			index = TID + (PNC+Addpt)*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			Addpt++;
            break;
        case 2: // 2 : ionization 1 -> 3 charged
			// second charged create
			index = TID + (PNC+Addpt)*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			Addpt++;
			// Third charged create
			index = info[1].St_num + TID + (PNC1+Addpt1)*Gsize;
			sp[index].CellID = TID+Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			n = (nvel-1)*curand_uniform(&LocalStates);
			dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			Addpt1++;
			// energy loss electron 
			dev_anewvel(engy,vel2,&VX,&VY,&VZ,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			index = TID + (PNC+Addpt)*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			Addpt++;
            break;
        default:
            break;
        }
		i-=Gsize;
	}
	data[TID].PtNumInCell = PNC + Addpt;
	data[Gsize + TID].PtNumInCell = PNC1 + Addpt1;
	states[TID]=LocalStates;
}
__device__ void Ar_Electron(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, int TnRct,float *MCCR, GGA *BG){
	int i,k,n,index;
	int PNMC,MPNC,Flag;
	int oldPNC;
	int Colltype;
	float mofm,R1;
	float VX,VY,VZ;
	float rengy,engy,dum,vel,vel2;		
	float SumSigma,SumEngyLoss;

    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	// Calculate total Collision probability
	i = TID + (MPNC-1)*Gsize;
	//printf("PNMC = %d\n",PNMC);
	for(k=0;k<PNMC;k++){
        // Calculate energy
		VX = sp[i].vx;
		VY = sp[i].vy;
		VZ = sp[i].vz;
        Flag = sp[i].CellID;
        dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		VX/=vel; VY/=vel; VZ/=vel;
		engy = info[0].Escale * dum;
        Colltype = 0;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Energy loss
		// 2 : ionization 1 -> 3 charged
		R1 = curand_uniform(&LocalStates);
        switch(Flag){
			case 0:{ // E + Ar
				mofm = info_CX[0].mofM;
				R1 *= sigv[0].val / vel;
				if(engy > info_CX[0].Th_e && R1<=(SumSigma=Argon_CrossSection(0, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1;
					MCCR[TID*TnRct]++;
				}else if(engy > info_CX[1].Th_e && R1<=(SumSigma += Argon_CrossSection(1, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1;
					engy-=info_CX[1].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+1]++;
				}else if(engy > info_CX[2].Th_e && R1<=(SumSigma += Argon_CrossSection(2, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1;
					engy-=info_CX[2].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+2]++;
				}else if(engy > info_CX[3].Th_e && R1<=(SumSigma += Argon_CrossSection(3, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2;
					engy-=info_CX[3].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+3]++;
				}				
				break;
			}
			case 1:{ // E + Ar*
				mofm = info_CX[4].mofM;
				R1 *= sigv[1].val / vel;
				if(engy > info_CX[4].Th_e && R1<=(SumSigma=Argon_CrossSection(4, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2;
					engy-=info_CX[4].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+4]++;
				}
				break;
			}
			default:{
            	break;
        	}
		} 
        switch (Colltype){
        case 0: // 0 : Null collision
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = TID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = sp[i].vx;
			sp[index].vy = sp[i].vy;
			sp[index].vz = sp[i].vz;
            break;
        case 1: // 1 : Energy loss
            dev_anewvel(engy,vel,&VX,&VY,&VZ,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = TID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
            break;
        case 2: // 2 : ionization 1 -> 3 charged
			// second charged create
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = TID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			// Third charged create
			oldPNC = atomicAdd(&data[TID+Gsize].PtNumInCell,1);
			index = info[1].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID+Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			n = (nvel-1)*curand_uniform(&LocalStates);
			dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			// energy loss electron 
			dev_anewvel(engy,vel2,&VX,&VY,&VZ,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = TID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
            break;
        default:
            break;
        }
		i-=Gsize;
	}
	states[TID]=LocalStates;
}
__device__ void Ar_Ar_ion(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArCollD *CX, int TnRct,float *MCCR, GGA *BG){
	int i,k,n,index;
	int ID,PNMC,MPNC;
	int oldPNC;
	int Colltype;
	float R1;
	float VX,VY,VZ;
	float engy,dum,vel;	
	float SumSigma,SumEngyLoss;
	float vneutx,vneuty,vneutz;

	ID = TID%Gsize;
    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;

    // Calculate total Collision probability
	i = info[1].St_num + ID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
        // Calculate energy
		n = (nvel-1)*curand_uniform(&LocalStates);
		dev_maxwellv(&vneutx,&vneuty,&vneutz,vsave[n],BG[ID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
        VX = sp[i].vx - vneutx;
		VY = sp[i].vy - vneuty;
		VZ = sp[i].vz - vneutz;
		dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		engy = info[1].Escale * dum;
        Colltype = 0;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Scattering
        // 2 : Charge exchange 
        R1 = curand_uniform(&LocalStates)*sigv[2].val / vel;
		if(engy > info_CX[5].Th_e &&R1<=(SumSigma=Argon_CrossSection(5, engy, N_LOGX, idLOGX, CX))){
			Colltype = 2; 
			MCCR[ID*TnRct+5]++;
		}else if(engy > info_CX[6].Th_e &&R1<=(SumSigma=Argon_CrossSection(6, engy, N_LOGX, idLOGX, CX))){
			Colltype = 1; 
			MCCR[ID*TnRct+6]++;
		}
		switch (Colltype){
        case 0: // 0 : Null collision
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[1].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = sp[i].vx;
			sp[index].vy = sp[i].vy;
			sp[index].vz = sp[i].vz;
            break;
        case 1: // 1 : Scattering
			dev_newvel_IONSC(&VX,&VY,&VZ,vel,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[1].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX+vneutx;
			sp[index].vy = VY+vneuty;
			sp[index].vz = VZ+vneutz;
            break;
        case 2: // 2 : Charge exchange o2+
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[1].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = vneutx;
			sp[index].vy = vneuty;
			sp[index].vz = vneutz;
			break;
        default:
            break;
        }
		i-=Gsize;
	}
	states[TID]=LocalStates;
}
__device__ float Argon_CrossSection(int R, float engy, int N_LOGX, float idLOGX, ArCollD *data){
	if(engy == 0) return 0.0;
	float lengy = log10(engy);
	float ee1, a1, a2;
	int ee2;
	ee1 = idLOGX * (lengy - data[0].xe);
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