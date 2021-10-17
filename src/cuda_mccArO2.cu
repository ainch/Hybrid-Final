#include "cuda_mccArO2.cuh"
__device__ void ArO2_Electron(int Gsize, int ngy, int ID, int MCCn, float dtm, float dx, float dy, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, GGA *BG, GFC *Fluid){
/*
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
    */
}
__device__ void ArO2_ArIon(int Gsize, int ngy, int ID, int MCCn, float dt, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, GGA *BG, GFC *Fluid){
	/*
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
    */
}
__device__ float ArO2_CrossSection(int R, float engy, int N_LOGX, float idLOGX, ArO2CollD *data){
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
        case 7 :
            if(lengy < data[0].xe){
			    return data[0].cx_7;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_7 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_7+a1*data[ee2+1].cx_7;
            break;
        case 8 :
            if(lengy < data[0].xe){
			    return data[0].cx_8;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_8 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_8+a1*data[ee2+1].cx_8;
            break;
        case 9 :
            if(lengy < data[0].xe){
			    return data[0].cx_9;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_9 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_9+a1*data[ee2+1].cx_9;
            break;
        case 10 :
            if(lengy < data[0].xe){
			    return data[0].cx_10;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_10 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_10+a1*data[ee2+1].cx_10;
            break;
        case 11 :
            if(lengy < data[0].xe){
			    return data[0].cx_11;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_11 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_11+a1*data[ee2+1].cx_11;
            break;
        case 12 :
            if(lengy < data[0].xe){
			    return data[0].cx_12;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_12 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_12+a1*data[ee2+1].cx_12;
            break;
        case 13 :
            if(lengy < data[0].xe){
			    return data[0].cx_13;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_13 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_13+a1*data[ee2+1].cx_13;
            break;
        case 14 :
            if(lengy < data[0].xe){
			    return data[0].cx_14;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_14 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_14+a1*data[ee2+1].cx_14;
            break;
        case 15 :
            if(lengy < data[0].xe){
			    return data[0].cx_15;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_15 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_15+a1*data[ee2+1].cx_15;
            break;
        case 16 :
            if(lengy < data[0].xe){
			    return data[0].cx_16;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_16 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_16+a1*data[ee2+1].cx_16;
            break;
        case 17 :
            if(lengy < data[0].xe){
			    return data[0].cx_17;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_17 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_17+a1*data[ee2+1].cx_17;
            break;
        case 18 :
            if(lengy < data[0].xe){
			    return data[0].cx_18;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_18 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_18+a1*data[ee2+1].cx_18;
            break;
        case 19 :
            if(lengy < data[0].xe){
			    return data[0].cx_19;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_19 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_19+a1*data[ee2+1].cx_19;
            break;
        case 20 :
            if(lengy < data[0].xe){
			    return data[0].cx_20;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_20 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_20+a1*data[ee2+1].cx_20;
            break;
        case 21 :
            if(lengy < data[0].xe){
			    return data[0].cx_21;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_21 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_21+a1*data[ee2+1].cx_21;
            break;
        case 22 :
            if(lengy < data[0].xe){
			    return data[0].cx_22;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_22 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_22+a1*data[ee2+1].cx_22;
            break;
        case 23 :
            if(lengy < data[0].xe){
			    return data[0].cx_23;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_23 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_23+a1*data[ee2+1].cx_23;
            break;
        case 24 :
            if(lengy < data[0].xe){
			    return data[0].cx_24;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_24 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_24+a1*data[ee2+1].cx_24;
            break;
        case 25 :
            if(lengy < data[0].xe){
			    return data[0].cx_25;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_25 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_25+a1*data[ee2+1].cx_25;
            break;
        case 26 :
            if(lengy < data[0].xe){
			    return data[0].cx_26;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_26 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_26+a1*data[ee2+1].cx_26;
            break;
        case 27 :
            if(lengy < data[0].xe){
			    return data[0].cx_27;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_27 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_27+a1*data[ee2+1].cx_27;
            break;
        case 28 :
            if(lengy < data[0].xe){
			    return data[0].cx_28;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_28 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_28+a1*data[ee2+1].cx_28;
            break;
        case 29 :
            if(lengy < data[0].xe){
			    return data[0].cx_29;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_29 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_29+a1*data[ee2+1].cx_29;
            break;
        case 30 :
            if(lengy < data[0].xe){
			    return data[0].cx_30;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_30 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_30+a1*data[ee2+1].cx_30;
            break;
        case 31 :
            if(lengy < data[0].xe){
			    return data[0].cx_31;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_31 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_31+a1*data[ee2+1].cx_31;
            break;
        case 32 :
            if(lengy < data[0].xe){
			    return data[0].cx_32;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_32 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_32+a1*data[ee2+1].cx_32;
            break;
        case 33 :
            if(lengy < data[0].xe){
			    return data[0].cx_33;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_33 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_33+a1*data[ee2+1].cx_33;
            break;
        case 34 :
            if(lengy < data[0].xe){
			    return data[0].cx_34;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_34 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_34+a1*data[ee2+1].cx_34;
            break;
        case 35 :
            if(lengy < data[0].xe){
			    return data[0].cx_35;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_35 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_35+a1*data[ee2+1].cx_35;
            break;
        case 36 :
            if(lengy < data[0].xe){
			    return data[0].cx_36;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_36 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_36+a1*data[ee2+1].cx_36;
            break;
        case 37 :
            if(lengy < data[0].xe){
			    return data[0].cx_37;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_37 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_37+a1*data[ee2+1].cx_37;
            break;
        case 38 :
            if(lengy < data[0].xe){
			    return data[0].cx_38;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_38 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_38+a1*data[ee2+1].cx_38;
            break;
        case 39 :
            if(lengy < data[0].xe){
			    return data[0].cx_39;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_39 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_39+a1*data[ee2+1].cx_39;
            break;
        case 40 :
            if(lengy < data[0].xe){
			    return data[0].cx_40;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_40 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_40+a1*data[ee2+1].cx_40;
            break;
        case 41 :
            if(lengy < data[0].xe){
			    return data[0].cx_41;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_41 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_41+a1*data[ee2+1].cx_41;
            break;
        case 42 :
            if(lengy < data[0].xe){
			    return data[0].cx_42;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_42 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_42+a1*data[ee2+1].cx_42;
            break;
        case 43 :
            if(lengy < data[0].xe){
			    return data[0].cx_43;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_43 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_43+a1*data[ee2+1].cx_43;
            break;
        case 44 :
            if(lengy < data[0].xe){
			    return data[0].cx_44;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_44 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_44+a1*data[ee2+1].cx_44;
            break;
        case 45 :
            if(lengy < data[0].xe){
			    return data[0].cx_45;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_45 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_45+a1*data[ee2+1].cx_45;
            break;
        case 46 :
            if(lengy < data[0].xe){
			    return data[0].cx_46;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_46 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_46+a1*data[ee2+1].cx_46;
            break;
        case 47 :
            if(lengy < data[0].xe){
			    return data[0].cx_47;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_47 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_47+a1*data[ee2+1].cx_47;
            break;
        case 48 :
            if(lengy < data[0].xe){
			    return data[0].cx_48;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_48 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_48+a1*data[ee2+1].cx_48;
            break;
        case 49 :
            if(lengy < data[0].xe){
			    return data[0].cx_49;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_49 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_49+a1*data[ee2+1].cx_49;
            break;
        case 50 :
            if(lengy < data[0].xe){
			    return data[0].cx_50;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_50 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_50+a1*data[ee2+1].cx_50;
            break;
        case 51 :
            if(lengy < data[0].xe){
			    return data[0].cx_51;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_51 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_51+a1*data[ee2+1].cx_51;
            break;
        case 52 :
            if(lengy < data[0].xe){
			    return data[0].cx_52;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_52 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_52+a1*data[ee2+1].cx_52;
            break;
        case 53 :
            if(lengy < data[0].xe){
			    return data[0].cx_53;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_53 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_53+a1*data[ee2+1].cx_53;
            break;
        case 54 :
            if(lengy < data[0].xe){
			    return data[0].cx_54;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_54 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_54+a1*data[ee2+1].cx_54;
            break;
        case 55 :
            if(lengy < data[0].xe){
			    return data[0].cx_55;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_55 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_55+a1*data[ee2+1].cx_55;
            break;
        case 56 :
            if(lengy < data[0].xe){
			    return data[0].cx_56;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_56 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_56+a1*data[ee2+1].cx_56;
            break;
        case 57 :
            if(lengy < data[0].xe){
			    return data[0].cx_57;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_57 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_57+a1*data[ee2+1].cx_57;
            break;
        default :
            printf("\nError : Call about cross section data in ArO2MCC.\n\n");
            return 0.0;
    }
}