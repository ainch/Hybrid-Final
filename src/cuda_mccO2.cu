#include "cuda_mccO2.cuh"
__device__ void O2_Electron(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, O2CollD *CX, int TnRct, float*MCCR,GGA *BG){
	int i,j,k,n,index,index2,index3;
	int PNMC,MPNC,Null,Flag;
	int Target,oldPNC;
	int Colltype;
	float mofm,R1,R2;
	float VX,VY,VZ;
	float engy,dum,vel;
	int Iz_isp1,Iz_isp2;		
	float SumSigma,SumEngyLoss;

    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	Null = 0;
    // Calculate total Collision probability
	i = info[0].St_num + TID + (MPNC-1)*Gsize;
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
        Colltype = 1;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Energy loss
        // 2 : Attachment using maxwellv
		// 3 : ionization 1 -> 3 charged
		// 4 : dissociative  recombination just delete
		// 5 : Detachment 
        switch(Flag){
			case 0:{
				mofm = info_CX[0].mofM;
                R1 = curand_uniform(&LocalStates) * sigv[0].val / vel;
				if(engy > info_CX[0].Th_e &&R1<=(SumSigma=O2_CrossSection(0, engy, N_LOGX, idLOGX, CX))){
                    // R0 Elastic
					MCCR[TID*TnRct]++;
				}else if(engy > info_CX[1].Th_e && R1<=(SumSigma += O2_CrossSection(1, engy, N_LOGX, idLOGX, CX))){
                    //"1.e+O2>e+O2*");
					engy-=info_CX[1].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+1]++;
				}else if(engy > info_CX[2].Th_e && R1<=(SumSigma += O2_CrossSection(2, engy, N_LOGX, idLOGX, CX))){
                    //"2.e+O2>e+O2*");
					engy-=info_CX[2].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+2]++;
				}else if(engy > info_CX[3].Th_e && R1<=(SumSigma += O2_CrossSection(3, engy, N_LOGX, idLOGX, CX))){
                    //"3.e+O2>e+O2A");
					engy-=info_CX[3].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+3]++;
				}else if(engy > info_CX[4].Th_e && R1<=(SumSigma += O2_CrossSection(4, engy, N_LOGX, idLOGX, CX))){
					//"4.e+O2>e+O2B");
                    engy-=info_CX[4].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+4]++;
				}else if(engy > info_CX[5].Th_e && R1<=(SumSigma += O2_CrossSection(5, engy, N_LOGX, idLOGX, CX))){
                    //"5.e+O2>e+O2*");
					engy-=info_CX[5].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+5]++;
				}else if(engy > info_CX[6].Th_e && R1<=(SumSigma += O2_CrossSection(6, engy, N_LOGX, idLOGX, CX))){
                    //"6.e+O2>OP+O-"
                    Colltype = 2;
					MCCR[TID*TnRct+6]++;
                }else if(engy > info_CX[7].Th_e && R1<=(SumSigma += O2_CrossSection(7, engy, N_LOGX, idLOGX, CX))){
                    //"7.e+O2>e+2OP");
					engy-=info_CX[7].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+7]++;
                }else if(engy > info_CX[8].Th_e && R1<=(SumSigma += O2_CrossSection(8, engy, N_LOGX, idLOGX, CX))){
                    //"8.e+O2>e+OP+OD");
					engy-=info_CX[8].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+8]++;
                }else if(engy > info_CX[9].Th_e && R1<=(SumSigma += O2_CrossSection(9, engy, N_LOGX, idLOGX, CX))){
                    //"9.e+O2>e+2OD");
					engy-=info_CX[9].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+9]++;
                }else if(engy > info_CX[10].Th_e && R1<=(SumSigma += O2_CrossSection(10, engy, N_LOGX, idLOGX, CX))){
                    //"10.e+O2>2e+O2^");
                    Colltype = 3;
					engy-=info_CX[10].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 1;
					MCCR[TID*TnRct+10]++;
                }else if(engy > info_CX[11].Th_e && R1<=(SumSigma += O2_CrossSection(11, engy, N_LOGX, idLOGX, CX))){
                    //"11.e+O2>e+OP+O*");
					engy-=info_CX[11].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+11]++;
                }else if(engy > info_CX[12].Th_e && R1<=(SumSigma += O2_CrossSection(12, engy, N_LOGX, idLOGX, CX))){
                    //"12.e+O2>e+O^+O-");
                    Colltype = 3;
					engy-=info_CX[12].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 2;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+12]++;
                }else if(engy > info_CX[13].Th_e && R1<=(SumSigma += O2_CrossSection(13, engy, N_LOGX, idLOGX, CX))){
                    //"13.e+O2>2e+O^+OP");  
                    Colltype = 3;
					engy-=info_CX[13].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+13]++;
				}else{
					Colltype = 0;
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[14].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[1].val / vel;
				if(engy > info_CX[14].Th_e &&R1<=(SumSigma=O2_CrossSection(14, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"14.e+O2A>2e+O2+");
					engy-=info_CX[14].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 1;
					MCCR[TID*TnRct+14]++;
				}else if(engy > info_CX[15].Th_e && R1<=(SumSigma += O2_CrossSection(15, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; //"15.e+O2A>OP+O-");
					MCCR[TID*TnRct+15]++;
				}else if(engy > info_CX[16].Th_e && R1<=(SumSigma += O2_CrossSection(16, engy, N_LOGX, idLOGX, CX))){
					engy+=0.977f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+16]++;
				}else if(engy > info_CX[17].Th_e && R1<=(SumSigma += O2_CrossSection(17, engy, N_LOGX, idLOGX, CX))){
					engy+=0.977f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+17]++;
				}else if(engy > info_CX[18].Th_e && R1<=(SumSigma += O2_CrossSection(18, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[18].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+18]++;
				}else if(engy > info_CX[19].Th_e && R1<=(SumSigma += O2_CrossSection(19, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[19].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+19]++;
				}else if(engy > info_CX[20].Th_e && R1<=(SumSigma += O2_CrossSection(20, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[20].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+20]++;
				}else if(engy > info_CX[21].Th_e && R1<=(SumSigma += O2_CrossSection(21, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"21.e+O2A>2e+O^+OP");
					engy-=info_CX[21].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+21]++;
				}else{
					Colltype = 0;
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[22].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[2].val / vel;
                if(engy > info_CX[22].Th_e &&R1<=(SumSigma=O2_CrossSection(22, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"22.e+O2B>2e+O2^");
					engy-=info_CX[22].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 1;
					MCCR[TID*TnRct+22]++;
				}else if(engy > info_CX[23].Th_e && R1<=(SumSigma += O2_CrossSection(23, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; //"23.e+O2B>OP+O-");
					MCCR[TID*TnRct+23]++;
				}else if(engy > info_CX[24].Th_e && R1<=(SumSigma += O2_CrossSection(24, engy, N_LOGX, idLOGX, CX))){
					engy+=1.627f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+24]++;
				}else if(engy > info_CX[25].Th_e && R1<=(SumSigma += O2_CrossSection(25, engy, N_LOGX, idLOGX, CX))){
					engy+=1.627f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+25]++;
				}else if(engy > info_CX[26].Th_e && R1<=(SumSigma += O2_CrossSection(26, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[26].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+26]++;
				}else if(engy > info_CX[27].Th_e && R1<=(SumSigma += O2_CrossSection(27, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[27].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+27]++;
				}else if(engy > info_CX[28].Th_e && R1<=(SumSigma += O2_CrossSection(28, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[28].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+28]++;
				}else if(engy > info_CX[29].Th_e && R1<=(SumSigma += O2_CrossSection(29, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"29.e+O2B>2e+O^+OP");
					engy-=info_CX[29].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+29]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[30].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[3].val / vel;
                if(engy > info_CX[30].Th_e && R1<=(SumSigma=O2_CrossSection(30, engy, N_LOGX, idLOGX, CX))){
					Colltype = 5; //"30.e+O->2e+OP");
					engy-=info_CX[30].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+30]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 4:{
				mofm = info_CX[31].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[4].val / vel;
                if(engy > info_CX[31].Th_e &&R1<=(SumSigma=O2_CrossSection(31, engy, N_LOGX, idLOGX, CX))){
					Colltype = 4; //"31.e+O2^>OP+OD");
					MCCR[TID*TnRct+31]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 5:{
				mofm = info_CX[32].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[5].val / vel;
                if(engy > info_CX[32].Th_e &&R1<=(SumSigma=O2_CrossSection(32, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[32].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+32]++;
				}else if(engy > info_CX[33].Th_e && R1<=(SumSigma += O2_CrossSection(33, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[33].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+33]++;
				}else if(engy > info_CX[34].Th_e && R1<=(SumSigma += O2_CrossSection(34, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[34].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+34]++;
				}else if(engy > info_CX[35].Th_e && R1<=(SumSigma += O2_CrossSection(35, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[35].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+35]++;
				}else if(engy > info_CX[36].Th_e && R1<=(SumSigma += O2_CrossSection(36, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[36].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+36]++;
				}else if(engy > info_CX[37].Th_e && R1<=(SumSigma += O2_CrossSection(37, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"37.e+OP>2e+O^");
					engy-=info_CX[37].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+37]++;
				}else if(engy > info_CX[38].Th_e && R1<=(SumSigma += O2_CrossSection(38, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[38].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+38]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 6:{
				mofm = info_CX[39].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[6].val / vel;
                if(engy > info_CX[39].Th_e &&R1<=(SumSigma=O2_CrossSection(39, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"39.e+OD>2e+O^");
					engy-=info_CX[39].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+39]++;
				}else if(engy > info_CX[40].Th_e && R1<=(SumSigma += O2_CrossSection(40, engy, N_LOGX, idLOGX, CX))){
					engy+=1.96f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+40]++;
				}else{
					Colltype = 0;
					Null++;
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
			index = info[0].St_num + TID + oldPNC*Gsize;
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
			index = info[0].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
            break;
        case 2: // 2 : Attachment using maxwellv
			oldPNC = atomicAdd(&data[TID+3*Gsize].PtNumInCell,1);
			index = info[3].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID+3*Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			n = (nvel-1)*curand_uniform(&LocalStates);
			dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
            break;
        case 3: // 3 : ionization 1 -> 3 charged
			// second charged create
			oldPNC = atomicAdd(&data[TID+Iz_isp1*Gsize].PtNumInCell,1);
			index = info[Iz_isp1].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID+Iz_isp1*Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			if(Iz_isp1 == 0){
				sp[index].vx = VX;
				sp[index].vy = VY;
				sp[index].vz = VZ;
			}else{
				n = (nvel-1)*curand_uniform(&LocalStates);
				dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			}
			// Third charged create
			oldPNC = atomicAdd(&data[TID+Iz_isp2*Gsize].PtNumInCell,1);
			index = info[Iz_isp2].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID+Iz_isp2*Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			if(Iz_isp2 == 0){
				sp[index].vx = VX;
				sp[index].vy = VY;
				sp[index].vz = VZ;
			}else{
				n = (nvel-1)*curand_uniform(&LocalStates);
				dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			}
			// energy loss electron 
			dev_anewvel(engy,vel,&VX,&VY,&VZ,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[0].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
            break;
        case 4: // 4 : dissociative  recombination just delete
            break;
		case 5: // 5 : // Detachment 
			//"30.e+O->2e+OP");
			// new electron
			//printf("Ecollision case 5\n");
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[0].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			// energy loss
			dev_anewvel(engy,vel,&VX,&VY,&VZ,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[0].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			// delete O-
			oldPNC = atomicAdd(&data[TID+3*Gsize].PtNumInCell,0);
			if(oldPNC>1){
				R2 = curand_uniform(&LocalStates);
				index = info[3].St_num + TID + oldPNC*Gsize;
				index2 = (int)((float)oldPNC * R2);
				index3 = info[3].St_num + TID + index2*Gsize;
				sp[index3].CellID = sp[index].CellID;
				sp[index3].x = sp[index].x;
				sp[index3].y = sp[index].y;
				sp[index3].vx = sp[index].vx;
				sp[index3].vy = sp[index].vy;
				sp[index3].vz = sp[index].vz;
				atomicAdd(&data[TID+3*Gsize].PtNumInCell,-1);
				//printf("2[%d][%d]: %g,%g,%g,%g,%g,\n",TID,sp[index].x,sp[index].y,sp[index].vx,sp[index].vy,sp[i].vz);
			}else if(oldPNC == 1){
				atomicAdd(&data[TID+3*Gsize].PtNumInCell,-1);
			}
            break;
        default:
            break;
        }
		i-=Gsize;
	}
	states[TID]=LocalStates;
	data[TID].PtNullMCCInCell = Null;
}
__device__ void O2_O2_ion(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, O2CollD *CX, int TnRct, float*MCCR,GGA *BG){
	int i,j,k,n,index,index2,index3;
	int ID,PNMC,MPNC,Null,Flag;
	int Target,oldPNC;
	int Colltype;
	float mofm,R1,R2;
	float VX,VY,VZ;
	float engy,dum,vel;	
	float SumSigma,SumEngyLoss;
	float vneutx,vneuty,vneutz;

	ID = TID%Gsize;
    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	Null = 0;
    // Calculate total Collision probability
	i = info[1].St_num + ID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
        // Calculate energy
		n = (nvel-1)*curand_uniform(&LocalStates);
		dev_maxwellv(&vneutx,&vneuty,&vneutz,vsave[n],BG[ID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
		VX = sp[i].vx;
		VY = sp[i].vy;
		VZ = sp[i].vz;
        Flag = sp[i].CellID;
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
        // 2 : Charge exchange O2+
		// 3 : Charge exchange O+
        switch(Flag){
			case 0:{
				mofm = info_CX[47].mofM;
                R1 = curand_uniform(&LocalStates) * sigv[12].val / vel;
				if(engy > info_CX[47].Th_e &&R1<=(SumSigma=O2_CrossSection(47, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+47]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[48].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[13].val / vel;
				if(engy > info_CX[48].Th_e &&R1<=(SumSigma=O2_CrossSection(48, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+48]++;
				}else if(engy > info_CX[49].Th_e &&R1<=(SumSigma=O2_CrossSection(49, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1; 
					MCCR[ID*TnRct+49]++;
				}else if(engy > info_CX[50].Th_e &&R1<=(SumSigma=O2_CrossSection(50, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+50]++;
				}else{
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[51].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[14].val / vel;
                if(engy > info_CX[51].Th_e &&R1<=(SumSigma=O2_CrossSection(51, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+51]++;
				}else{
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[52].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[15].val / vel;
                if(engy > info_CX[52].Th_e && R1<=(SumSigma=O2_CrossSection(52, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+52]++;
				}else{
					Null++;
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
        case 3: // 3 : Charge exchange o+
			oldPNC = atomicAdd(&data[TID+Gsize].PtNumInCell,1);
			index = info[2].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID+Gsize;
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
	data[TID].PtNullMCCInCell = Null;
}
__device__ void O2_O_ion(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, O2CollD *CX, int TnRct, float*MCCR,GGA *BG){
	int i,j,k,n,index,index2,index3;
	int ID,PNMC,MPNC,Null,Flag;
	int Target,oldPNC;
	int Colltype;
	float mofm,R1,R2;
	float VX,VY,VZ;
	float engy,dum,vel;	
	float SumSigma,SumEngyLoss;
	float vneutx,vneuty,vneutz;

	ID = TID%Gsize;
    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	Null = 0;
    // Calculate total Collision probability
	i = info[2].St_num + ID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
        // Calculate energy
		n = (nvel-1)*curand_uniform(&LocalStates);
		dev_maxwellv(&vneutx,&vneuty,&vneutz,vsave[n],BG[ID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
		VX = sp[i].vx;
		VY = sp[i].vy;
		VZ = sp[i].vz;
        Flag = sp[i].CellID;
        VX = sp[i].vx - vneutx;
		VY = sp[i].vy - vneuty;
		VZ = sp[i].vz - vneutz;
		dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		engy = info[2].Escale * dum;
        Colltype = 0;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Scattering
        // 2 : Charge exchange O2+
		// 3 : Charge exchange O+
        switch(Flag){
			case 0:{
				mofm = info_CX[53].mofM;
                R1 = curand_uniform(&LocalStates) * sigv[16].val / vel;
				if(engy > info_CX[53].Th_e &&R1<=(SumSigma=O2_CrossSection(53, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+53]++;
				}else if(engy > info_CX[54].Th_e &&R1<=(SumSigma=O2_CrossSection(54, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1; 
					MCCR[ID*TnRct+54]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[55].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[17].val / vel;
				if(engy > info_CX[55].Th_e &&R1<=(SumSigma=O2_CrossSection(55, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+55]++;
				}else{
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[56].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[18].val / vel;
                if(engy > info_CX[56].Th_e &&R1<=(SumSigma=O2_CrossSection(56, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+56]++;
				}else{
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[57].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[19].val / vel;
                if(engy > info_CX[57].Th_e && R1<=(SumSigma=O2_CrossSection(57, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+57]++;
				}else{
					Null++;
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
			index = info[2].St_num + ID + oldPNC*Gsize;
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
			index = info[2].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX+vneutx;
			sp[index].vy = VY+vneuty;
			sp[index].vz = VZ+vneutz;
            break;
        case 2: // 2 : Charge exchange o2+
			oldPNC = atomicAdd(&data[TID-Gsize].PtNumInCell,1);
			index = info[1].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID-Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = vneutx;
			sp[index].vy = vneuty;
			sp[index].vz = vneutz;
			break;
        case 3: // 3 : Charge exchange o+
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[2].St_num + ID + oldPNC*Gsize;
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
	data[TID].PtNullMCCInCell = Null;
}
__device__ void O2_O_negative(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, O2CollD *CX, int TnRct, float*MCCR, GGA *BG){
	int i,j,k,n,index,index2,index3;
	int ID,PNMC,MPNC,Null,Flag;
	int Target,oldPNC;
	int Colltype;
	float mofm,R1,R2;
	float VX,VY,VZ;
	float engy,dum,vel;
	int Iz_isp1,Iz_isp2;		
	float SumSigma,SumEngyLoss;

	ID = TID%Gsize;
    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	Null = 0;
    // Calculate total Collision probability
	i = info[3].St_num + ID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
        // Calculate energy
		VX = sp[i].vx;
		VY = sp[i].vy;
		VZ = sp[i].vz;
        Flag = sp[i].CellID;
		//if(Flag !=0) printf("\n[%d] : Flag = %d \n\n",TID,Flag);
        dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		VX/=vel; VY/=vel; VZ/=vel;
		engy = info[3].Escale * dum;
        Colltype = 0;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Scattering
        // 2 : Detachment using maxwellv
		// 3 : dissociative  recombination just delete
        switch(Flag){
			case 0:{
				mofm = info_CX[41].mofM;
                R1 = curand_uniform(&LocalStates) * sigv[9].val / vel;
				if(engy > info_CX[41].Th_e &&R1<=(SumSigma=O2_CrossSection(41, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1; 
					MCCR[ID*TnRct+41]++;
				}else if(engy > info_CX[42].Th_e && R1<=(SumSigma += O2_CrossSection(42, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+42]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[43].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[10].val / vel;
				if(engy > info_CX[43].Th_e &&R1<=(SumSigma=O2_CrossSection(43, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+43]++;
				}else{
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[44].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[11].val / vel;
                if(engy > info_CX[44].Th_e &&R1<=(SumSigma=O2_CrossSection(44, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+44]++;
				}else{
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[45].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[12].val / vel;
                if(engy > info_CX[45].Th_e && R1<=(SumSigma=O2_CrossSection(45, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+45]++;
				}else{
					Null++;
				}
                break;
			}
			case 4:{
				mofm = info_CX[46].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[13].val / vel;
                if(engy > info_CX[46].Th_e &&R1<=(SumSigma=O2_CrossSection(46, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+46]++;
				}else{
					Null++;
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
			//printf("\n[%d][%d] : oldPNC = %d \n\n",TID,Gsize,oldPNC);
			index = info[3].St_num + ID + oldPNC*Gsize;
			//printf("0[%d][%d]: %g,%g,%g,%g,%g,\n",TID,ID,sp[i].x,sp[i].y,sp[i].vx,sp[i].vy,sp[i].vz);
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = sp[i].vx;
			sp[index].vy = sp[i].vy;
			sp[index].vz = sp[i].vz;
            break;
        case 1: // 1 : Scattering
			//printf("1[%d][%d]: %g,%g,%g,%g,%g,\n",TID,ID,sp[i].x,sp[i].y,sp[i].vx,sp[i].vy,sp[i].vz);
			dev_newvel_IONSC(&VX,&VY,&VZ,vel,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[3].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
            break;
        case 2: // 2 : Detachment using maxwellv
			//printf("2[%d][%d]: %g,%g,%g,%g,%g,\n",TID,ID,sp[i].x,sp[i].y,sp[i].vx,sp[i].vy,sp[i].vz);
			oldPNC = atomicAdd(&data[ID].PtNumInCell,1);
			index = info[0].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = ID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
			break;
        case 3: // 3 : dissociative  recombination just delete
            break;
        default:
            break;
        }
		i-=Gsize;
	}
	states[TID]=LocalStates;
	data[TID].PtNullMCCInCell = Null;
}
__device__ void  O2Collision_Check(int Gsize, int Csize, int ngy, int TID, float dt, int MCCn, float dtm, float dx, float dy,
                                        curandState *states, Species *info, GPG *data, GCP *sp, MCC_sigmav *sigv, GGA *BG, GFC *Fluid){
	int i,j,k,index,Randn;
	int ID,isp,CID,PNMC,MPNC;
    int PNC,Flag;
	int nx,ny,ngx;
	float Tprob,Prob1,Prob2,Prob3,Prob4,Prob5,Prob6,Prob7;
	float R1;
	ID = TID%Gsize;
    isp = TID/Gsize;
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
    switch (isp){
    case 0: // Electron
        Prob1 = 1.0f - exp(-1*dtm*sigv[0].val*BG[ID].BackDen1);
	    Prob2 = Prob1 + 1.0f - exp(-1*dtm*sigv[1].val*Fluid[CID].ave_den);
        Prob3 = Prob2 + 1.0f - exp(-1*dtm*sigv[2].val*Fluid[CID+Csize].ave_den);
        Prob4 = Prob3 + 1.0f - exp(-1*dtm*sigv[3].val*data[ID+3*Gsize].den*info[3].np2c*dx*dy);
        Prob5 = Prob4 + 1.0f - exp(-1*dtm*sigv[4].val*data[ID+Gsize].den*info[1].np2c*dx*dy);
        Prob6 = Prob5 + 1.0f - exp(-1*dtm*sigv[5].val*Fluid[CID+2*Csize].ave_den);
        Prob7 = Prob6 + 1.0f - exp(-1*dtm*sigv[6].val*Fluid[CID+3*Csize].ave_den);
	    Tprob = Prob7;
		Randn = MCCn;
        break;
    case 1: // O2+
        Prob1 = 1.0 - exp(-1*dtm*sigv[12].val*Fluid[CID+2*Csize].ave_den);
	    Prob2 = Prob1 + 1.0 - exp(-1*dtm*sigv[13].val*BG[ID].BackDen1);
	    Prob3 = Prob2 + 1.0 - exp(-1*dtm*sigv[14].val*Fluid[CID].ave_den);
	    Prob4 = Prob3 + 1.0 - exp(-1*dtm*sigv[15].val*Fluid[CID+Csize].ave_den);
        Tprob = Prob4;
		Randn = MCCn;
        break;
    case 2: // O+
        Prob1 = 1.0 - exp(-1*dt*sigv[16].val*BG[ID].BackDen1);
	    Prob2 = Prob1 + 1.0 - exp(-1*dt*sigv[17].val*Fluid[CID+2*Csize].ave_den);
	    Prob3 = Prob2 + 1.0 - exp(-1*dt*sigv[18].val*Fluid[CID].ave_den);
	    Prob4 = Prob3 + 1.0 - exp(-1*dt*sigv[19].val*Fluid[CID+Csize].ave_den);
	    Tprob = Prob4;
		Randn = 1;
        break;
    case 3: // O-
        Prob1 = 1.0 - exp(-1*dt*sigv[7].val*BG[ID].BackDen1);
	    Prob2 = Prob1 + 1.0 - exp(-1*dt*sigv[8].val*Fluid[CID+2*Csize].ave_den);
	    Prob3 = Prob2 + 1.0 - exp(-1*dt*sigv[9].val*data[ID+Gsize].den*info[1].np2c*dx*dy);
	    Prob4 = Prob3 + 1.0 - exp(-1*dt*sigv[10].val*data[ID+2*Gsize].den*info[2].np2c*dx*dy);
	    Prob5 = Prob4 + 1.0 - exp(-1*dt*sigv[11].val*Fluid[CID].ave_den);
	    Tprob = Prob5;
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
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
                else if(R1 <= Prob4)	Flag = (int)3;
                else if(R1 <= Prob5)	Flag = (int)4;
                else if(R1 <= Prob6)	Flag = (int)5;
		        else			        Flag = (int)6;
                break;
            case 1:
                if(R1 <= Prob1)	        Flag = (int)0;
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
		        else			        Flag = (int)3;
                break;
            case 2:
                if(R1 <= Prob1)	        Flag = (int)0;
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
		        else			        Flag = (int)3;
                break;
            case 3:
                if(R1 <= Prob1)	        Flag = (int)0;
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
                else if(R1 <= Prob4)	Flag = (int)3;
		        else			        Flag = (int)4;
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
__device__ float O2_CrossSection(int R, float engy, int N_LOGX, float idLOGX, O2CollD *data){
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
            printf("\nError : Call about cross section data in O2MCC.\n\n");
            return 0.0;
    }
}