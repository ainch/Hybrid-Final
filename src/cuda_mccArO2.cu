#include "cuda_mccArO2.cuh"
__device__ void ArO2_Electron(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
	int i,j,k,n,index,index2,index3;
	int PNMC,MPNC,Null,Flag;
	int Target,oldPNC;
	int Colltype;
	float mofm,R1,R2;
	float VX,VY,VZ;
	float rengy,engy,dum,vel,vel2;
	int Iz_isp1,Iz_isp2;		
	float SumSigma,SumEngyLoss;

    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	Null = 0;
	// Calculate total Collision probability
	i = info[0].St_num + TID + (MPNC-1)*Gsize;
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
		//printf("[%d]:[%3.3e, %3.3e, %3.3e] = [%g], case[%d] engy = %g \n",k,sp[i].vx,sp[i].vy,sp[i].vz,dum,Flag,engy);
        Colltype = 1;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Energy loss
        // 2 : Attachment using maxwellv
		// 3 : ionization 1 -> 3 charged
		// 4 : dissociative  recombination just delete
		// 5 : Detachment 
		R1 = curand_uniform(&LocalStates);
        switch(Flag){
			case 0:{ // E + Ar
				mofm = info_CX[0].mofM;
				R1 *= sigv[0].val / vel;
				//printf("k[%d]:c[%d]:engy[%g]:R1[%g] \n",k,Flag,engy,R1);
				//printf("k[%d]:CX0[%1.3e]\n",k,ArO2_CrossSection(0, engy, N_LOGX, idLOGX, CX));
				if(engy > info_CX[0].Th_e && R1<=(SumSigma=ArO2_CrossSection(0, engy, N_LOGX, idLOGX, CX))){
					MCCR[TID*TnRct]++;
				}else if(engy > info_CX[1].Th_e && R1<=(SumSigma += ArO2_CrossSection(1, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[1].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+1]++;
				}else if(engy > info_CX[2].Th_e && R1<=(SumSigma += ArO2_CrossSection(2, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[2].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+2]++;
				}else if(engy > info_CX[3].Th_e && R1<=(SumSigma += ArO2_CrossSection(3, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3;
					engy-=info_CX[3].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 1;
					MCCR[TID*TnRct+3]++;
				}else{
					Colltype = 0;
					Null++;
				}
				
				break;
			}
			case 1:{ // E + Ar*
				mofm = info_CX[4].mofM;
				R1 *= sigv[1].val / vel;
				if(engy > info_CX[4].Th_e && R1<=(SumSigma=ArO2_CrossSection(4, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3;
					engy-=info_CX[4].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 1;
					MCCR[TID*TnRct+4]++;
				}else{
					Colltype = 0;
					Null++;
				}
				break;
			}
			case 2:{ // E + O2
				mofm = info_CX[5].mofM;
                R1 *= sigv[2].val / vel;
				if(engy > info_CX[5].Th_e &&R1<=(SumSigma=ArO2_CrossSection(5, engy, N_LOGX, idLOGX, CX))){
                    // R0 Elastic
					MCCR[TID*TnRct+5]++;
				}else if(engy > info_CX[6].Th_e && R1<=(SumSigma += ArO2_CrossSection(6, engy, N_LOGX, idLOGX, CX))){
                    //"6.e+O2>e+O2*");
					engy-=info_CX[6].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+6]++;
				}else if(engy > info_CX[7].Th_e && R1<=(SumSigma += ArO2_CrossSection(7, engy, N_LOGX, idLOGX, CX))){
                    //"7.e+O2>e+O2*");
					engy-=info_CX[7].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+7]++;
				}else if(engy > info_CX[8].Th_e && R1<=(SumSigma += ArO2_CrossSection(8, engy, N_LOGX, idLOGX, CX))){
                    //"8.e+O2>e+O2A");
					engy-=info_CX[8].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+8]++;
				}else if(engy > info_CX[9].Th_e && R1<=(SumSigma += ArO2_CrossSection(9, engy, N_LOGX, idLOGX, CX))){
					//"9.e+O2>e+O2B");
                    engy-=info_CX[9].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+9]++;
				}else if(engy > info_CX[10].Th_e && R1<=(SumSigma += ArO2_CrossSection(10, engy, N_LOGX, idLOGX, CX))){
                    //"10.e+O2>e+O2*");
					engy-=info_CX[10].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+10]++;
				}else if(engy > info_CX[11].Th_e && R1<=(SumSigma += ArO2_CrossSection(11, engy, N_LOGX, idLOGX, CX))){
                    //"11.e+O2>OP+O-"
                    Colltype = 2;
					MCCR[TID*TnRct+11]++;
                }else if(engy > info_CX[12].Th_e && R1<=(SumSigma += ArO2_CrossSection(12, engy, N_LOGX, idLOGX, CX))){
                    //"12.e+O2>e+2OP");
					engy-=info_CX[12].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+12]++;
                }else if(engy > info_CX[13].Th_e && R1<=(SumSigma += ArO2_CrossSection(13, engy, N_LOGX, idLOGX, CX))){
                    //"13.e+O2>e+OP+OD");
					engy-=info_CX[13].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+13]++;
                }else if(engy > info_CX[14].Th_e && R1<=(SumSigma += ArO2_CrossSection(14, engy, N_LOGX, idLOGX, CX))){
                    //"14.e+O2>e+2OD");
					engy-=info_CX[14].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+14]++;
                }else if(engy > info_CX[15].Th_e && R1<=(SumSigma += ArO2_CrossSection(15, engy, N_LOGX, idLOGX, CX))){
                    //"15.e+O2>2e+O2^");
                    Colltype = 3;
					engy-=info_CX[15].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+15]++;
                }else if(engy > info_CX[16].Th_e && R1<=(SumSigma += ArO2_CrossSection(16, engy, N_LOGX, idLOGX, CX))){
                    //"16.e+O2>e+OP+O*");
					engy-=info_CX[11].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+16]++;
                }else if(engy > info_CX[17].Th_e && R1<=(SumSigma += ArO2_CrossSection(17, engy, N_LOGX, idLOGX, CX))){
                    //"17.e+O2>e+O^+O-");
                    Colltype = 3;
					engy-=info_CX[17].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 3;
					Iz_isp2 = 4;
					MCCR[TID*TnRct+17]++;
                }else if(engy > info_CX[18].Th_e && R1<=(SumSigma += ArO2_CrossSection(18, engy, N_LOGX, idLOGX, CX))){
                    //"18.e+O2>2e+O^+OP");  
                    Colltype = 3;
					engy-=info_CX[18].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+18]++;
				}else{
					Colltype = 0;
					Null++;
				}
				break;
			}
			case 3:{ // E + O2A
				mofm = info_CX[19].mofM;
                R1 *= sigv[3].val / vel;
				if(engy > info_CX[19].Th_e &&R1<=(SumSigma=ArO2_CrossSection(19, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"19.e+O2A>2e+O2+");
					engy-=info_CX[19].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+19]++;
				}else if(engy > info_CX[20].Th_e && R1<=(SumSigma += ArO2_CrossSection(20, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; //"20.e+O2A>OP+O-");
					MCCR[TID*TnRct+20]++;
				}else if(engy > info_CX[21].Th_e && R1<=(SumSigma += ArO2_CrossSection(21, engy, N_LOGX, idLOGX, CX))){
					engy+=0.977f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+21]++;
				}else if(engy > info_CX[22].Th_e && R1<=(SumSigma += ArO2_CrossSection(22, engy, N_LOGX, idLOGX, CX))){
					engy+=0.977f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+22]++;
				}else if(engy > info_CX[23].Th_e && R1<=(SumSigma += ArO2_CrossSection(23, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[23].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+23]++;
				}else if(engy > info_CX[24].Th_e && R1<=(SumSigma += ArO2_CrossSection(24, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[24].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+24]++;
				}else if(engy > info_CX[25].Th_e && R1<=(SumSigma += ArO2_CrossSection(25, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[25].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+25]++;
				}else if(engy > info_CX[26].Th_e && R1<=(SumSigma += ArO2_CrossSection(26, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"26.e+O2A>2e+O^+OP");
					engy-=info_CX[26].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+26]++;
				}else{
					Colltype = 0;
					Null++;
				}
				break;
			}
			case 4:{ // E + O2B
				mofm = info_CX[27].mofM;
                R1 *= sigv[4].val / vel;
                if(engy > info_CX[27].Th_e && R1<=(SumSigma=ArO2_CrossSection(27, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"22.e+O2B>2e+O2^");
					engy-=info_CX[27].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+27]++;
				}else if(engy > info_CX[28].Th_e && R1<=(SumSigma += ArO2_CrossSection(28, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; //"23.e+O2B>OP+O-");
					MCCR[TID*TnRct+28]++;
				}else if(engy > info_CX[29].Th_e && R1<=(SumSigma += ArO2_CrossSection(29, engy, N_LOGX, idLOGX, CX))){
					engy+=1.627f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+29]++;
				}else if(engy > info_CX[30].Th_e && R1<=(SumSigma += ArO2_CrossSection(30, engy, N_LOGX, idLOGX, CX))){
					engy+=1.627f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+30]++;
				}else if(engy > info_CX[31].Th_e && R1<=(SumSigma += ArO2_CrossSection(31, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[31].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+31]++;
				}else if(engy > info_CX[32].Th_e && R1<=(SumSigma += ArO2_CrossSection(32, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[32].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+32]++;
				}else if(engy > info_CX[33].Th_e && R1<=(SumSigma += ArO2_CrossSection(33, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[33].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+33]++;
				}else if(engy > info_CX[34].Th_e && R1<=(SumSigma += ArO2_CrossSection(34, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"34.e+O2B>2e+O^+OP");
					engy-=info_CX[34].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+34]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 5:{ // E + O-
				mofm = info_CX[35].mofM;
                R1 *= sigv[5].val / vel;
                if(engy > info_CX[35].Th_e && R1<=(SumSigma=ArO2_CrossSection(35, engy, N_LOGX, idLOGX, CX))){
					Colltype = 5; //"30.e+O->2e+OP");
					engy-=info_CX[35].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+35]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 6:{ // E + O2+
				mofm = info_CX[36].mofM;
                R1 *= sigv[6].val / vel;
                if(engy > info_CX[36].Th_e &&R1<=(SumSigma=ArO2_CrossSection(36, engy, N_LOGX, idLOGX, CX))){
					Colltype = 4; //"36.e+O2^>OP+OD");
					MCCR[TID*TnRct+36]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 7:{ // E + OP
				mofm = info_CX[37].mofM;
                R1 *= sigv[7].val / vel;
                if(engy > info_CX[37].Th_e &&R1<=(SumSigma=ArO2_CrossSection(37, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[37].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+37]++;
				}else if(engy > info_CX[38].Th_e && R1<=(SumSigma += ArO2_CrossSection(38, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[38].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+38]++;
				}else if(engy > info_CX[39].Th_e && R1<=(SumSigma += ArO2_CrossSection(39, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[39].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+39]++;
				}else if(engy > info_CX[40].Th_e && R1<=(SumSigma += ArO2_CrossSection(40, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[40].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+40]++;
				}else if(engy > info_CX[41].Th_e && R1<=(SumSigma += ArO2_CrossSection(41, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[41].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+41]++;
				}else if(engy > info_CX[42].Th_e && R1<=(SumSigma += ArO2_CrossSection(42, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"42.e+OP>2e+O^");
					engy-=info_CX[42].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+42]++;
				}else if(engy > info_CX[43].Th_e && R1<=(SumSigma += ArO2_CrossSection(43, engy, N_LOGX, idLOGX, CX))){
					engy-=info_CX[43].Th_e;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+43]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 8:{ // E + OD
				mofm = info_CX[44].mofM;
                R1 *= sigv[8].val / vel;
                if(engy > info_CX[44].Th_e &&R1<=(SumSigma=ArO2_CrossSection(44, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; //"44.e+OD>2e+O^");
					engy-=info_CX[44].Th_e;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+44]++;
				}else if(engy > info_CX[45].Th_e && R1<=(SumSigma += ArO2_CrossSection(45, engy, N_LOGX, idLOGX, CX))){
					engy+=1.96f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+45]++;
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
			oldPNC = atomicAdd(&data[TID+4*Gsize].PtNumInCell,1);
			index = info[4].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID+4*Gsize;
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
				dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
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
				dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			}else{
				n = (nvel-1)*curand_uniform(&LocalStates);
				dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			}
			// energy loss electron 
			dev_anewvel(engy,vel2,&VX,&VY,&VZ,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
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
			oldPNC = atomicAdd(&data[TID+4*Gsize].PtNumInCell,0);
			if(oldPNC>1){
				R2 = curand_uniform(&LocalStates);
				index = info[4].St_num + TID + oldPNC*Gsize;
				index2 = (int)((float)oldPNC * R2);
				index3 = info[4].St_num + TID + index2*Gsize;
				sp[index3].CellID = sp[index].CellID;
				sp[index3].x = sp[index].x;
				sp[index3].y = sp[index].y;
				sp[index3].vx = sp[index].vx;
				sp[index3].vy = sp[index].vy;
				sp[index3].vz = sp[index].vz;
				atomicAdd(&data[TID+4*Gsize].PtNumInCell,-1);
				//printf("2[%d][%d]: %g,%g,%g,%g,%g,\n",TID,sp[index].x,sp[index].y,sp[index].vx,sp[index].vy,sp[i].vz);
			}else if(oldPNC == 1){
				atomicAdd(&data[TID+4*Gsize].PtNumInCell,-1);
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
__device__ void ArO2_Ar_ion(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
	int i,j,k,n,index,index2,index3;
	int ID,PNMC,MPNC,Null,Flag;
	int Target,oldPNC;
	int Colltype;
	float mofm,R1,R2;
	float VX,VY,VZ;
	float engy,dum,vel;	
	float SumSigma,SumEngyLoss;
	float vneut,vneutx,vneuty,vneutz;

	ID = TID%Gsize;
    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	Null = 0;
    // Calculate total Collision probability
	i = info[1].St_num + ID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
        // Calculate energy
		Flag = sp[i].CellID;
		if(Flag == 0) vneut = BG[ID].BackVel1;
		else  vneut = BG[ID].BackVel2;
		n = (nvel-1)*curand_uniform(&LocalStates);
		dev_maxwellv(&vneutx,&vneuty,&vneutz,vsave[n],vneut,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
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
		R1 = curand_uniform(&LocalStates);
        switch(Flag){
			case 0:{
				mofm = info_CX[65].mofM;
        		R1 = curand_uniform(&LocalStates)*sigv[23].val / vel;
				if(engy > info_CX[65].Th_e &&R1<=(SumSigma=ArO2_CrossSection(65, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+65]++;
				}else if(engy > info_CX[66].Th_e &&R1<=(SumSigma=ArO2_CrossSection(66, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1; 
					MCCR[ID*TnRct+66]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[67].mofM;
        		R1 = curand_uniform(&LocalStates)*sigv[24].val / vel;
				if(engy > info_CX[67].Th_e &&R1<=(SumSigma=ArO2_CrossSection(67, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1; 
					MCCR[ID*TnRct+67]++;
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
        case 2: // 2 : Charge exchange 
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
	data[TID].PtNullMCCInCell = Null;
}
__device__ void ArO2_O2_ion(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
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
		// 4 : Charge exchange AR+
        switch(Flag){
			case 0:{
				mofm = info_CX[52].mofM;
                R1 = curand_uniform(&LocalStates) * sigv[14].val / vel;
				if(engy > info_CX[52].Th_e &&R1<=(SumSigma=ArO2_CrossSection(52, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+52]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[53].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[15].val / vel;
				if(engy > info_CX[53].Th_e &&R1<=(SumSigma=ArO2_CrossSection(53, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+53]++;
				}else if(engy > info_CX[54].Th_e &&R1<=(SumSigma=ArO2_CrossSection(54, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1; 
					MCCR[ID*TnRct+54]++;
				}else if(engy > info_CX[55].Th_e &&R1<=(SumSigma=ArO2_CrossSection(55, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+55]++;
				}else{
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[56].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[16].val / vel;
                if(engy > info_CX[56].Th_e &&R1<=(SumSigma=ArO2_CrossSection(56, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+56]++;
				}else{
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[57].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[17].val / vel;
                if(engy > info_CX[57].Th_e && R1<=(SumSigma=ArO2_CrossSection(57, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+57]++;
				}else{
					Null++;
				}
                break;
			}
			case 4:{
				mofm = info_CX[58].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[18].val / vel;
                if(engy > info_CX[58].Th_e && R1<=(SumSigma=ArO2_CrossSection(58, engy, N_LOGX, idLOGX, CX))){
					Colltype = 4; 
					MCCR[ID*TnRct+58]++;
				}else if(engy > info_CX[59].Th_e && R1<=(SumSigma=ArO2_CrossSection(59, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+59]++;
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
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[2].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = vneutx;
			sp[index].vy = vneuty;
			sp[index].vz = vneutz;
			break;
        case 3: // 3 : Charge exchange o+
			oldPNC = atomicAdd(&data[TID+Gsize].PtNumInCell,1);
			index = info[3].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID+Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = vneutx;
			sp[index].vy = vneuty;
			sp[index].vz = vneutz;
            break;
		case 4: // 4 : Charge exchange AR+
			oldPNC = atomicAdd(&data[TID-Gsize].PtNumInCell,1);
			index = info[1].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID-Gsize;
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
__device__ void ArO2_O_ion(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
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
	i = info[3].St_num + ID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
        // Calculate energy
		n = (nvel-1)*curand_uniform(&LocalStates);
		dev_maxwellv(&vneutx,&vneuty,&vneutz,vsave[n],BG[ID].BackVel2,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
		VX = sp[i].vx;
		VY = sp[i].vy;
		VZ = sp[i].vz;
        Flag = sp[i].CellID;
        VX = sp[i].vx - vneutx;
		VY = sp[i].vy - vneuty;
		VZ = sp[i].vz - vneutz;
		dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		engy = info[3].Escale * dum;
        Colltype = 0;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Scattering
        // 2 : Charge exchange O2+
		// 3 : Charge exchange O+
        switch(Flag){
			case 0:{
				mofm = info_CX[60].mofM;
                R1 = curand_uniform(&LocalStates) * sigv[19].val / vel;
				if(engy > info_CX[60].Th_e &&R1<=(SumSigma=ArO2_CrossSection(60, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+60]++;
				}else if(engy > info_CX[61].Th_e &&R1<=(SumSigma=ArO2_CrossSection(61, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1; 
					MCCR[ID*TnRct+61]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[62].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[20].val / vel;
				if(engy > info_CX[62].Th_e &&R1<=(SumSigma=ArO2_CrossSection(62, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+62]++;
				}else{
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[63].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[21].val / vel;
                if(engy > info_CX[63].Th_e &&R1<=(SumSigma=ArO2_CrossSection(63, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+63]++;
				}else{
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[64].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[22].val / vel;
                if(engy > info_CX[64].Th_e && R1<=(SumSigma=ArO2_CrossSection(64, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+64]++;
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
			index = info[3].St_num + ID + oldPNC*Gsize;
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
			index = info[3].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX+vneutx;
			sp[index].vy = VY+vneuty;
			sp[index].vz = VZ+vneutz;
            break;
        case 2: // 2 : Charge exchange o2+
			oldPNC = atomicAdd(&data[TID-Gsize].PtNumInCell,1);
			index = info[2].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID-Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = vneutx;
			sp[index].vy = vneuty;
			sp[index].vz = vneutz;
			break;
        case 3: // 3 : Charge exchange o+
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[3].St_num + ID + oldPNC*Gsize;
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
__device__ void ArO2_O_negative(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
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
	i = info[4].St_num + ID + (MPNC-1)*Gsize;
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
		engy = info[4].Escale * dum;
        Colltype = 0;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Scattering
        // 2 : Detachment using maxwellv
		// 3 : dissociative  recombination just delete
        switch(Flag){
			case 0:{
				mofm = info_CX[46].mofM;
                R1 = curand_uniform(&LocalStates) * sigv[9].val / vel;
				if(engy > info_CX[46].Th_e &&R1<=(SumSigma=ArO2_CrossSection(46, engy, N_LOGX, idLOGX, CX))){
					Colltype = 1; 
					MCCR[ID*TnRct+46]++;
				}else if(engy > info_CX[47].Th_e && R1<=(SumSigma += ArO2_CrossSection(47, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+47]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[48].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[10].val / vel;
				if(engy > info_CX[48].Th_e &&R1<=(SumSigma=ArO2_CrossSection(48, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+48]++;
				}else{
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[49].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[11].val / vel;
                if(engy > info_CX[49].Th_e &&R1<=(SumSigma=ArO2_CrossSection(49, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+49]++;
				}else{
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[50].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[12].val / vel;
                if(engy > info_CX[50].Th_e && R1<=(SumSigma=ArO2_CrossSection(50, engy, N_LOGX, idLOGX, CX))){
					Colltype = 3; 
					MCCR[ID*TnRct+50]++;
				}else{
					Null++;
				}
                break;
			}
			case 4:{
				mofm = info_CX[51].mofM;
                R1 = curand_uniform(&LocalStates)*sigv[13].val / vel;
                if(engy > info_CX[51].Th_e &&R1<=(SumSigma=ArO2_CrossSection(51, engy, N_LOGX, idLOGX, CX))){
					Colltype = 2; 
					MCCR[ID*TnRct+51]++;
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
			index = info[4].St_num + ID + oldPNC*Gsize;
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
			index = info[4].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
            break;
        case 2: // 2 : Detachment using maxwellv
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
__device__ void  ArO2_Collision_Check(int Gsize, int Csize, int ngy, int TID, float dt, int MCCn, float dtm, float dx, float dy,
                                        curandState *states, Species *info, GPG *data, GCP *sp, MCC_sigmav *sigv, GGA *BG, GFC *Fluid){
	int i,j,k,index,Randn;
	int ID,isp,CID,PNMC,MPNC;
    int PNC,Flag;
	int nx,ny,ngx;
	float Tprob,Prob1,Prob2,Prob3,Prob4,Prob5,Prob6,Prob7,Prob8,Prob9;
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
		Prob1 = 1.0f - exp(-1*dtm*sigv[0].val*BG[ID].BackDen1);  // E + Ar
		Prob2 = Prob1 + 1.0f - exp(-1*dtm*sigv[1].val*Fluid[CID].ave_den);  // E + Ar*
        Prob3 = Prob2 + 1.0f - exp(-1*dtm*sigv[2].val*BG[ID].BackDen2);  // E + O2
	    Prob4 = Prob3 + 1.0f - exp(-1*dtm*sigv[3].val*Fluid[CID+Csize].ave_den);  // E + O2A
        Prob5 = Prob4 + 1.0f - exp(-1*dtm*sigv[4].val*Fluid[CID+2*Csize].ave_den);  // E + O2B
        Prob6 = Prob5 + 1.0f - exp(-1*dtm*sigv[5].val*data[ID+4*Gsize].den*info[4].np2c*dx*dy);  // E + O-
        Prob7 = Prob6 + 1.0f - exp(-1*dtm*sigv[6].val*data[ID+2*Gsize].den*info[2].np2c*dx*dy);  // E + O2+
        Prob8 = Prob7 + 1.0f - exp(-1*dtm*sigv[7].val*Fluid[CID+3*Csize].ave_den);  // E + OP
        Prob9 = Prob8 + 1.0f - exp(-1*dtm*sigv[8].val*Fluid[CID+4*Csize].ave_den);  // E + OD
	    Tprob = Prob9; 
		//printf("Case[%d] : 1[%1.2e] 2[%1.2e] 3[%1.2e] 4[%1.2e] 5[%1.2e] 6[%1.2e] 7[%1.2e] 8[%1.2e] 9[%1.2e]\n",isp,Prob1, Prob2, Prob3, Prob4, Prob5, Prob6, Prob7, Prob8, Prob9);
		//printf("Case[%d] : Ar[%1.2e] O2[%1.2e]\n",isp,BG[ID].BackDen1,BG[ID].BackDen2);
		//printf("Case[%d] : Ar*[%1.2e] O2A[%1.2e] O2B[%1.2e] OP[%1.2e] OD[%1.2e]\n",isp,Fluid[CID].ave_den,Fluid[CID+1*Csize].ave_den,Fluid[CID+2*Csize].ave_den,Fluid[CID+3*Csize].ave_den,Fluid[CID+4*Csize].ave_den);
		Randn = MCCn;
        break;
	case 1: // Ar+
		Prob1 = 1.0 - exp(-1*dt*sigv[23].val*BG[ID].BackDen1); // Ar+ + Ar
	    Prob2 = Prob1 + 1.0 - exp(-1*dt*sigv[24].val*BG[ID].BackDen2); // Ar+ + O2
		Tprob = Prob2;
		//printf("Case[%d] : 1[%1.2e] 2[%1.2e]\n",isp,Prob1, Prob2);
		Randn = 1;
		break;
    case 2: // O2+
        Prob1 = 1.0 - exp(-1*dt*sigv[14].val*Fluid[CID+3*Csize].ave_den); // O2+ + OP
	    Prob2 = Prob1 + 1.0 - exp(-1*dt*sigv[15].val*BG[ID].BackDen2); // O2+ + O2
	    Prob3 = Prob2 + 1.0 - exp(-1*dt*sigv[16].val*Fluid[CID+Csize].ave_den); // O2+ + O2A
	    Prob4 = Prob3 + 1.0 - exp(-1*dt*sigv[17].val*Fluid[CID+2*Csize].ave_den); // O2+ + O2B
		Prob5 = Prob4 + 1.0 - exp(-1*dt*sigv[18].val*BG[ID].BackDen1); // O2+ + AR
        Tprob = Prob5;
		//printf("Case[%d] : 1[%1.2e] 2[%1.2e] 3[%1.2e] 4[%1.2e] 5[%1.2e]\n",isp,Prob1, Prob2, Prob3, Prob4, Prob5);
		Randn = 1;
        break;
    case 3: // O+
        Prob1 = 1.0 - exp(-1*dt*sigv[19].val*BG[ID].BackDen2); // O+ + O2
	    Prob2 = Prob1 + 1.0 - exp(-1*dt*sigv[20].val*Fluid[CID+3*Csize].ave_den); // O+ + OP
	    Prob3 = Prob2 + 1.0 - exp(-1*dt*sigv[21].val*Fluid[CID+Csize].ave_den); // O+ + O2A
	    Prob4 = Prob3 + 1.0 - exp(-1*dt*sigv[22].val*Fluid[CID+2*Csize].ave_den); // O+ + O2B
	    Tprob = Prob4;
		//printf("Case[%d] : 1[%1.2e] 2[%1.2e] 3[%1.2e] 4[%1.2e]\n",isp,Prob1, Prob2, Prob3, Prob4);
		Randn = 1;
        break;
    case 4: // O-
        Prob1 = 1.0 - exp(-1*dt*sigv[9].val*BG[ID].BackDen2); // O- + O2
	    Prob2 = Prob1 + 1.0 - exp(-1*dt*sigv[10].val*Fluid[CID+3*Csize].ave_den); // O- + OP
	    Prob3 = Prob2 + 1.0 - exp(-1*dt*sigv[11].val*data[ID+2*Gsize].den*info[2].np2c*dx*dy); // O- + O2+
	    Prob4 = Prob3 + 1.0 - exp(-1*dt*sigv[12].val*data[ID+3*Gsize].den*info[3].np2c*dx*dy); // O- + O+
	    Prob5 = Prob4 + 1.0 - exp(-1*dt*sigv[13].val*Fluid[CID+Csize].ave_den); // O- + O2A
	    Tprob = Prob5;
		//printf("Case[%d] : 1[%1.2e] 2[%1.2e] 3[%1.2e] 4[%1.2e] 5[%1.2e]\n",isp,Prob1, Prob2, Prob3, Prob4, Prob5);
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
			//printf("k[%d], PNMC[%d], R1[%g], Tprob[%g]\n",k,PNMC,R1,Tprob);
            switch (isp){
            case 0:
                if(R1 <= Prob1)	        Flag = (int)0;
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
                else if(R1 <= Prob4)	Flag = (int)3;
                else if(R1 <= Prob5)	Flag = (int)4;
                else if(R1 <= Prob6)	Flag = (int)5;
				else if(R1 <= Prob7)	Flag = (int)6;
				else if(R1 <= Prob8)	Flag = (int)7;
		        else			        Flag = (int)8;
				//printf("k[%d], PNMC[%d], R1[%g], Tprob[%g], Flag[%d]\n",k,PNMC,R1,Tprob,Flag);
                break;
            case 1:
                if(R1 <= Prob1)	        Flag = (int)0;
		        else			        Flag = (int)1;
                break;
            case 2:
                if(R1 <= Prob1)	        Flag = (int)0;
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
				else if(R1 <= Prob4)	Flag = (int)3;
		        else			        Flag = (int)4;
                break;
            case 3:
                if(R1 <= Prob1)	        Flag = (int)0;
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
		        else			        Flag = (int)3;
                break;
            case 4:
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
__device__ void ArO2_Electron_TEST(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
	int i,j,k,n,index,index2,index3;
	int PNMC,MPNC,Null,Flag;
	int Target,oldPNC;
	int Colltype;
	float mofm,R1,R2;
	float VX,VY,VZ;
	float rengy,engy,dum,vel,vel2;
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
		R1 = curand_uniform(&LocalStates);
        switch(Flag){
			case 0:{ // E + Ar
				mofm = info_CX[0].mofM;
				if(R1<=0.2){
					MCCR[TID*TnRct]++;
				}else if(R1<=0.4){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+1]++;
				}else if(R1<=0.6){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+2]++;
				}else if(R1<=0.8){
					Colltype = 3;
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 1;
					MCCR[TID*TnRct+3]++;
				}else{
					Colltype = 0;
					Null++;
				}
				
				break;
			}
			case 1:{ // E + Ar*
				mofm = info_CX[4].mofM;
				if(R1<=0.5){
					Colltype = 3;
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 1;
					MCCR[TID*TnRct+4]++;
				}else{
					Colltype = 0;
					Null++;
				}
				break;
			}
			case 2:{ // E + O2
				mofm = info_CX[5].mofM;
				if(R1<=0.0625){
                    // R0 Elastic
					MCCR[TID*TnRct+5]++;
				}else if(R1<=0.125){
                    //"6.e+O2>e+O2*");
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+6]++;
				}else if(R1<=0.1875){
                    //"7.e+O2>e+O2*");
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+7]++;
				}else if(R1<=0.25){
                    //"8.e+O2>e+O2A");
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+8]++;
				}else if(R1<=0.3125){
					//"9.e+O2>e+O2B");
                    engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+9]++;
				}else if(R1<=0.375){
                    //"10.e+O2>e+O2*");
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+10]++;
				}else if(R1<=0.4375){
                    //"11.e+O2>OP+O-"
                    Colltype = 2;
					MCCR[TID*TnRct+11]++;
                }else if(R1<=0.5){
                    //"12.e+O2>e+2OP");
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+12]++;
                }else if(R1<=0.5625){
                    //"13.e+O2>e+OP+OD");
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+13]++;
                }else if(R1<=0.625){
                    //"14.e+O2>e+2OD");
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+14]++;
                }else if(R1<=0.6875){
                    //"15.e+O2>2e+O2^");
                    Colltype = 3;
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+15]++;
                }else if(R1<=0.75){
                    //"16.e+O2>e+OP+O*");
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+16]++;
                }else if(R1<=0.8125){
                    //"17.e+O2>e+O^+O-");
                    Colltype = 3;
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 3;
					Iz_isp2 = 4;
					MCCR[TID*TnRct+17]++;
                }else if(R1<=0.875){
                    //"18.e+O2>2e+O^+OP");  
                    Colltype = 3;
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+18]++;
				}else{
					Colltype = 0;
					Null++;
				}
				break;
			}
			case 3:{ // E + O2A
				mofm = info_CX[19].mofM;
				if(R1<=0.1){
					Colltype = 3; //"19.e+O2A>2e+O2+");
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+19]++;
				}else if(R1<=0.2){
					Colltype = 2; //"20.e+O2A>OP+O-");
					MCCR[TID*TnRct+20]++;
				}else if(R1<=0.3){
					engy+=0.977f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+21]++;
				}else if(R1<=0.4){
					engy+=0.977f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+22]++;
				}else if(R1<=0.5){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+23]++;
				}else if(R1<=0.6){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+24]++;
				}else if(R1<=0.7){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+25]++;
				}else if(R1<=0.8){
					Colltype = 3; //"26.e+O2A>2e+O^+OP");
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+26]++;
				}else{
					Colltype = 0;
					Null++;
				}
				break;
			}
			case 4:{ // E + O2B
				mofm = info_CX[27].mofM;
                if(R1<=0.1){
					Colltype = 3; //"22.e+O2B>2e+O2^");
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 2;
					MCCR[TID*TnRct+27]++;
				}else if(R1<=0.2){
					Colltype = 2; //"23.e+O2B>OP+O-");
					MCCR[TID*TnRct+28]++;
				}else if(R1<=0.3){
					engy+=1.627f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+29]++;
				}else if(R1<=0.4){
					engy+=1.627f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+30]++;
				}else if(R1<=0.5){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+31]++;
				}else if(R1<=0.6){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+32]++;
				}else if(R1<=0.7){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+33]++;
				}else if(R1<=0.8){
					Colltype = 3; //"34.e+O2B>2e+O^+OP");
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+34]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 5:{ // E + O-
				mofm = info_CX[35].mofM;
                if(R1<=0.5){
					Colltype = 5; //"30.e+O->2e+OP");
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+35]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 6:{ // E + O2+
				mofm = info_CX[36].mofM;
                if(R1<=0.5){
					Colltype = 4; //"36.e+O2^>OP+OD");
					MCCR[TID*TnRct+36]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 7:{ // E + OP
				mofm = info_CX[37].mofM;
                if(R1<=0.125){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+37]++;
				}else if(R1<=0.25){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+38]++;
				}else if(R1<=0.375){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+39]++;
				}else if(R1<=0.5){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+40]++;
				}else if(R1<=0.625){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+41]++;
				}else if(R1<=0.75){
					Colltype = 3; //"42.e+OP>2e+O^");
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+42]++;
				}else if(R1<=0.875){
					engy/=2.0f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+43]++;
				}else{
					Colltype = 0;
					Null++;
				}
                break;
			}
			case 8:{ // E + OD
				mofm = info_CX[44].mofM;
                if(R1<=0.3333){
					Colltype = 3; //"44.e+OD>2e+O^");
					engy/=2.0f;
					rengy=10.0*__tanf(curand_uniform(&LocalStates)*atan(engy/20.0));
					engy-=rengy;
					vel = sqrt(fabs(rengy)/info[0].Escale);
					vel2 = sqrt(fabs(engy)/info[0].Escale);
					Iz_isp1 = 0;
					Iz_isp2 = 3;
					MCCR[TID*TnRct+44]++;
				}else if(R1<=0.6666){
					engy+=1.96f;
					vel=sqrt(fabs(engy)/info[0].Escale);
					MCCR[TID*TnRct+45]++;
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
			oldPNC = atomicAdd(&data[TID+4*Gsize].PtNumInCell,1);
			index = info[4].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID+4*Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			n = (nvel-1)*curand_uniform(&LocalStates);
			dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
            break;
        case 3: // 3 : ionization 1 -> 3 charged
			// second charged create
			/*
			oldPNC = atomicAdd(&data[TID+Iz_isp1*Gsize].PtNumInCell,1);
			index = info[Iz_isp1].St_num + TID + oldPNC*Gsize;
			sp[index].CellID = TID+Iz_isp1*Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			if(Iz_isp1 == 0){
				sp[index].vx = VX;
				sp[index].vy = VY;
				sp[index].vz = VZ;
				dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			}else{
				n = (nvel-1)*curand_uniform(&LocalStates);
				dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			}*/
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
				dev_anewvel(rengy,vel,&sp[index].vx,&sp[index].vy,&sp[index].vz,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			}else{
				n = (nvel-1)*curand_uniform(&LocalStates);
				dev_maxwellv(&sp[index].vx,&sp[index].vy,&sp[index].vz,vsave[n],BG[TID].BackVel1,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
			}
			// energy loss electron 
			dev_anewvel(engy,vel2,&VX,&VY,&VZ,0,mofm,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
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
			oldPNC = atomicAdd(&data[TID+4*Gsize].PtNumInCell,0);
			if(oldPNC>1){
				R2 = curand_uniform(&LocalStates);
				index = info[4].St_num + TID + oldPNC*Gsize;
				index2 = (int)((float)oldPNC * R2);
				index3 = info[4].St_num + TID + index2*Gsize;
				sp[index3].CellID = sp[index].CellID;
				sp[index3].x = sp[index].x;
				sp[index3].y = sp[index].y;
				sp[index3].vx = sp[index].vx;
				sp[index3].vy = sp[index].vy;
				sp[index3].vz = sp[index].vz;
				atomicAdd(&data[TID+4*Gsize].PtNumInCell,-1);
				//printf("2[%d][%d]: %g,%g,%g,%g,%g,\n",TID,sp[index].x,sp[index].y,sp[index].vx,sp[index].vy,sp[i].vz);
			}else if(oldPNC == 1){
				atomicAdd(&data[TID+4*Gsize].PtNumInCell,-1);
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
__device__ void ArO2_Ar_ion_TEST(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
	int i,j,k,n,index,index2,index3;
	int ID,PNMC,MPNC,Null,Flag;
	int Target,oldPNC;
	int Colltype;
	float mofm,R1,R2;
	float VX,VY,VZ;
	float engy,dum,vel;	
	float SumSigma,SumEngyLoss;
	float vneut,vneutx,vneuty,vneutz;

	ID = TID%Gsize;
    curandState LocalStates = states[TID];
	PNMC = data[TID].PtNumMCCInCell;
	MPNC = data[TID].MaxPtNumInCell;
	Null = 0;
    // Calculate total Collision probability
	i = info[1].St_num + ID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
        // Calculate energy
		Flag = sp[i].CellID;
		if(Flag == 0) vneut = BG[ID].BackVel1;
		else  vneut = BG[ID].BackVel2;
		n = (nvel-1)*curand_uniform(&LocalStates);
		dev_maxwellv(&vneutx,&vneuty,&vneutz,vsave[n],vneut,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
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
		R1 = curand_uniform(&LocalStates);
        switch(Flag){
			case 0:{
				mofm = info_CX[65].mofM;
				if(engy > info_CX[65].Th_e&& R1<=0.33333){
					Colltype = 2; 
					MCCR[ID*TnRct+65]++;
				}else if(engy > info_CX[66].Th_e && R1<=0.66666){
					Colltype = 1; 
					MCCR[ID*TnRct+66]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[67].mofM;
				if(engy > info_CX[67].Th_e && R1<=0.5){
					Colltype = 1; 
					MCCR[ID*TnRct+67]++;
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
        case 2: // 2 : Charge exchange 
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
	data[TID].PtNullMCCInCell = Null;
}
__device__ void ArO2_O2_ion_TEST(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
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
		// 4 : Charge exchange AR+
        switch(Flag){
			case 0:{
				mofm = info_CX[52].mofM;
                R1 = curand_uniform(&LocalStates);
				if(R1<=0.5){
					Colltype = 3; 
					MCCR[ID*TnRct+52]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[53].mofM;
                R1 = curand_uniform(&LocalStates);
				if(R1<=0.25){
					Colltype = 2; 
					MCCR[ID*TnRct+53]++;
				}else if(R1<=0.5){
					Colltype = 1; 
					MCCR[ID*TnRct+54]++;
				}else if(R1<=0.75){
					Colltype = 3; 
					MCCR[ID*TnRct+55]++;
				}else{
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[56].mofM;
                R1 = curand_uniform(&LocalStates);
                if(R1<=0.5){
					Colltype = 2; 
					MCCR[ID*TnRct+56]++;
				}else{
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[57].mofM;
                R1 = curand_uniform(&LocalStates);
                if(R1<=0.5){
					Colltype = 2; 
					MCCR[ID*TnRct+57]++;
				}else{
					Null++;
				}
                break;
			}
			case 4:{
				mofm = info_CX[58].mofM;
                R1 = curand_uniform(&LocalStates);
                if(R1<=0.33333){
					Colltype = 4; 
					MCCR[ID*TnRct+58]++;
				}else if(R1<=0.66666){
					Colltype = 2; 
					MCCR[ID*TnRct+59]++;
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
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[2].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = vneutx;
			sp[index].vy = vneuty;
			sp[index].vz = vneutz;
			break;
        case 3: // 3 : Charge exchange o+
			oldPNC = atomicAdd(&data[TID+Gsize].PtNumInCell,1);
			index = info[3].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID+Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = vneutx;
			sp[index].vy = vneuty;
			sp[index].vz = vneutz;
            break;
		case 4: // 4 : Charge exchange AR+
			oldPNC = atomicAdd(&data[TID-Gsize].PtNumInCell,1);
			index = info[1].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID-Gsize;
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
__device__ void ArO2_O_ion_TEST(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
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
	i = info[3].St_num + ID + (MPNC-1)*Gsize;
	for(k=0;k<PNMC;k++){
        // Calculate energy
		n = (nvel-1)*curand_uniform(&LocalStates);
		dev_maxwellv(&vneutx,&vneuty,&vneutz,vsave[n],BG[ID].BackVel2,curand_uniform(&LocalStates),curand_uniform(&LocalStates));
		VX = sp[i].vx;
		VY = sp[i].vy;
		VZ = sp[i].vz;
        Flag = sp[i].CellID;
        VX = sp[i].vx - vneutx;
		VY = sp[i].vy - vneuty;
		VZ = sp[i].vz - vneutz;
		dum = VX*VX+VY*VY+VZ*VZ;
		vel = sqrt(dum);
		engy = info[3].Escale * dum;
        Colltype = 0;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Scattering
        // 2 : Charge exchange O2+
		// 3 : Charge exchange O+
        switch(Flag){
			case 0:{
				mofm = info_CX[60].mofM;
                R1 = curand_uniform(&LocalStates);
				if(engy > info_CX[60].Th_e && R1<=0.33333){
					Colltype = 2; 
					MCCR[ID*TnRct+60]++;
				}else if(engy > info_CX[61].Th_e && R1<=0.66666){
					Colltype = 1; 
					MCCR[ID*TnRct+61]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[62].mofM;
                R1 = curand_uniform(&LocalStates);
				if(engy > info_CX[62].Th_e && R1<=0.5){
					Colltype = 3; 
					MCCR[ID*TnRct+62]++;
				}else{
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[63].mofM;
                R1 = curand_uniform(&LocalStates);
                if(engy > info_CX[63].Th_e && R1<=0.5){
					Colltype = 2; 
					MCCR[ID*TnRct+63]++;
				}else{
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[64].mofM;
                R1 = curand_uniform(&LocalStates);
                if(engy > info_CX[64].Th_e && R1<=0.5){
					Colltype = 2; 
					MCCR[ID*TnRct+64]++;
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
			index = info[3].St_num + ID + oldPNC*Gsize;
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
			index = info[3].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX+vneutx;
			sp[index].vy = VY+vneuty;
			sp[index].vz = VZ+vneutz;
            break;
        case 2: // 2 : Charge exchange o2+
			oldPNC = atomicAdd(&data[TID-Gsize].PtNumInCell,1);
			index = info[2].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID-Gsize;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = vneutx;
			sp[index].vy = vneuty;
			sp[index].vz = vneutz;
			break;
        case 3: // 3 : Charge exchange o+
			oldPNC = atomicAdd(&data[TID].PtNumInCell,1);
			index = info[3].St_num + ID + oldPNC*Gsize;
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
__device__ void ArO2_O_negative_TEST(int Gsize, int ngy, int TID, int nvel, float *vsave, curandState *states, 
											Species *info, GPG *data, GCP *sp, int N_LOGX, float idLOGX, 
											MCC_sigmav *sigv, CollF *info_CX, ArO2CollD *CX, int TnRct,float *MCCR, GGA *BG){
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
	i = info[4].St_num + ID + (MPNC-1)*Gsize;
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
		engy = info[4].Escale * dum;
        Colltype = 0;
        //Start
        // Colltype
        // 0 : Null collision
        // 1 : Scattering
        // 2 : Detachment using maxwellv
		// 3 : dissociative  recombination just delete
        switch(Flag){
			case 0:{
				mofm = info_CX[46].mofM;
                R1 = curand_uniform(&LocalStates);
				if(R1<=0.3333){
					Colltype = 1; 
					MCCR[ID*TnRct+46]++;
				}else if(R1<=0.6666){
					Colltype = 2; 
					MCCR[ID*TnRct+47]++;
				}else{
					Null++;
				}
				break;
			}
			case 1:{
				mofm = info_CX[48].mofM;
                R1 = curand_uniform(&LocalStates);
				if(R1<=0.5){
					Colltype = 2; 
					MCCR[ID*TnRct+48]++;
				}else{
					Null++;
				}
				break;
			}
			case 2:{
				mofm = info_CX[49].mofM;
                R1 = curand_uniform(&LocalStates);
                if(R1<=0.5){
					Colltype = 3; 
					MCCR[ID*TnRct+49]++;
				}else{
					Null++;
				}
                break;
			}
			case 3:{
				mofm = info_CX[50].mofM;
                R1 = curand_uniform(&LocalStates);
                if(R1<=0.5){
					Colltype = 3; 
					MCCR[ID*TnRct+50]++;
				}else{
					Null++;
				}
                break;
			}
			case 4:{
				mofm = info_CX[51].mofM;
                R1 = curand_uniform(&LocalStates);
                if(R1<=0.5){
					Colltype = 2; 
					MCCR[ID*TnRct+51]++;
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
			index = info[4].St_num + ID + oldPNC*Gsize;
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
			index = info[4].St_num + ID + oldPNC*Gsize;
			sp[index].CellID = TID;
			sp[index].x = sp[i].x;
			sp[index].y = sp[i].y;
			sp[index].vx = VX;
			sp[index].vy = VY;
			sp[index].vz = VZ;
            break;
        case 2: // 2 : Detachment using maxwellv
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
__device__ void  ArO2_Collision_Check_TEST(int Gsize, int Csize, int ngy, int TID, float dt, int MCCn, float dtm, float dx, float dy,
                                        curandState *states, Species *info, GPG *data, GCP *sp, MCC_sigmav *sigv, GGA *BG, GFC *Fluid){
	int i,j,k,index,Randn;
	int ID,isp,CID,PNMC,MPNC;
    int PNC,Flag;
	int nx,ny,ngx;
	float Tprob,Prob1,Prob2,Prob3,Prob4,Prob5,Prob6,Prob7,Prob8,Prob9;
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
		Prob1 = 0.01;
		Prob2 = Prob1 + 0.01;
        Prob3 = Prob2 + 0.01;
	    Prob4 = Prob3 + 0.01;
        Prob5 = Prob4 + 0.01;
        Prob6 = Prob5 + 0.01;
        Prob7 = Prob6 + 0.01;
        Prob8 = Prob7 + 0.01;
        Prob9 = Prob8 + 0.01;
	    Tprob = Prob9; 
		Randn = MCCn;
        break;
	case 1: // Ar+
		Prob1 = 0.1;
	    Prob2 = Prob1 + 0.1;
		Tprob = Prob2;
		Randn = MCCn;
		break;
    case 2: // O2+
        Prob1 = 0.01;
	    Prob2 = Prob1 + 0.01;
	    Prob3 = Prob2 + 0.01;
	    Prob4 = Prob3 + 0.01;
		Prob5 = Prob4 + 0.01;
        Tprob = Prob5;
		Randn = MCCn;
        break;
    case 3: // O+
        Prob1 = 0.01;
	    Prob2 = Prob1 + 0.01;
	    Prob3 = Prob2 + 0.01;
	    Prob4 = Prob3 + 0.01;
	    Tprob = Prob4;
		Randn = MCCn;
        break;
    case 4: // O-
        Prob1 = 0.01;
	    Prob2 = Prob1 + 0.01;
	    Prob3 = Prob2 + 0.01;
	    Prob4 = Prob3 + 0.01;
	    Prob5 = Prob4 + 0.01;
	    Tprob = Prob5;
		Randn = MCCn;
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
			//printf("k[%d], PNMC[%d], R1[%g], Tprob[%g]\n",k,PNMC,R1,Tprob);
            switch (isp){
            case 0:
                if(R1 <= Prob1)	        Flag = (int)0;
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
                else if(R1 <= Prob4)	Flag = (int)3;
                else if(R1 <= Prob5)	Flag = (int)4;
                else if(R1 <= Prob6)	Flag = (int)5;
				else if(R1 <= Prob7)	Flag = (int)6;
				else if(R1 <= Prob8)	Flag = (int)7;
		        else			        Flag = (int)8;
				//printf("k[%d], PNMC[%d], R1[%g], Tprob[%g], Flag[%d]\n",k,PNMC,R1,Tprob,Flag);
                break;
            case 1:
                if(R1 <= Prob1)	        Flag = (int)0;
		        else			        Flag = (int)1;
                break;
            case 2:
                if(R1 <= Prob1)	        Flag = (int)0;
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
				else if(R1 <= Prob4)	Flag = (int)3;
		        else			        Flag = (int)4;
                break;
            case 3:
                if(R1 <= Prob1)	        Flag = (int)0;
                else if(R1 <= Prob2)	Flag = (int)1;
                else if(R1 <= Prob3)	Flag = (int)2;
		        else			        Flag = (int)3;
                break;
            case 4:
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
__device__ float ArO2_CrossSection(int R, float engy, int N_LOGX, float idLOGX, ArO2CollD *data){
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
        case 58 :
            if(lengy < data[0].xe){
			    return data[0].cx_58;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_58 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_58+a1*data[ee2+1].cx_58;
            break;
        case 59 :
            if(lengy < data[0].xe){
			    return data[0].cx_59;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_59 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_59+a1*data[ee2+1].cx_59;
            break;
        case 60 :
            if(lengy < data[0].xe){
			    return data[0].cx_60;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_60 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_60+a1*data[ee2+1].cx_60;
            break;
        case 61 :
            if(lengy < data[0].xe){
			    return data[0].cx_61;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_61 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_61+a1*data[ee2+1].cx_61;
            break;
        case 62 :
            if(lengy < data[0].xe){
			    return data[0].cx_62;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_62 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_62+a1*data[ee2+1].cx_62;
            break;
        case 63 :
            if(lengy < data[0].xe){
			    return data[0].cx_63;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_63 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_63+a1*data[ee2+1].cx_63;
            break;
        case 64 :
            if(lengy < data[0].xe){
			    return data[0].cx_64;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_64 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_64+a1*data[ee2+1].cx_64;
            break;
        case 65 :
            if(lengy < data[0].xe){
			    return data[0].cx_65;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_65 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_65+a1*data[ee2+1].cx_65;
            break;
        case 66 :
            if(lengy < data[0].xe){
			    return data[0].cx_66;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_66 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_66+a1*data[ee2+1].cx_66;
            break;
        case 67 :
            if(lengy < data[0].xe){
			    return data[0].cx_67;
		    }else if(lengy > data[N_LOGX-1].xe){
			    return data[N_LOGX-1].cx_67 * 0.1 * exp(-1 * (lengy - data[N_LOGX-1].xe));
		    }
		    return a2*data[ee2].cx_67+a1*data[ee2+1].cx_67;
            break;
        default :
            printf("\nError : Call about cross section data in ArO2MCC.\n\n");
            return 0.0;
    }
}