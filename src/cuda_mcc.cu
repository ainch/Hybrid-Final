#include "cuda_mcc.cuh"
void MCC_Ar_cuda(){
	//MCC_Ar_Basic<<<MCC_GRID2, MCC_BLOCK2>>>(Gsize, ngy, devStates, dev_info_sp, dev_sp, dev_G_sp, dev_GvecSet);
    void *kernelArgs[] = {
        (void*)&Gsize,(void*)&ngy,(void*)&DT_MCCn,(void*)&dt_mcc,
		(void*)&idx,(void*)&idy,(void*)&h_nvel,(void*)&dev_vsave,
		(void*)&devStates, (void*)&dev_SigmaV,(void*)&dev_Coll_Flag,(void*)&dev_ArCX,
		(void*)&dev_FG, (void*)&dev_C_F,(void*)&dev_GvecSet, 
		(void*)&dev_info_sp,  (void*)&dev_G_sp, (void*)&dev_sp,
    };
    cudaLaunchCooperativeKernel((void *)MCC_Ar_cooper,MCC_GRID,MCC_BLOCK,kernelArgs, sMemSize_MCC, NULL);
    cudaDeviceSynchronize();
	exit(1);
}
__global__ void MCC_Ar_cooper(int Gsize, int ngy, float dt, int MCCn, float dtm,float idx,float idy, int nvel, float *vsave,
											curandState *states, MCC_sigmav *sigv, CollF *CollP, ArCollD *CX, 
											Fluid *infoF, GFC *Fluid, GGA *BG, Species *info, GPG *data, GCP *sp){
	cg::thread_block cta = cg::this_thread_block();
    cg::grid_group grid = cg::this_grid();

	for (int i=grid.thread_rank(); i < 2*Gsize; i+= grid.size()){
		printf("I = %d\n",i);
	}
}
void Set_NullCollisionTime_cuda(){
    // This function calculates the following variables :
    // 1. int DT_MCCn
    // 2. float dt_mcc
    // 3. dev_Coll_Flag [CollF] - [TnRct]
    // 4. dev_**CX [**CollD] - [N_LOGX]
    // 5. Max_sigma_v[num_a] ,Ttarget[num_a] ColProb[num_b]
    //    Argon - [3][2]
    //    Oxygen - [20][4]
    //    ArO2 - [25][5]
    int i;
    int num_a; // Number of (Projectile + target)
    float engy;
	float *sigma;
    float *Ttarget,*ColProb;
    float ratio, prob_cut = 0.05; // = 5%. Collision rate per time step
    //Cross section Data copy  CPU >> GPU
    checkCudaErrors(cudaMalloc((void**)&dev_Coll_Flag, TnRct * sizeof(CollF)));
    checkCudaErrors(cudaMemcpy(dev_Coll_Flag, Coll_Flag, TnRct * sizeof(CollF), cudaMemcpyHostToDevice));
    if(MainGas == ARGON){
        checkCudaErrors(cudaMalloc((void**)&dev_ArCX, N_LOGX * sizeof(ArCollD)));
        checkCudaErrors(cudaMemcpy(dev_ArCX, Ar_Data, N_LOGX * sizeof(ArCollD), cudaMemcpyHostToDevice));
        num_a = 3; 			
        Host_SigmaV = (MCC_sigmav *)malloc(num_a * sizeof(MCC_sigmav));
        for(i=0;i<num_a;i++) Host_SigmaV[i].val = 0.0;
		sigma = VFMalloc(num_a);
		Ttarget    = VFMalloc(num_a);
		VFInit(Ttarget, 0.0, num_a);
		ColProb    = VFMalloc(nsp);
		VFInit(ColProb, 0.0, nsp);
		for(i=0;i<N_LOGX;i++){
			// Electron
			engy = Ar_Data[i].xee;
			sigma[0] = Ar_Data[i].cx_0+Ar_Data[i].cx_1+Ar_Data[i].cx_2+Ar_Data[i].cx_3; // e + Ar
			Host_SigmaV[0].val = max(Host_SigmaV[0].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[0]);
			sigma[1] = Ar_Data[i].cx_4;    // e + Ar*
			Host_SigmaV[1].val = max(Host_SigmaV[1].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[1]);
			// Ar ion
			engy = Ar_Data[i].xee;
			sigma[2] = Ar_Data[i].cx_5+Ar_Data[i].cx_6;
			Host_SigmaV[2].val = max(Host_SigmaV[2].val,sqrt(2*1.602e-19*engy/SP[1].mass)*sigma[2]);
		}
		Ttarget[0] = 1.0-exp(-1 * dt * Host_SigmaV[0].val * BG[0].InitDens);
		Ttarget[1] = 1.0-exp(-1 * dt * Host_SigmaV[1].val * FG[0].InitDens);
		Ttarget[2] = 1.0-exp(-1 * dt * Host_SigmaV[2].val * BG[0].InitDens);
		ColProb[0] = Ttarget[0] + Ttarget[1];
		ColProb[1] = Ttarget[2];
		if(ColProb[0] > prob_cut) {
			DT_MCCn = (int)ceil(ColProb[0]/prob_cut);
			ratio = 1/(float)DT_MCCn;
			dt_mcc = dt*ratio;
		}
		else {
			DT_MCCn = 1;
			dt_mcc = dt;
		}
		// Null Method Information
		fprintf(stderr,"--------------<Null Collision Information>-------------\n");
		fprintf(stderr, " - Total Collision probability per Time step\n");
		fprintf(stderr, "   Electon - %2.4f %\n",ColProb[0] * 100);
		fprintf(stderr, "   Ar ion  - %2.4f %\n",ColProb[1] * 100);
		fprintf(stderr, " - Number of Electron MCC Cycle\n");
		fprintf(stderr, "   Cycle : %d, dt_mcc : %g \n",DT_MCCn, dt_mcc);
		fprintf(stderr,"-------------------------------------------------------\n");
    }else if(MainGas == OXYGEN){
        checkCudaErrors(cudaMalloc((void**)&dev_O2CX, N_LOGX * sizeof(O2CollD)));
        checkCudaErrors(cudaMemcpy(dev_O2CX, O2_Data, N_LOGX * sizeof(O2CollD), cudaMemcpyHostToDevice));
		num_a = 20; 			
		Host_SigmaV = (MCC_sigmav *)malloc(num_a * sizeof(MCC_sigmav));
        for(i=0;i<num_a;i++) Host_SigmaV[i].val = 0.0;
		sigma = VFMalloc(num_a);
		Ttarget    = VFMalloc(num_a);
		VFInit(Ttarget, 0.0, num_a);
		ColProb    = VFMalloc(nsp);
		VFInit(ColProb, 0.0, nsp);
		for(i=0;i<N_LOGX;i++){
			engy = O2_Data[i].xee;
			sigma[0] = O2_Data[i].cx_0+O2_Data[i].cx_1+O2_Data[i].cx_2+O2_Data[i].cx_3+O2_Data[i].cx_4;
			sigma[0] += O2_Data[i].cx_5+O2_Data[i].cx_6+O2_Data[i].cx_7+O2_Data[i].cx_8+O2_Data[i].cx_9;
			sigma[0] += O2_Data[i].cx_10+O2_Data[i].cx_11+O2_Data[i].cx_12+O2_Data[i].cx_13;
			Host_SigmaV[0].val = max(Host_SigmaV[0].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[0]);
			sigma[1] = O2_Data[i].cx_14+O2_Data[i].cx_15+O2_Data[i].cx_16+O2_Data[i].cx_17+O2_Data[i].cx_18;
			sigma[1] += O2_Data[i].cx_19+O2_Data[i].cx_20+O2_Data[i].cx_21;
			Host_SigmaV[1].val = max(Host_SigmaV[1].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[1]);
			sigma[2] = O2_Data[i].cx_22+O2_Data[i].cx_23+O2_Data[i].cx_24+O2_Data[i].cx_25+O2_Data[i].cx_26;
			sigma[2] += O2_Data[i].cx_27+O2_Data[i].cx_28+O2_Data[i].cx_29;
			Host_SigmaV[2].val = max(Host_SigmaV[2].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[2]);
			sigma[3] = O2_Data[i].cx_30;
			Host_SigmaV[3].val = max(Host_SigmaV[3].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[3]);
			sigma[4] = O2_Data[i].cx_31;
			Host_SigmaV[4].val = max(Host_SigmaV[4].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[4]);
			sigma[5] = O2_Data[i].cx_32+O2_Data[i].cx_33+O2_Data[i].cx_34+O2_Data[i].cx_35+O2_Data[i].cx_36;
			sigma[5] += O2_Data[i].cx_37+O2_Data[i].cx_38;
			Host_SigmaV[5].val = max(Host_SigmaV[5].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[5]);
			sigma[6] = O2_Data[i].cx_39+O2_Data[i].cx_40;
			Host_SigmaV[6].val = max(Host_SigmaV[6].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[6]);
			sigma[7] = O2_Data[i].cx_41 + O2_Data[i].cx_42;
			Host_SigmaV[7].val = max(Host_SigmaV[7].val,sqrt(2*1.602e-19*engy/SP[3].mass)*sigma[7]);
			sigma[8] = O2_Data[i].cx_43;
			Host_SigmaV[8].val = max(Host_SigmaV[8].val,sqrt(2*1.602e-19*engy/SP[3].mass)*sigma[8]);
			sigma[9] = O2_Data[i].cx_44;
			Host_SigmaV[9].val = max(Host_SigmaV[9].val,sqrt(2*1.602e-19*engy/SP[3].mass)*sigma[9]);
			sigma[10] = O2_Data[i].cx_45;
			Host_SigmaV[10].val = max(Host_SigmaV[10].val,sqrt(2*1.602e-19*engy/SP[3].mass)*sigma[10]);
			sigma[11] = O2_Data[i].cx_46;
			Host_SigmaV[11].val = max(Host_SigmaV[11].val,sqrt(2*1.602e-19*engy/SP[3].mass)*sigma[11]);
			sigma[12] = O2_Data[i].cx_47;
			Host_SigmaV[12].val = max(Host_SigmaV[12].val,sqrt(2*1.602e-19*engy/SP[1].mass)*sigma[12]);
			sigma[13] = O2_Data[i].cx_48 + O2_Data[i].cx_49 + O2_Data[i].cx_50;
			Host_SigmaV[13].val = max(Host_SigmaV[13].val,sqrt(2*1.602e-19*engy/SP[1].mass)*sigma[13]);
			sigma[14] = O2_Data[i].cx_51;
			Host_SigmaV[14].val = max(Host_SigmaV[14].val,sqrt(2*1.602e-19*engy/SP[1].mass)*sigma[14]);
			sigma[15] = O2_Data[i].cx_52;
			Host_SigmaV[15].val = max(Host_SigmaV[15].val,sqrt(2*1.602e-19*engy/SP[1].mass)*sigma[15]);
			sigma[16] = O2_Data[i].cx_53 + O2_Data[i].cx_54;
			Host_SigmaV[16].val = max(Host_SigmaV[16].val,sqrt(2*1.602e-19*engy/SP[2].mass)*sigma[16]);
			sigma[17] = O2_Data[i].cx_55;
			Host_SigmaV[17].val = max(Host_SigmaV[17].val,sqrt(2*1.602e-19*engy/SP[2].mass)*sigma[17]);
			sigma[18] = O2_Data[i].cx_56;
			Host_SigmaV[18].val = max(Host_SigmaV[18].val,sqrt(2*1.602e-19*engy/SP[2].mass)*sigma[18]);
			sigma[19] = O2_Data[i].cx_57;
			Host_SigmaV[19].val = max(Host_SigmaV[19].val,sqrt(2*1.602e-19*engy/SP[2].mass)*sigma[19]);
		}
		Ttarget[0] = 1.0-exp(-1 * dt * Host_SigmaV[0].val * BG[0].InitDens);// E + O2
		Ttarget[1] = 1.0-exp(-1 * dt * Host_SigmaV[1].val * FG[0].InitDens);// E + O2A
		Ttarget[2] = 1.0-exp(-1 * dt * Host_SigmaV[2].val * FG[1].InitDens);// E + O2B
		Ttarget[3] = 1.0-exp(-1 * dt * Host_SigmaV[3].val * SP[3].InitDens);// E + O-
		Ttarget[4] = 1.0-exp(-1 * dt * Host_SigmaV[4].val * SP[1].InitDens);// E + O2^
		Ttarget[5] = 1.0-exp(-1 * dt * Host_SigmaV[5].val * FG[2].InitDens);// E + OP
		Ttarget[6] = 1.0-exp(-1 * dt * Host_SigmaV[6].val * FG[3].InitDens);// E + OD
		Ttarget[7] = 1.0-exp(-1 * dt * Host_SigmaV[7].val * BG[0].InitDens);// O- + O2
		Ttarget[8] = 1.0-exp(-1 * dt * Host_SigmaV[8].val * FG[2].InitDens);// O- + OP
		Ttarget[9] = 1.0-exp(-1 * dt * Host_SigmaV[9].val* SP[1].InitDens);// O- + O2^
		Ttarget[10] = 1.0-exp(-1 * dt * Host_SigmaV[10].val * SP[2].InitDens);// O- + O^
		Ttarget[11] = 1.0-exp(-1 * dt * Host_SigmaV[11].val * FG[0].InitDens);// O- + O2A
		Ttarget[12] = 1.0-exp(-1 * dt * Host_SigmaV[12].val * FG[2].InitDens);// O2^ + OP
		Ttarget[13] = 1.0-exp(-1 * dt * Host_SigmaV[13].val * BG[0].InitDens);// O2^ + O2
		Ttarget[14] = 1.0-exp(-1 * dt * Host_SigmaV[14].val * FG[0].InitDens);// O2^ + O2A
		Ttarget[15] = 1.0-exp(-1 * dt * Host_SigmaV[15].val * FG[1].InitDens);// O2^ + O2B
		Ttarget[16] = 1.0-exp(-1 * dt * Host_SigmaV[16].val * BG[0].InitDens);// O^ + O2
		Ttarget[17] = 1.0-exp(-1 * dt * Host_SigmaV[17].val * FG[2].InitDens);// O^ + OP
		Ttarget[18] = 1.0-exp(-1 * dt * Host_SigmaV[18].val* FG[0].InitDens);// O^ + O2A
		Ttarget[19] = 1.0-exp(-1 * dt * Host_SigmaV[19].val * FG[1].InitDens);// O^ + O2B
		ColProb[0] = Ttarget[0] + Ttarget[1] + Ttarget[2] + Ttarget[3] + Ttarget[4] + Ttarget[5] + Ttarget[6];// Electron
		ColProb[1] = Ttarget[7] + Ttarget[8] + Ttarget[9] + Ttarget[10] + Ttarget[11];// O-
		ColProb[2] = Ttarget[12] + Ttarget[13] + Ttarget[14] + Ttarget[15];// O2+
		ColProb[3] = Ttarget[16] + Ttarget[17] + Ttarget[18] + Ttarget[19];// O+
		if(ColProb[0] > prob_cut) {
			DT_MCCn = (int)ceil(ColProb[0]/prob_cut);
			ratio = 1/(float)DT_MCCn;
			dt_mcc = dt*ratio;
		}
		else {
			DT_MCCn = 1;
			dt_mcc = dt;
		}
		// Null Method Information
		fprintf(stderr,"--------------<Null Collision Information>-------------\n");
		fprintf(stderr, " - Total Collision probability per Time step\n");
		fprintf(stderr, "   Electon - %2.4f %\n",ColProb[0] * 100);
		fprintf(stderr, "   O- ion  - %2.4f %\n",ColProb[1] * 100);
		fprintf(stderr, "   O2+ion  - %2.4f %\n",ColProb[2] * 100);
		fprintf(stderr, "   O+ ion  - %2.4f %\n",ColProb[3] * 100);
		fprintf(stderr, " - Number of Electron MCC Cycle\n");
		fprintf(stderr, "   Cycle : %d, dt_mcc : %g \n",DT_MCCn, dt_mcc);
		fprintf(stderr,"-------------------------------------------------------\n");
    }else if(MainGas == ARO2){
        checkCudaErrors(cudaMalloc((void**)&dev_ArO2CX, N_LOGX * sizeof(ArO2CollD)));
        checkCudaErrors(cudaMemcpy(dev_ArO2CX, ArO2_Data, N_LOGX * sizeof(ArO2CollD), cudaMemcpyHostToDevice));
		num_a = 25; 			
		Host_SigmaV = (MCC_sigmav *)malloc(num_a * sizeof(MCC_sigmav));
        for(i=0;i<num_a;i++) Host_SigmaV[i].val = 0.0;
		sigma = VFMalloc(num_a);
		Ttarget    = VFMalloc(num_a);
		VFInit(Ttarget, 0.0, num_a);
		ColProb    = VFMalloc(nsp);
		VFInit(ColProb, 0.0, nsp);
		for(i=0;i<N_LOGX;i++){
			engy = ArO2_Data[i].xee;
			sigma[0] = ArO2_Data[i].cx_0+ArO2_Data[i].cx_1+ArO2_Data[i].cx_2+ArO2_Data[i].cx_3; // e + Ar
			Host_SigmaV[0].val = max(Host_SigmaV[0].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[0]);
			sigma[1] = ArO2_Data[i].cx_4;    // e + Ar*
			Host_SigmaV[1].val = max(Host_SigmaV[1].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[1]);
			sigma[2] = ArO2_Data[i].cx_5+ArO2_Data[i].cx_6+ArO2_Data[i].cx_7+ArO2_Data[i].cx_8+ArO2_Data[i].cx_9;
			sigma[2] += ArO2_Data[i].cx_10+ArO2_Data[i].cx_11+ArO2_Data[i].cx_12+ArO2_Data[i].cx_13+ArO2_Data[i].cx_14;
			sigma[2] += ArO2_Data[i].cx_15+ArO2_Data[i].cx_16+ArO2_Data[i].cx_17+ArO2_Data[i].cx_18;
			Host_SigmaV[2].val = max(Host_SigmaV[2].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[2]);// e + O2
			sigma[3] = ArO2_Data[i].cx_19+ArO2_Data[i].cx_20+ArO2_Data[i].cx_21+ArO2_Data[i].cx_22+ArO2_Data[i].cx_23;
			sigma[3] += ArO2_Data[i].cx_24+ArO2_Data[i].cx_25+ArO2_Data[i].cx_26;
			Host_SigmaV[3].val = max(Host_SigmaV[3].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[3]);// e + O2A
			sigma[4] = ArO2_Data[i].cx_27+ArO2_Data[i].cx_28+ArO2_Data[i].cx_29+ArO2_Data[i].cx_30+ArO2_Data[i].cx_31;
			sigma[4] += ArO2_Data[i].cx_32+ArO2_Data[i].cx_33+ArO2_Data[i].cx_34;
			Host_SigmaV[4].val = max(Host_SigmaV[4].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[4]);// e + O2B
			sigma[5] = ArO2_Data[i].cx_35;
			Host_SigmaV[5].val = max(Host_SigmaV[5].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[5]);
			sigma[6] = ArO2_Data[i].cx_36;
			Host_SigmaV[6].val = max(Host_SigmaV[6].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[6]);
			sigma[7] = ArO2_Data[i].cx_37+ArO2_Data[i].cx_38+ArO2_Data[i].cx_39+ArO2_Data[i].cx_40+ArO2_Data[i].cx_41;
			sigma[7] += ArO2_Data[i].cx_42+ArO2_Data[i].cx_43;
			Host_SigmaV[7].val = max(Host_SigmaV[7].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[7]);
			sigma[8] = ArO2_Data[i].cx_44+ArO2_Data[i].cx_45;
			Host_SigmaV[8].val = max(Host_SigmaV[8].val,sqrt(2*1.602e-19*engy/SP[0].mass)*sigma[8]);
			sigma[9] = ArO2_Data[i].cx_46 + ArO2_Data[i].cx_47;
			Host_SigmaV[9].val = max(Host_SigmaV[9].val,sqrt(2*1.602e-19*engy/SP[4].mass)*sigma[9]);
			sigma[10] = ArO2_Data[i].cx_48;
			Host_SigmaV[10].val = max(Host_SigmaV[10].val,sqrt(2*1.602e-19*engy/SP[4].mass)*sigma[10]);
			sigma[11] = ArO2_Data[i].cx_49;
			Host_SigmaV[11].val = max(Host_SigmaV[11].val,sqrt(2*1.602e-19*engy/SP[4].mass)*sigma[11]);
			sigma[12] = ArO2_Data[i].cx_50;
			Host_SigmaV[12].val = max(Host_SigmaV[12].val,sqrt(2*1.602e-19*engy/SP[4].mass)*sigma[12]);
			sigma[13] = ArO2_Data[i].cx_51;
			Host_SigmaV[13].val = max(Host_SigmaV[13].val,sqrt(2*1.602e-19*engy/SP[4].mass)*sigma[13]);
			sigma[14] = ArO2_Data[i].cx_52;
			Host_SigmaV[14].val = max(Host_SigmaV[14].val,sqrt(2*1.602e-19*engy/SP[2].mass)*sigma[14]);
			sigma[15] = ArO2_Data[i].cx_53 + ArO2_Data[i].cx_54 + ArO2_Data[i].cx_55;
			Host_SigmaV[15].val = max(Host_SigmaV[15].val,sqrt(2*1.602e-19*engy/SP[2].mass)*sigma[15]);
			sigma[16] = ArO2_Data[i].cx_56;
			Host_SigmaV[16].val = max(Host_SigmaV[16].val,sqrt(2*1.602e-19*engy/SP[2].mass)*sigma[16]);
			sigma[17] = ArO2_Data[i].cx_57;
			Host_SigmaV[17].val = max(Host_SigmaV[17].val,sqrt(2*1.602e-19*engy/SP[2].mass)*sigma[17]);
			sigma[18] = ArO2_Data[i].cx_58 + ArO2_Data[i].cx_59;
			Host_SigmaV[18].val = max(Host_SigmaV[18].val,sqrt(2*1.602e-19*engy/SP[2].mass)*sigma[18]);
			sigma[19] = ArO2_Data[i].cx_60 + ArO2_Data[i].cx_61;
			Host_SigmaV[19].val = max(Host_SigmaV[19].val,sqrt(2*1.602e-19*engy/SP[3].mass)*sigma[19]);
			sigma[20] = ArO2_Data[i].cx_62;
			Host_SigmaV[20].val = max(Host_SigmaV[20].val,sqrt(2*1.602e-19*engy/SP[3].mass)*sigma[20]);
			sigma[21] = ArO2_Data[i].cx_63;
			Host_SigmaV[21].val = max(Host_SigmaV[21].val,sqrt(2*1.602e-19*engy/SP[3].mass)*sigma[21]);
			sigma[22] = ArO2_Data[i].cx_64;
			Host_SigmaV[22].val = max(Host_SigmaV[22].val,sqrt(2*1.602e-19*engy/SP[3].mass)*sigma[22]);
			sigma[23] = ArO2_Data[i].cx_65+ArO2_Data[i].cx_66;
			Host_SigmaV[23].val = max(Host_SigmaV[23].val,sqrt(2*1.602e-19*engy/SP[1].mass)*sigma[23]);
			sigma[24] = ArO2_Data[i].cx_67;
			Host_SigmaV[24].val = max(Host_SigmaV[24].val,sqrt(2*1.602e-19*engy/SP[1].mass)*sigma[24]);
		}
		Ttarget[0] = 1.0-exp(-1 * dt * Host_SigmaV[0].val * BG[0].InitDens);// E + Ar
		Ttarget[1] = 1.0-exp(-1 * dt * Host_SigmaV[1].val * FG[0].InitDens);// E + Ar*
		Ttarget[2] = 1.0-exp(-1 * dt * Host_SigmaV[2].val * BG[1].InitDens);// E + O2
		Ttarget[3] = 1.0-exp(-1 * dt * Host_SigmaV[3].val * FG[1].InitDens);// E + O2A
		Ttarget[4] = 1.0-exp(-1 * dt * Host_SigmaV[4].val * FG[2].InitDens);// E + O2B
		Ttarget[5] = 1.0-exp(-1 * dt * Host_SigmaV[5].val * SP[4].InitDens);// E + O-
		Ttarget[6] = 1.0-exp(-1 * dt * Host_SigmaV[6].val * SP[2].InitDens);// E + O2^
		Ttarget[7] = 1.0-exp(-1 * dt * Host_SigmaV[7].val * FG[3].InitDens);// E + OP
		Ttarget[8] = 1.0-exp(-1 * dt * Host_SigmaV[8].val * FG[4].InitDens);// E + OD
		Ttarget[9] = 1.0-exp(-1 * dt * Host_SigmaV[9].val * BG[1].InitDens);// O- + O2
		Ttarget[10] = 1.0-exp(-1 * dt * Host_SigmaV[10].val * FG[3].InitDens);// O- + OP
		Ttarget[11] = 1.0-exp(-1 * dt * Host_SigmaV[11].val* SP[2].InitDens);// O- + O2^
		Ttarget[12] = 1.0-exp(-1 * dt * Host_SigmaV[12].val * SP[3].InitDens);// O- + O^
		Ttarget[13] = 1.0-exp(-1 * dt * Host_SigmaV[13].val * FG[1].InitDens);// O- + O2A
		Ttarget[14] = 1.0-exp(-1 * dt * Host_SigmaV[14].val * FG[3].InitDens);// O2^ + OP
		Ttarget[15] = 1.0-exp(-1 * dt * Host_SigmaV[15].val * BG[1].InitDens);// O2^ + O2
		Ttarget[16] = 1.0-exp(-1 * dt * Host_SigmaV[16].val * FG[1].InitDens);// O2^ + O2A
		Ttarget[17] = 1.0-exp(-1 * dt * Host_SigmaV[17].val * FG[2].InitDens);// O2^ + O2B
		Ttarget[18] = 1.0-exp(-1 * dt * Host_SigmaV[18].val * BG[0].InitDens);// O2^ + Ar
		Ttarget[19] = 1.0-exp(-1 * dt * Host_SigmaV[19].val * BG[1].InitDens);// O^ + O2
		Ttarget[20] = 1.0-exp(-1 * dt * Host_SigmaV[20].val * FG[3].InitDens);// O^ + OP
		Ttarget[21] = 1.0-exp(-1 * dt * Host_SigmaV[21].val* FG[1].InitDens);// O^ + O2A
		Ttarget[22] = 1.0-exp(-1 * dt * Host_SigmaV[22].val * FG[2].InitDens);// O^ + O2B
		Ttarget[23] = 1.0-exp(-1 * dt * Host_SigmaV[23].val * BG[0].InitDens);// Ar^ + Ar
		Ttarget[24] = 1.0-exp(-1 * dt * Host_SigmaV[24].val * BG[1].InitDens);// Ar^ + O2
		ColProb[0] = Ttarget[0] + Ttarget[1] + Ttarget[2] + Ttarget[3] + Ttarget[4] + Ttarget[5] + Ttarget[6] + Ttarget[7] + Ttarget[8];// Electron
		ColProb[1] = Ttarget[23] + Ttarget[24];// Ar+
		ColProb[2] = Ttarget[14] + Ttarget[15] + Ttarget[16] + Ttarget[17] + Ttarget[18];// O2+
		ColProb[3] = Ttarget[19] + Ttarget[20] + Ttarget[21] + Ttarget[22];// O+
		ColProb[4] = Ttarget[9] + Ttarget[10] + Ttarget[11] + Ttarget[12] + Ttarget[13];// O-
		if(ColProb[0] > prob_cut) {
			DT_MCCn = (int)ceil(ColProb[0]/prob_cut);
			ratio = 1/(float)DT_MCCn;
			dt_mcc = dt*ratio;
		}
		else {
			DT_MCCn = 1;
			dt_mcc = dt;
		}
		// Null Method Information
		fprintf(stderr,"--------------<Null Collision Information>-------------\n");
		fprintf(stderr, " - Total Collision probability per Time step\n");
		fprintf(stderr, "   Electon - %2.4f %\n",ColProb[0] * 100);
		fprintf(stderr, "   Ar+ion  - %2.4f %\n",ColProb[1] * 100);
		fprintf(stderr, "   O2+ion  - %2.4f %\n",ColProb[2] * 100);
		fprintf(stderr, "   O+ ion  - %2.4f %\n",ColProb[3] * 100);
		fprintf(stderr, "   O- ion  - %2.4f %\n",ColProb[4] * 100);
		fprintf(stderr, " - Number of Electron MCC Cycle\n");
		fprintf(stderr, "   Cycle : %d, dt_mcc : %g \n",DT_MCCn, dt_mcc);
		fprintf(stderr,"-------------------------------------------------------\n");
    }else{
        printf("Error : MainGas = %d\n",MainGas);
        exit(1);
    }
    //Maximum sigma*v copy  CPU >> GPU
    checkCudaErrors(cudaMalloc((void**)&dev_SigmaV, num_a * sizeof(MCC_sigmav)));
    checkCudaErrors(cudaMemcpy(dev_SigmaV, Host_SigmaV, num_a * sizeof(MCC_sigmav), cudaMemcpyHostToDevice));
    printf("MCC Initializing Complete!\n");
}

__global__ void MCC_Ar_Basic(int Gsize, int ngy, curandState *states, Species *info, GCP *sp, GPG *data, GGA *Field){
	int TID = threadIdx.x + blockIdx.x * blockDim.x;
	int isp,ID;
 	isp = (int)TID/Gsize; //species number [< nsp]
    ID = (int)TID%Gsize; // Grid ID [< Gsize]
    if(TID>Gsize*info[isp].spnum) return;
	curandState LocalStates;
	LocalStates = states[TID];
	if(TID < 5 ){
		float a,b;
		a = curand_uniform(&LocalStates);
		b = curand_uniform(&LocalStates);
		printf("[%d] 1 = %g, 2 = %g\n",TID,a,b);
	}
	
	states[TID]=LocalStates;
}	