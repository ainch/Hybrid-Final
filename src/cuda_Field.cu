#include "cuda_Field.cuh"
void Efield_cuda_Basic(){ // No circuit
    // Function
    int i,j,isp;
    float K_t;
    float q_conv;   
	cudaMemcpy(CondCharge, dev_CondCharge, nsp * CondNUMR * sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemset((void *) dev_CondCharge, 0, nsp * CondNUMR * sizeof(float));
    // Get voltage 
    for(i=0;i<CondNUMR;i++) {
		phi_cond[i] = 0.0f;
		if(Efield_Flag[i]){
			for(j=0;j<SrcNUM;j++){
                if(SrcM_ID[j] == i+1){
                    phi_cond[i] += SrcDC[j] + SrcAC[j]*sin(Src2piFREQ[j]*t+SrcRPHASE[j]);
					V_t[i] = phi_cond[i];
                }
            }
		}
    }
    // Potential summation
	cudaMemset((void *) TotPotential, 0.0, Gsize * sizeof(float));
    for(i=0;i<CondNUMR;i++){
        LoadAT2D<<<FIELD_GRID2,FIELD_BLOCK2>>>(Lap_PHI_Sol, pitch, i, dev_phi_buf, Gsize);
		VectorSum<<<FIELD_GRID2, FIELD_BLOCK2>>>(Gsize, TotPotential, phi_cond[i], dev_phi_buf);
		//cudaDeviceSynchronize();
    }
	// TotPotential = laplace * votage
    cudaMemcpy(LapPotential, TotPotential, Gsize * sizeof(float),cudaMemcpyDeviceToDevice);
	// TotPotential = laplace * votage 
	// LapPotential = laplace * votage 
	VectorSum<<<FIELD_GRID2, FIELD_BLOCK2>>>(Gsize, TotPotential, 1, dev_phi);
	// TotPotential = laplace * votage + Poisson's
	// LapPotential = laplace * votage 
    GGACopy_Potential<<<FIELD_GRID2, FIELD_BLOCK2>>>(Gsize, dev_GvecSet, LapPotential, dev_phi);
	Cond_Sigma<<<FIELD_GRID2,FIELD_BLOCK2>>>(ngx, ngy, hdx, hdy, idxy, zlength, dev_GvecSet, dev_CvecSet, dev_info_sp, dev_G_sp, dev_phi, dev_Sigma);
	//Cond_Sigma_v2<<<FIELD_GRID2,FIELD_BLOCK2>>>(ngx, ngy, hdx, hdy, idxy, zlength, dev_GvecSet, dev_CvecSet, dev_info_sp, dev_G_sp, dev_phi, dev_Sigma);
	VtoEfield<<<FIELD_GRID2, FIELD_BLOCK2>>>(ngx,ngy,dx,dy,hdx,hdy,idx,idy, dev_Sigma, TotPotential, dev_G_sp, dev_CvecSet, dev_GvecSet);
    //cudaMemcpy(vec_G, dev_GvecSet, Gsize * sizeof(GGA), cudaMemcpyDeviceToHost); // for TEST
    //Main_Variable_printorSave(); // for TEST
	cudaDeviceSynchronize();
}
void Efield_cuda(){ // External circuit
    // Function
    int i,j,isp;
    float K_t;
    float q_conv;   
	cudaMemcpy(CondCharge, dev_CondCharge, nsp * CondNUMR * sizeof(float),cudaMemcpyDeviceToHost);
	cudaMemset((void *) dev_CondCharge, 0, nsp * CondNUMR * sizeof(float));
    // Get voltage 
    for(i=0;i<CondNUMR;i++) {
		V_t[i] = 0.0;
		if(Efield_Flag[i]){
			for(j=0;j<SrcNUM;j++){
                if(SrcM_ID[j] == i+1){
                    V_t[i] += SrcDC[j] + SrcAC[j]*sin(Src2piFREQ[j]*t+SrcRPHASE[j]);
                }
            }
		}
		K_t = (CC_a[i][1]*extq[i] + CC_a[i][2]*extq_1[i] + CC_a[i][3]*extq_2[i] + CC_a[i][4]*extq_3[i])/CC_a[i][0];
        //convective charge
		q_conv=0;
		for(isp=0;isp<nsp;isp++) {
			q_conv-=SP[isp].q_density*CondCharge[isp*CondNUMR+i];
		}
        b_t[i]=(V_t[i]-K_t)/CC_a[i][0] - extq[i] - q_conv + Surf_charge[i] - Pois_SIG_Sol[i];
		for(j=0;j<CondNUMR;j++){
			if(i==j) AM[i][j] = Lap_SIG_Sol[j][i]+1/CC_a[i][0];
			else AM[i][j] = Lap_SIG_Sol[j][i];
		}
    }
    cofactor(AM, CondNUMR);
    for (i = 0; i < CondNUMR; i++) {
		phi_cond[i] = 0.0;
		for (j = 0; j < CondNUMR; j++) {
			phi_cond[i] += AM[i][j] * b_t[j];
		}
        extq_3[i] = extq_2[i];
		extq_2[i] = extq_1[i];
		extq_1[i] = extq[i];
		extq[i] = (V_t[i] - phi_cond[i]) / CC_a[i][0] - K_t;
    }
    // Potential summation
	cudaMemset((void *) TotPotential, 0.0, Gsize * sizeof(float));
    for(i=0;i<CondNUMR;i++){
        LoadAT2D<<<FIELD_GRID2,FIELD_BLOCK2>>>(Lap_PHI_Sol, pitch, i, dev_phi_buf, Gsize);
		VectorSum<<<FIELD_GRID2, FIELD_BLOCK2>>>(Gsize, TotPotential, phi_cond[i], dev_phi_buf);
		cudaDeviceSynchronize();
    }
	// TotPotential = laplace * votage
    cudaMemcpy(LapPotential, TotPotential, Gsize * sizeof(float),cudaMemcpyDeviceToDevice);
	// TotPotential = laplace * votage 
	// LapPotential = laplace * votage 
	VectorSum<<<FIELD_GRID2, FIELD_BLOCK2>>>(Gsize, TotPotential, 1, dev_phi);
	// TotPotential = laplace * votage + Poisson's
	// LapPotential = laplace * votage 
    GGACopy_Potential<<<FIELD_GRID2, FIELD_BLOCK2>>>(Gsize, dev_GvecSet, LapPotential, dev_phi);
	Cond_Sigma<<<FIELD_GRID2,FIELD_BLOCK2>>>(ngx, ngy,  hdx, hdy, idxy, zlength, dev_GvecSet, dev_CvecSet, dev_info_sp, dev_G_sp, dev_phi, dev_Sigma);
    VtoEfield<<<FIELD_GRID2, FIELD_BLOCK2>>>(ngx,ngy,dx,dy,hdx,hdy,idx,idy, dev_Sigma, TotPotential, dev_G_sp, dev_CvecSet, dev_GvecSet);
    //cudaMemcpy(vec_G, dev_GvecSet, Gsize * sizeof(GGA), cudaMemcpyDeviceToHost); // for TEST
    //Main_Variable_printorSave(); // for TEST
	cudaDeviceSynchronize();
}

static float * Pois_SIG_Sol_tmp;

static void __global__ update
(
	float * Pois_SIG_Sol_tmp,
	GGA * const Field, 
	float * const dev_phi_buf, 
	int Gsize,
	int CondNUMR
)
{
	cg::thread_block group_block  = cg::this_thread_block();
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	__shared__ float a_shd[128];
	
	int id_cond = 0;
	float area;
	float potential;
	if(tid < Gsize)
	{
		id_cond = Field[tid].CondID;
		area = Field[tid].Area;
		potential = dev_phi_buf[tid];
		if(id_cond != 0)
			atomicAdd(Pois_SIG_Sol_tmp + id_cond - 1, area * potential);
	}
}

void PCG_SOLVER(){
    int i;
	void *kernelArgs[] = {
        (void*)&dev_Ai,(void*)&dev_Aj,(void*)&dev_A,(void*)&dev_X,
        (void*)&dev_M, (void*)&dev_AP,(void*)&dev_P,(void*)&dev_R,
        (void*)&dev_Z, (void*)&N,     (void*)&nz,   (void*)&PCGtol2,
        (void*)&FIter, (void*)&dot_result2,
    };
	checkCudaErrors(cudaMemset((void *) dot_result2,  0, 1e3 * sizeof(float)));
	cudaLaunchCooperativeKernel((void *)PCG_float,FIELD_GRID,FIELD_BLOCK, kernelArgs, sMemSize, NULL);
    cudaDeviceSynchronize();
    PCG_Deposit<<<FIELD_GRID2,FIELD_BLOCK2>>>(Gsize, dev_A_idx, dev_GvecSet, dev_X, dev_phi);
    cudaMemset((void *) dev_phi_buf, 0.0, Gsize * sizeof(float));
    Cond_Sigma<<<FIELD_GRID2,FIELD_BLOCK2>>>(ngx, ngy, hdx, hdy, idxy, zlength, dev_GvecSet, dev_CvecSet, dev_info_sp, dev_G_sp, dev_phi, dev_phi_buf);
	cudaMemcpy(Host_G_buf, dev_phi_buf, Gsize * sizeof(float),cudaMemcpyDeviceToHost);
	VFInit(Pois_SIG_Sol,0.0,CondNUMR); 
	cudaDeviceSynchronize();
	
	if(false)
	{
		for (i = 0; i < Gsize; i++) 
		if (vec_G[i].CondID) 
		Pois_SIG_Sol[vec_G[i].CondID - 1] += Host_G_buf[i] * vec_G[i].Area;
		for(int j=0;j<CondNUMR;j++)
			printf("%.12f ", Pois_SIG_Sol[j]);
		puts("");
	}

	if(true){
		dim3 dim_num_block = dim3(Gsize / 128 + 1);
		dim3 dim_size_block = dim3(128);
		checkCudaErrors(cudaMemset((void *) Pois_SIG_Sol_tmp, 0.0f, CondNUMR * sizeof(float)));
		
		update<<<dim_num_block, dim_size_block>>>
		(
			Pois_SIG_Sol_tmp, dev_GvecSet,
			dev_phi_buf, Gsize, CondNUMR
		);
	
		cudaDeviceSynchronize();
		checkCudaErrors
		(
			cudaMemcpy
			(
				Pois_SIG_Sol, Pois_SIG_Sol_tmp, CondNUMR * sizeof(float), cudaMemcpyDeviceToHost
			)
		);
		
		//printf("%.12f %.12f %.12f\n", Pois_SIG_Sol[0],Pois_SIG_Sol[1],Pois_SIG_Sol[2]);
	}

}
void PCG_SOLVER_Laplace(){
    int j,k;
    cudaEvent_t start, stop; // SPEED TEST
    float gputime = 0.0f; // SPEED TEST
    void *kernelArgs[] = {
        (void*)&dev_Ai,(void*)&dev_Aj,(void*)&dev_A,(void*)&dev_X,
        (void*)&dev_M, (void*)&dev_AP,(void*)&dev_P,(void*)&dev_R,
        (void*)&dev_Z, (void*)&N,     (void*)&nz,   (void*)&PCGtol2,
        (void*)&FIter, (void*)&dot_result2,
    };
    for (k = 0; k < CondNUMR; k++) {
        printf(" Laplace Solution %d",k);
        checkCudaErrors(cudaMemcpy(dev_R, cond_b[k], N * sizeof(float),cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemset((void *) dev_X, 0, N * sizeof(float)));
        cudaEventCreate(&start); cudaEventCreate(&stop);
	    cudaEventRecord( start, 0 );
		checkCudaErrors(cudaMemset((void *) dot_result2,  0, 1e4 * sizeof(float)));
        checkCudaErrors(cudaLaunchCooperativeKernel((void *)PCG_float,FIELD_GRID,FIELD_BLOCK, kernelArgs, sMemSize, NULL));
        checkCudaErrors(cudaDeviceSynchronize());
        cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    cudaEventElapsedTime( &gputime, start, stop );
	    cudaEventDestroy( start );cudaEventDestroy( stop );
        printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k+1);
        printf(" - Iter = %d, time = %2.3g (ms), rsold^2 = %g\n",*FIter,gputime,*dot_result2);
        //
        checkCudaErrors(cudaMemset((void *) dev_phi, 0.0, Gsize * sizeof(float)));
        PCG_Deposit_Lap<<<FIELD_GRID2,FIELD_BLOCK2>>>(Gsize, dev_A_idx, dev_GvecSet, k, dev_X, dev_phi);
        //
        checkCudaErrors(cudaMemset((void *) dev_phi_buf, 0.0, Gsize * sizeof(float)));
        Cond_Sigma_Lap<<<FIELD_GRID2,FIELD_BLOCK2>>>(ngx, ngy, dx, dy, zlength, dev_GvecSet, dev_CvecSet, dev_phi, dev_phi_buf);
        VFInit(Host_G_buf,0.0,Gsize);
        checkCudaErrors(cudaMemcpy(Host_G_buf, dev_phi_buf, Gsize * sizeof(float),cudaMemcpyDeviceToHost));
		for (j = 0; j < Gsize; j++) {
			if (vec_G[j].CondID) Lap_SIG_Sol[k][vec_G[j].CondID - 1] += Host_G_buf[j] * vec_G[j].Area;
		}
        for (j = 0; j < CondNUMR; j++) printf(" - Lap_SIG_Sol[%d][%d]= %g\n", k, j, Lap_SIG_Sol[k][j]);
        SaveAT2D<<<FIELD_GRID2,FIELD_BLOCK2>>>(Lap_PHI_Sol, pitch, k, dev_phi, Gsize);
    }
	checkCudaErrors(cudaMemset((void *) dev_phi_buf, 0.0, Gsize * sizeof(float)));
	if(DumpFlag == 0){
    	printf("/***********Calculate temperature distribution**********/\n");
    	checkCudaErrors(cudaMemcpy(dev_R, dev_Tb, N * sizeof(float),cudaMemcpyDeviceToDevice));
    	checkCudaErrors(cudaMemset((void *) dev_X, 0, N * sizeof(float)));
		
		checkCudaErrors(cudaMemset((void *) dot_result2,  0, 1e4 * sizeof(float)));
		checkCudaErrors(cudaLaunchCooperativeKernel((void *)PCG_float,FIELD_GRID,FIELD_BLOCK, kernelArgs, sMemSize, NULL));
		printf(" - Iter = %d, rsold^2 = %g\n",*FIter,*dot_result2);
		printf("/*******************************************************/\n");
    	PCG_Deposit_Temp<<<FIELD_GRID2,FIELD_BLOCK2>>>(Gsize, dev_A_idx, dev_X, dev_GvecSet);
    	if(MainGas == ARGON || MainGas == OXYGEN) Calculate_1GasPara<<<FIELD_GRID2,FIELD_BLOCK2>>>(Gsize, BG[0].mass, BG[0].Pres, dev_GvecSet); 
    	else if(MainGas == ARO2) Calculate_2GasPara<<<FIELD_GRID2,FIELD_BLOCK2>>>(Gsize, BG[0].mass, BG[0].Pres, BG[1].mass, BG[1].Pres, dev_GvecSet);
    	checkCudaErrors(cudaMemcpy(vec_G, dev_GvecSet, Gsize * sizeof(GGA), cudaMemcpyDeviceToHost));
	}
	checkCudaErrors(cudaMemset((void *) dev_X, 0, N * sizeof(float)));
}
void Set_MatrixPCG_cuda(){
    int i,j;
    N = A_size;
    nz = 5 * A_size;
    printf(" Field Solver : [GPU] Preconditioned Conjugate Gradient\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Matrix Size = %d, ngx x ngy = %d X %d = %d\n", N, ngx, ngy, Gsize);
    // Real Solution
    // Laplace Solution
    checkCudaErrors(cudaMallocPitch(&Lap_PHI_Sol, &pitch, Gsize * sizeof(float), CondNUMR)); // for Laplace Solution
    checkCudaErrors(cudaMalloc((void**) &dev_phi, Gsize * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_phi, 0.0, Gsize * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_phi_buf, Gsize * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_phi_buf, 0.0, Gsize * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &Pois_SIG_Sol_tmp, CondNUMR * sizeof(float)));
    Lap_SIG_Sol = MFMalloc(CondNUMR,CondNUMR);
    MFInit(Lap_SIG_Sol,0.0,CondNUMR,CondNUMR);
    Pois_SIG_Sol = VFMalloc(CondNUMR);
    VFInit(Pois_SIG_Sol,0.0,CondNUMR);
    // Allocate
    checkCudaErrors(cudaMalloc((void**) &dev_A, nz * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_Aj, nz * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_Ai, (N + 1) * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_b,  N * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_b, 0.0, N * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_X,  N * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_X, 0.0, N * sizeof(float)));
    //
    checkCudaErrors(cudaMalloc((void**) &dev_TA, nz * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_Tb,  N * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_Tb, 0.0, N * sizeof(float)));
    //
    checkCudaErrors(cudaMalloc((void**) &dev_AP,  N * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_AP, 0.0, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_R,  N * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_R, 0.0, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_P,  N * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_P, 0.0, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_M,  N * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_M, 0.0, N * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_Z,  N * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_Z, 0.0, N * sizeof(float)));
    // Data cpu > gpu
    checkCudaErrors(cudaMemcpy(dev_A, A_val, nz * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(dev_TA,TA_val,nz*sizeof(float),cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Aj, Aj, nz * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Ai, Ai, (N + 1) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(dev_Tb,temp_b, N * sizeof(float),cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(dev_M, MatM, N * sizeof(float), cudaMemcpyHostToDevice));
    // GGA GCA copy
    checkCudaErrors(cudaMalloc((void**)&dev_GvecSet, Gsize * sizeof(GGA)));
    checkCudaErrors(cudaMemcpy(dev_GvecSet, vec_G, Gsize * sizeof(GGA), cudaMemcpyHostToDevice));
    for(i=0;i<Csize;i++){
      vec_C[i].eps = EPS0 * vec_C[i].eps_r;
    }
    checkCudaErrors(cudaMalloc((void**)&dev_CvecSet, Csize * sizeof(GCA)));
    checkCudaErrors(cudaMemcpy(dev_CvecSet, vec_C, Csize * sizeof(GCA), cudaMemcpyHostToDevice));
    //Unified memory value for Field residual
    cudaMallocManaged((void **)&dot_result, sizeof(double));
    *dot_result = 0.0;
	cudaMallocManaged((void **)&dot_result2, sizeof(float) * 10000);
    *dot_result2 = 0.0;
    cudaMallocManaged((void **)&FIter, sizeof(int));
    *FIter = 0;
    //
    vec_A_idx = (int *) malloc(Gsize * sizeof(int));
    for (i = 0; i < ngx; i++) {
		for (j = 0; j < ngy; j++) {
			vec_A_idx[j + i * ngy] = A_idx[i][j];
		} 
	}
    checkCudaErrors(cudaMalloc((void**) &dev_A_idx, Gsize * sizeof(int)));
    checkCudaErrors(cudaMemcpy(dev_A_idx, vec_A_idx, Gsize * sizeof(int),cudaMemcpyHostToDevice));
    // dev_Sigma, dev_Source, TotPotential
    checkCudaErrors(cudaMalloc((void**) &dev_Sigma, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *) dev_Sigma, 0, Gsize * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_Source, Gsize * sizeof(float))); // Total charge --> dev_b or dev_R
	checkCudaErrors(cudaMemset((void *) dev_Source, 0, Gsize * sizeof(float))); // Total charge --> dev_b or dev_R
    checkCudaErrors(cudaMalloc((void**) &TotPotential, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *) TotPotential, 0, Gsize * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &LapPotential, Gsize * sizeof(float)));
	checkCudaErrors(cudaMemset((void *) LapPotential, 0, Gsize * sizeof(float)));
	// Condunctor 
	checkCudaErrors(cudaMalloc((void**)&dev_CondCharge, nsp * CondNUMR * sizeof(float)));
	checkCudaErrors(cudaMemset((void *)dev_CondCharge, 0, nsp * CondNUMR * sizeof(float)));
}
__global__ void GGACopy_Potential(int Gsize, GGA *vecG, float *V1, float *V2){
	int TID = threadIdx.x + blockIdx.x * blockDim.x;
	if(TID>=Gsize) return;
	vecG[TID].Lap_Pot=V1[TID];
    vecG[TID].Pois_Pot=V2[TID];
}
__global__ void VectorSum(int Gsize,float *TotPhi,float V,float *Phi){
	int TID = threadIdx.x + blockIdx.x * blockDim.x;
	if(TID>=Gsize) return;
	//if(TID<5) printf(" VectorSum:TID[%d]\n",TID);
	TotPhi[TID]+=V*Phi[TID];
}
__device__ void Mat_x_Vec(int *I, int *J, float *val, int nnz, int num_rows, float alpha, float *inputVecX, 
                        float *outputVecY, cg::thread_block &cta, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < num_rows; i+= grid.size())    {
        int row_elem = I[i];
        int next_row_elem = I[i+1];
        float output = 0.0;
        for (int j=row_elem-1; j < next_row_elem-1; j++){
            //if(i==0) printf("val[%d][]\n",j);
            // I or J or val arrays - can be put in shared memory 
            // as the access is random and reused in next calls of gpuSpMV function.
            output +=  alpha*val[j] * inputVecX[J[j]-1];
            //if(i==0) printf("val[%d][%d] = %g, %g, %g\n",j,J[j]-1,val[j],inputVecX[J[j]-1],output);
        }
        outputVecY[i] = output;
    }
}
__device__ void A_x_X_p_Y(float a, float *x, float *y, int size, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < size; i+= grid.size()) y[i] = a*x[i] + y[i];
}
__device__ void A_x_Y_p_X(float a, float *x, float *y, int size, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < size; i+= grid.size()) y[i] = a*y[i] + x[i];
}
__device__ void Vec_Dot_Sum(float *vecA, float *vecB, double *result, int size, const cg::thread_block &cta, const cg::grid_group &grid)
{
   __shared__ double tmp[THREADS_PER_BLOCK];
    double temp_sum = 0.0;
    for (int i=grid.thread_rank(); i < size; i+=grid.size()){
        temp_sum += (double) (vecA[i] * vecB[i]);
    }
    tmp[cta.thread_rank()] = temp_sum;
    cg::sync(cta);
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);
    double beta  = temp_sum;
    double temp;
    for (int i = tile32.size() / 2; i > 0; i >>= 1) {
        if (tile32.thread_rank() < i) {
            temp       = tmp[cta.thread_rank() + i];
            beta       += temp;
            tmp[cta.thread_rank()] = beta;
        }
        cg::sync(tile32);
    }
    cg::sync(cta);
    if (cta.thread_rank() == 0) {
        beta  = 0.0;
        for (int i = 0; i < cta.size(); i += tile32.size()) {
            beta  += tmp[i];
        }
        atomicAdd(result, beta);
    }
}
__device__ void CopyVector(float *srcA, float *destB, int size, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < size; i+= grid.size()) destB[i] = srcA[i];
}
__device__ void Vec_x_Vec(float *vecA, float *vecB, float *vecC, int size, const cg::thread_block &cta, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < size; i+=grid.size()) 
		vecC[i] = (vecA[i] * vecB[i]);
}
__global__ void PCG(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, double *d_result){
    //Jacovi diagonal preconditioner version
    cg::thread_block cta = cg::this_thread_block();
    cg::grid_group grid = cg::this_grid();
    //int TID = blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    int max_iter = 100000;
    float a = 1.0;
    float na = -1.0;
    float rsold,rnew,Temp;
    float nalpha,alpha,beta;
    rsold = 0.0;
    if (threadIdx.x == 0 && blockIdx.x == 0){
        *Iter = 0;
        *d_result = 0.0;  
    } 
    Mat_x_Vec(I, J, val, nnz, N, a, x, Ax, cta, grid); 
    A_x_X_p_Y(na, Ax, r, N, grid); 
    Vec_x_Vec(M, r, Z, N, cta, grid);
    CopyVector(Z, p, N, grid);
    cg::sync(grid);
    Vec_Dot_Sum(r, Z, d_result, N, cta, grid); 
    cg::sync(grid);
    rsold = *d_result;
    while (rsold > tol2 && *Iter <= max_iter){
        Mat_x_Vec(I, J, val, nnz, N, a, p, Ax, cta, grid);
        if (threadIdx.x == 0 && blockIdx.x == 0){
            *Iter = *Iter + 1;
            *d_result = 0.0;  
        } 
        cg::sync(grid);
        Vec_Dot_Sum(p, Ax, d_result, N, cta, grid);
        cg::sync(grid);
        Temp = *d_result;
        alpha = (Temp)? rsold/Temp:0.0f;
		//if(threadIdx.x == 0 && blockIdx.x == 0) printf("Iter = %d, temp = %g,  Res = %g AL = %g\n",*Iter,Temp,rsold,alpha);
        A_x_X_p_Y(alpha, p, x, N, grid);
        nalpha = -alpha;
        A_x_X_p_Y(nalpha, Ax, r, N, grid);
        Vec_x_Vec(M, r, Z, N, cta, grid);
        if (threadIdx.x == 0 && blockIdx.x == 0){
            *d_result = 0.0;  
        } 
        cg::sync(grid);
        Vec_Dot_Sum(r, Z, d_result, N, cta, grid);
        cg::sync(grid);
        rnew = *d_result;
        beta = (rsold) ? rnew/rsold: 0.0f;
        A_x_Y_p_X(beta, Z, p, N, grid);
        rsold = rnew;
        //rnew = 0.0;
    }
    //if(threadIdx.x == 0 && blockIdx.x == 0 ) printf("End Iter = %d, Res = %g, b = %g, a = %g\n",*Iter,Temp,alpha,beta,rsold);
}
__global__ void Cond_Sigma_Lap(int ngx, int ngy, float dx, float dy, float zlength, GGA *vecG, GCA *vecC, float *Phi, float *Sigma)
{
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
	int x, y;
	int Phi_left, Phi_right, Phi_up, Phi_down;
	int EPS_TID, EPS_left, EPS_down, EPS_cross, EPS_xover, EPS_yover;

	if(TID>=ngx*ngy) return;
	x=TID/ngy; y=TID%ngy;
	//Calculate surface charge
	EPS_TID=x*(ngy-1)+y;
	EPS_left=(x) ? EPS_TID-ngy+1:EPS_TID;
	EPS_down=(y) ? EPS_TID-1:EPS_TID;
	EPS_cross=EPS_left+EPS_down-EPS_TID;
	EPS_xover=(x==ngx-1) ? ngy-1: 0;
	EPS_yover=(y==ngy-1) ? 1: 0;
	EPS_TID-=(EPS_xover+EPS_yover);
	EPS_left-=EPS_yover;
	EPS_down-=EPS_xover;

	Phi_left =(x) ? TID-ngy:TID;
	Phi_right=(x==ngx-1) ? TID:TID+ngy;
	Phi_down =(y) ? TID-1:TID;
	Phi_up =(y==ngy-1) ? TID:TID+1;

	if((vecG[TID].Boundary==CONDUCTOR || vecG[TID].Boundary==DIRICHLET) && vecG[TID].Face!=NO_FACE) {
		if(vecG[TID].Face==UP) {
			Sigma[TID] = 0.5*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[TID]-Phi[Phi_up])/dy;
		}
		else if(vecG[TID].Face==DOWN) {
			Sigma[TID] = 0.5*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_down])/dy;
		}
		else if(vecG[TID].Face==LEFT) {
			Sigma[TID] = 0.5*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_left])/dx;
		}
		else if(vecG[TID].Face==RIGHT) {
			Sigma[TID] = 0.5*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[TID]-Phi[Phi_right])/dx;
		}
		else if(vecG[TID].Face==UL_CORN) {
			Sigma[TID] = 0.5*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[TID]-Phi[Phi_up])/dy+0.5*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_left])/dx;
		}
		else if(vecG[TID].Face==UR_CORN) {
			Sigma[TID] = 0.5*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[TID]-Phi[Phi_up])/dy+0.5*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[TID]-Phi[Phi_right])/dx;
		}
		else if(vecG[TID].Face==LL_CORN) {
			Sigma[TID] = 0.5*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_down])/dy+0.5*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_left])/dx;
		}
		else if(vecG[TID].Face==LR_CORN) {
			Sigma[TID] = 0.5*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_down])/dy+0.5*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[TID]-Phi[Phi_right])/dx;
		}
	}
}
__global__ void Cond_Sigma(int ngx, int ngy, float hdx, float hdy, float idxy, float zlength, 
					GGA *vecG, GCA *vecC, Species *info, GPG *data, float *Phi, float *Sigma)
{
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
	int x, y;
    int i;
    float SumSig;
	int Phi_left, Phi_right, Phi_up, Phi_down;
	int EPS_TID, EPS_left, EPS_down, EPS_cross, EPS_xover, EPS_yover;

	if(TID>=ngx*ngy) return;
	//if(TID<5) printf("Cond_Sigma:TID[%d]\n",TID);
	x=TID/ngy; y=TID%ngy;
	//Calculate surface charge
	EPS_TID=x*(ngy-1)+y;
	EPS_left=(x) ? EPS_TID-ngy+1:EPS_TID;
	EPS_down=(y) ? EPS_TID-1:EPS_TID;
	EPS_cross=EPS_left+EPS_down-EPS_TID;
	EPS_xover=(x==ngx-1) ? ngy-1: 0;
	EPS_yover=(y==ngy-1) ? 1: 0;
	EPS_TID-=(EPS_xover+EPS_yover);
	EPS_left-=EPS_yover;
	EPS_down-=EPS_xover;

	Phi_left =(x) ? TID-ngy:TID;
	Phi_right=(x==ngx-1) ? TID:TID+ngy;
	Phi_down =(y) ? TID-1:TID;
	Phi_up =(y==ngy-1) ? TID:TID+1;

    SumSig = 0;
    for(i=0;i<info[0].spnum;i++) SumSig += data[TID + i*ngx*ngy].den * info[i].q_density;

	if((vecG[TID].Boundary==CONDUCTOR || vecG[TID].Boundary==DIRICHLET) && vecG[TID].Face!=NO_FACE) {
		if(vecG[TID].Face==UP) {
			Sigma[TID] = hdy*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[TID]-Phi[Phi_up])-hdx*SumSig;
		}
		else if(vecG[TID].Face==DOWN) {
			Sigma[TID] = hdy*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_down])-hdx*SumSig;
		}
		else if(vecG[TID].Face==LEFT) {
			Sigma[TID] = hdx*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_left])-hdy*SumSig;
		}
		else if(vecG[TID].Face==RIGHT) {
			Sigma[TID] = hdx*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[TID]-Phi[Phi_right])-hdy*SumSig;
		}
		else if(vecG[TID].Face==UL_CORN) {
			Sigma[TID] = hdy*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[TID]-Phi[Phi_up])+hdx*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_left])-idxy*SumSig;
		}
		else if(vecG[TID].Face==UR_CORN) {
			Sigma[TID] = hdy*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[TID]-Phi[Phi_up])+hdx*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[TID]-Phi[Phi_right])-idxy*SumSig;
		}
		else if(vecG[TID].Face==LL_CORN) {
			Sigma[TID] = hdy*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_down])+hdx*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_left])-idxy*SumSig;
		}
		else if(vecG[TID].Face==LR_CORN) {
			Sigma[TID] = hdy*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_down])+hdx*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[TID]-Phi[Phi_right])-idxy*SumSig;
		}
	}
}
__global__ void Cond_Sigma_v2(int ngx, int ngy, float hdx, float hdy, float idxy, float zlength, 
					GGA *vecG, GCA *vecC, Species *info, GPG *data, float *Phi, float *Sigma)
{
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
	int x, y;
    int i;
    float SumSig;
	int Phi_left, Phi_right, Phi_up, Phi_down;
	int EPS_TID, EPS_left, EPS_down, EPS_cross, EPS_xover, EPS_yover;
	float a,b,c;
	if(TID>=ngx*ngy) return;
	//if(TID<5) printf("Cond_Sigma:TID[%d]\n",TID);
	x=TID/ngy; y=TID%ngy;
	//Calculate surface charge
	EPS_TID=x*(ngy-1)+y;
	EPS_left=(x) ? EPS_TID-ngy+1:EPS_TID;
	EPS_down=(y) ? EPS_TID-1:EPS_TID;
	EPS_cross=EPS_left+EPS_down-EPS_TID;
	EPS_xover=(x==ngx-1) ? ngy-1: 0;
	EPS_yover=(y==ngy-1) ? 1: 0;
	EPS_TID-=(EPS_xover+EPS_yover);
	EPS_left-=EPS_yover;
	EPS_down-=EPS_xover;

	Phi_left =(x) ? TID-ngy:TID;
	Phi_right=(x==ngx-1) ? TID:TID+ngy;
	Phi_down =(y) ? TID-1:TID;
	Phi_up =(y==ngy-1) ? TID:TID+1;

    SumSig = 0;
    for(i=0;i<info[0].spnum;i++) SumSig += data[TID + i*ngx*ngy].den * info[i].q_density;

	if((vecG[TID].Boundary==CONDUCTOR || vecG[TID].Boundary==DIRICHLET) && vecG[TID].Face!=NO_FACE) {
		if(vecG[TID].Face==UP) {
			a = hdy*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[TID]-Phi[Phi_up]);
			b = 0.0f;
			c = hdx*SumSig;
		}
		else if(vecG[TID].Face==DOWN) {
			a = hdy*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_down]);
			b = 0.0f;
			c = hdx*SumSig;
		}
		else if(vecG[TID].Face==LEFT) {
			a = hdx*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_left]);
			b = 0.0f;
			c = hdy*SumSig;
		}
		else if(vecG[TID].Face==RIGHT) {
			a = hdx*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[TID]-Phi[Phi_right]);
			b = 0.0f;
			c = hdy*SumSig;
		}
		else if(vecG[TID].Face==UL_CORN) {
			a = hdy*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[TID]-Phi[Phi_up]);
			b = hdx*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_left]);
			c = idxy*SumSig;
		}
		else if(vecG[TID].Face==UR_CORN) {
			a = hdy*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[TID]-Phi[Phi_up]);
			b = hdx*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[TID]-Phi[Phi_right]);
			c = idxy*SumSig;
		}
		else if(vecG[TID].Face==LL_CORN) {
			a = hdy*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_down]);
			b = hdx*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_left]);
			c = idxy*SumSig;
		}
		else if(vecG[TID].Face==LR_CORN) {
			a = hdy*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[TID]-Phi[Phi_down]);
			b = hdx*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[TID]-Phi[Phi_right]);
			c = idxy*SumSig;
		}else{
			return;
		}
	}else{
		return;
	}
	Sigma[TID] = a + b - c;
}
__global__ void PCG_Deposit(int Gsize, int *IDX, GGA *vecG, float *X, float *PHI){
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
	if(TID>=Gsize) return;
	//if(TID<5) printf("PCG_Deposit:TID[%d]\n",TID);
	if(IDX[TID]) PHI[TID]=X[IDX[TID]-1];
    else PHI[TID] = 0;
}
__global__ void PCG_Deposit_Lap(int Gsize, int *IDX, GGA *vecG, int k, float *X, float *PHI){
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
	if(TID>=Gsize) return;
	if(IDX[TID])
		PHI[TID]=X[IDX[TID]-1];
	else if(vecG[TID].CondID == k+1)
		PHI[TID]=1.0;   
}
__global__ void PCG_Deposit_Temp(int Gsize, int *IDX, float *X, GGA *vecG)
{
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;

	if(TID>=Gsize) return;
	else if(IDX[TID]) vecG[TID].Temp = X[IDX[TID]-1];
}
__global__ void Calculate_1GasPara(int Gsize, float mass, float press, GGA *vecG){
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
	if(TID>=Gsize) return;
	if(vecG[TID].DensRegion){
    	vecG[TID].BackDen1 = NperTORR*press/(vecG[TID].Temp*8.6142e-5+DBL_MIN);
		vecG[TID].BackVel1 = sqrt(vecG[TID].Temp*1.38e-23/mass);
	}
}
__global__ void Calculate_2GasPara(int Gsize, float mass1, float press1, float mass2, float press2, GGA *vecG){
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
	if(TID>=Gsize) return;
	if(vecG[TID].DensRegion){
    	vecG[TID].BackDen1 = NperTORR*press1/(vecG[TID].Temp*8.6142e-5+DBL_MIN);
		vecG[TID].BackVel1 = sqrt(vecG[TID].Temp*1.38e-23/mass1);
    	vecG[TID].BackDen2 = NperTORR*press2/(vecG[TID].Temp*8.6142e-5+DBL_MIN);
		vecG[TID].BackVel2 = sqrt(vecG[TID].Temp*1.38e-23/mass2);
	}
}
__global__ void SaveAT2D(float *A, size_t pitch, int height, float *PHI, int n){
    // High save and load for Matrix type variable 
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;

	if(TID>=n) return;
	float *row=(float *)((char *)A+height*pitch);

	row[TID]=PHI[TID];
}
__global__ void LoadAT2D(float *A, size_t pitch, int height, float *PHI, int n){
    // High save and load for Matrix type variable
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;

	if(TID>=n) return;
	//if(TID<5) printf("LoadAT2D:TID[%d]\n",TID);
	float *row=(float *)((char *)A+height*pitch);

	PHI[TID]=row[TID];
}
__global__ void VtoEfield(int ngx,int ngy,float dx,float dy,float hdx,float hdy,float idx,
                            float idy, float *Sigma, float *Phi, GPG *data, GCA *vecC, GGA *vecG){
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x ;
	int x, y;
	int Phi_left, Phi_right, Phi_up, Phi_down;
	int EPS_TID, EPS_left, EPS_down, EPS_cross, EPS_xover, EPS_yover;

	if(TID>=ngx*ngy) return;
	//if(TID<5) printf("VtoEfield:TID[%d]\n",TID);

	x=TID/ngy; y=TID%ngy;
	Phi_left =(x) ? TID-ngy:TID;
	Phi_right=(x==ngx-1) ? TID:TID+ngy;
	Phi_down =(y) ? TID-1:TID;
	Phi_up =(y==ngy-1) ? TID:TID+1;
	EPS_TID=x*(ngy-1)+y;
	EPS_left=(x) ? EPS_TID-ngy+1:EPS_TID;
	EPS_down=(y) ? EPS_TID-1:EPS_TID;
	EPS_cross=EPS_left+EPS_down-EPS_TID;
	EPS_xover=(x==ngx-1) ? ngy-1: 0;
	EPS_yover=(y==ngy-1) ? 1: 0;
	EPS_TID-=(EPS_xover+EPS_yover);
	EPS_left-=EPS_yover;
	EPS_down-=EPS_xover;

	if(vecG[TID].Face == LEFT) {
		if(vecG[TID].Boundary==CONDUCTOR || vecG[TID].Boundary==NEUMANN || vecG[TID].Boundary==DIRICHLET) {
			vecG[TID].Ex = -2.0*Sigma[TID]/(vecC[EPS_left].eps+vecC[EPS_cross].eps);
			vecG[TID].Ey = 0.0;
		}
		else if(vecG[TID].Boundary==DIELECTRIC || vecG[TID].Boundary>=100) {
			vecG[TID].Ex = idx*(Phi[Phi_left]-Phi[TID]);
			vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
		}
		else {
			vecG[TID].Ex = 0.5*(vecC[EPS_TID].eps+vecC[EPS_down].eps)*(Phi[Phi_right]-Phi[TID])*idx;
			vecG[TID].Ex-= (0.5*dx/dy)*vecC[EPS_TID].eps*(Phi[TID]-Phi[Phi_up])*idy;
			vecG[TID].Ex-= (0.5*dx/dy)*vecC[EPS_down].eps*(Phi[TID]-Phi[Phi_down])*idy;
			vecG[TID].Ex/= -0.5*(vecC[EPS_left].eps+vecC[EPS_cross].eps);
			vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
		}
	}
	else if(vecG[TID].Face == RIGHT) {
		if(vecG[TID].Boundary==CONDUCTOR ||  vecG[TID].Boundary==NEUMANN || vecG[TID].Boundary==DIRICHLET) {
			vecG[TID].Ex = 2.0*Sigma[TID]/(vecC[EPS_TID].eps+vecC[EPS_down].eps);
			vecG[TID].Ey = 0.0;
		}
		else if(vecG[TID].Boundary==DIELECTRIC || vecG[TID].Boundary>=100) {
			vecG[TID].Ex = idx*(Phi[TID]-Phi[Phi_right]);
			vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
		}
		else {
			vecG[TID].Ex = 0.5*(vecC[EPS_left].eps+vecC[EPS_cross].eps)*(Phi[Phi_left]-Phi[TID])*idx;
			vecG[TID].Ex-= (0.5*dx/dy)*vecC[EPS_left].eps*(Phi[TID]-Phi[Phi_up])*idy;
			vecG[TID].Ex-= (0.5*dx/dy)*vecC[EPS_cross].eps*(Phi[TID]-Phi[Phi_down])*idy;
			vecG[TID].Ex/= 0.5*(vecC[EPS_TID].eps+vecC[EPS_down].eps);

			vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
		}
	}
	else if(vecG[TID].Face == UP) {
		if(vecG[TID].Boundary==CONDUCTOR ||  vecG[TID].Boundary==NEUMANN || vecG[TID].Boundary==DIRICHLET) {
			vecG[TID].Ex = 0.0;
			vecG[TID].Ey = 2.0*Sigma[TID]/(vecC[EPS_TID].eps+vecC[EPS_left].eps);
		}
		else if(vecG[TID].Boundary==DIELECTRIC || vecG[TID].Boundary>=100) {
			vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
			vecG[TID].Ey = idy*(Phi[TID]-Phi[Phi_up]);
		}
		else {
			vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
			vecG[TID].Ey = 0.5*(vecC[EPS_down].eps+vecC[EPS_cross].eps)*(Phi[Phi_down]-Phi[TID])*idy;
			vecG[TID].Ey-= (0.5*dy/dx)*vecC[EPS_down].eps*(Phi[TID]-Phi[Phi_right])*idx;
			vecG[TID].Ey-= (0.5*dy/dx)*vecC[EPS_cross].eps*(Phi[TID]-Phi[Phi_left])*idx;
			vecG[TID].Ey/= 0.5*(vecC[EPS_TID].eps+vecC[EPS_left].eps);
		}
	}
	else if(vecG[TID].Face == DOWN) {
		if(vecG[TID].Boundary==CONDUCTOR ||  vecG[TID].Boundary==NEUMANN || vecG[TID].Boundary==DIRICHLET) {
			vecG[TID].Ex = 0.0;
			vecG[TID].Ey = -2.0*Sigma[TID]/(vecC[EPS_down].eps+vecC[EPS_cross].eps);
		}
		else if(vecG[TID].Boundary==DIELECTRIC || vecG[TID].Boundary>=100) {
			vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
			vecG[TID].Ey = idy*(Phi[Phi_down]-Phi[TID]);
		}
		else {
			vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
			vecG[TID].Ey = 0.5*(vecC[EPS_TID].eps+vecC[EPS_left].eps)*(Phi[Phi_up]-Phi[TID])*idy;
			vecG[TID].Ey-= (0.5*dy/dx)*vecC[EPS_TID].eps*(Phi[TID]-Phi[Phi_right])*idx;
			vecG[TID].Ey-= (0.5*dy/dx)*vecC[EPS_left].eps*(Phi[TID]-Phi[Phi_left])*idx;
			vecG[TID].Ey/= -0.5*(vecC[EPS_down].eps+vecC[EPS_cross].eps);
		}
	}
	else if(vecG[TID].Face == UL_CORN) {
		if(vecG[TID].Boundary==CONDUCTOR || vecG[TID].Boundary==DIRICHLET || vecG[TID].Boundary==DIELECTRIC) {
			vecG[TID].Ex = idx*(Phi[Phi_left]-Phi[TID]);
			vecG[TID].Ey = idy*(Phi[TID]-Phi[Phi_up]);
		}
		else {
			vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
			vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
		}
	}
	else if(vecG[TID].Face == UR_CORN) {
		if(vecG[TID].Boundary==CONDUCTOR || vecG[TID].Boundary==DIRICHLET || vecG[TID].Boundary==DIELECTRIC) {
			vecG[TID].Ex = idx*(Phi[TID]-Phi[Phi_right]);
			vecG[TID].Ey = idy*(Phi[TID]-Phi[Phi_up]);
		}
		else {
			vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
			vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
		}
	}
	else if(vecG[TID].Face == LL_CORN) {
		if(vecG[TID].Boundary==CONDUCTOR || vecG[TID].Boundary==DIRICHLET || vecG[TID].Boundary==DIELECTRIC) {
			vecG[TID].Ex = idx*(Phi[Phi_left]-Phi[TID]);
			vecG[TID].Ey = idy*(Phi[Phi_down]-Phi[TID]);
		}
		else {
			vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
			vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
		}
	}
	else if(vecG[TID].Face == LR_CORN) {
		if(vecG[TID].Boundary==CONDUCTOR || vecG[TID].Boundary==DIRICHLET || vecG[TID].Boundary==DIELECTRIC) {
			vecG[TID].Ex = idx*(Phi[TID]-Phi[Phi_right]);
			vecG[TID].Ey = idy*(Phi[Phi_down]-Phi[TID]);
		}
		else {
			vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
			vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
		}
	}
	else if(vecG[TID].Boundary==NEUMANN) {
		if(x==0 || x==ngx-1) {
			vecG[TID].Ex = 0.0;
			vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
		}
		else if(y==0 || y==ngy-1) {
			vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
			vecG[TID].Ey = 0.0;
		}
	}
	else {
		vecG[TID].Ex = hdx*(Phi[Phi_left]-Phi[Phi_right]);
		vecG[TID].Ey = hdy*(Phi[Phi_down]-Phi[Phi_up]);
	}
}
void transpose(float **matrix, float matrix_cofactor[N_MAX][N_MAX],
		float size) {
	int i, j;
	float m_transpose[N_MAX][N_MAX], m_inverse[N_MAX][N_MAX], d;

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			m_transpose[i][j] = matrix_cofactor[j][i];
		}
	}
	d = determinant(matrix, size);
	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			m_inverse[i][j] = m_transpose[i][j] / d;
		}
	}

	for (i = 0; i < size; i++) {
		for (j = 0; j < size; j++) {
			matrix[i][j] = m_inverse[i][j];
		}
	}
}
float determinant(float **matrix, float size) {
	float s = 1, det = 0;
	float **m_minor;
	int i, j, m, n, c;

	if (size == 1) {
		return (matrix[0][0]);
	} else {
		m_minor = (float **) malloc(size * sizeof(float *));
		for (i = 0; i < size; i++)
			m_minor[i] = (float *) malloc(size * sizeof(float));
		det = 0;
		for (c = 0; c < size; c++) {
			m = 0;
			n = 0;
			for (i = 0; i < size; i++) {
				for (j = 0; j < size; j++) {
					m_minor[i][j] = 0;
					if (i != 0 && j != c) {
						m_minor[m][n] = matrix[i][j];
						if (n < (size - 2))
							n++;
						else {
							n = 0;
							m++;
						}
					}
				}
			}
			det = det + s * (matrix[0][c] * determinant(m_minor, size - 1));
			s = -1 * s;
		}
	}
	for (i = 0; i < size; i++)
		free(m_minor[i]);
	free(m_minor);
	return (det);
}
void cofactor(float **matrix, float size) {
	float matrix_cofactor[N_MAX][N_MAX];
	float **m_cofactor;
	int p, q, m, n, i, j;

	m_cofactor = (float **) malloc(size * sizeof(float *));
	for (i = 0; i < size; i++)
		m_cofactor[i] = (float *) malloc(size * sizeof(float));

	for (q = 0; q < size; q++) {
		for (p = 0; p < size; p++) {
			m = 0;
			n = 0;
			for (i = 0; i < size; i++) {
				for (j = 0; j < size; j++) {
					if (i != q && j != p) {
						m_cofactor[m][n] = matrix[i][j];
						if (n < (size - 2))
							n++;
						else {
							n = 0;
							m++;
						}
					}
				}
			}
			matrix_cofactor[q][p] = pow(-1, q + p) * determinant(m_cofactor, size - 1);
		}
	}
	transpose(matrix, matrix_cofactor, size);
	for (i = 0; i < size; i++)
		free(m_cofactor[i]);
	free(m_cofactor);
}
__global__ void PCG_float2(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, float *d_result){
    //Jacovi diagonal preconditioner version
    cg::thread_block cta = cg::this_thread_block();
    cg::grid_group grid = cg::this_grid();
    //int TID = blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    int max_iter = 100000;
    float a = 1.0f;
    float na = -1.0f;
    float rsold,rnew,Temp;
    float nalpha,alpha,beta;
    rsold = 0.0f;
    if (threadIdx.x == 0 && blockIdx.x == 0){
        *Iter = 0;
        *d_result = 0.0f;  
    } 
    Mat_x_Vec(I, J, val, nnz, N, a, x, Ax, cta, grid); 
    A_x_X_p_Y(na, Ax, r, N, grid); 
    Vec_x_Vec(M, r, Z, N, cta, grid);
    CopyVector(Z, p, N, grid);
    cg::sync(grid);
    Vec_Dot_Sum_F(r, Z, d_result, N, cta, grid); 
    cg::sync(grid);
    rsold = *d_result;
    while (rsold > tol2 && *Iter <= max_iter){
        Mat_x_Vec(I, J, val, nnz, N, a, p, Ax, cta, grid);
        if (threadIdx.x == 0 && blockIdx.x == 0){
            *Iter = *Iter + 1;
            *d_result = 0.0f;  
        }
        cg::sync(grid);
        Vec_Dot_Sum_F(p, Ax, d_result, N, cta, grid);
        cg::sync(grid);
        Temp = *d_result;
        alpha = (Temp)? rsold/Temp:0.0f;
		A_x_X_p_Y(alpha, p, x, N, grid);
        nalpha = -alpha;
        A_x_X_p_Y(nalpha, Ax, r, N, grid);
        Vec_x_Vec(M, r, Z, N, cta, grid);
        if (threadIdx.x == 0 && blockIdx.x == 0){
            *d_result = 0.0f;  
        } 
        cg::sync(grid);
        Vec_Dot_Sum_F(r, Z, d_result, N, cta, grid);
        cg::sync(grid);
        rnew = *d_result;
        beta = (rsold) ? rnew/rsold: 0.0f;
        A_x_Y_p_X(beta, Z, p, N, grid);
        rsold = rnew;
        //rnew = 0.0;
    }
    //if(threadIdx.x == 0 && blockIdx.x == 0 ) printf("End Iter = %d, Res = %g, b = %g, a = %g\n",*Iter,Temp,alpha,beta,rsold);
}
__global__ void PCG_float(int *I, int *J, float *val, float *x, float *M, float *Ax, float *p, float *r, float *Z, 
            int N, int nnz, float tol2, int *Iter, float *d_result){
    //Jacovi diagonal preconditioner version
    cg::thread_block cta = cg::this_thread_block();
    cg::grid_group grid = cg::this_grid();
    int max_iter = 100000;
    float a = 1.0f;
    float na = -1.0f;
    float rsold,rnew,Temp;
    float nalpha,alpha,beta;
    rsold = 0.0f;
    if (threadIdx.x == 0 && blockIdx.x == 0){
        *Iter = 0;
    } 
    Mat_x_Vec(I, J, val, nnz, N, a, x, Ax, cta, grid); 
    A_x_X_p_Y(na, Ax, r, N, grid); 
    Vec_x_Vec(M, r, Z, N, cta, grid);
    CopyVector(Z, p, N, grid);
    cg::sync(grid);
    Vec_Dot_Sum_F(r, Z, d_result, N, cta, grid); 
    cg::sync(grid);
    rsold = *d_result;

	float p_local, r_local, z_local;
	int iter_local = 0;
	
	int cnt;

	float val_reg_0, val_reg_1, val_reg_2, val_reg_3, val_reg_4;
	int J_reg_0, J_reg_1, J_reg_2, J_reg_3, J_reg_4;
	float M_reg_0;

	cnt = 0;
	for (int i=grid.thread_rank(); i < N; i+= grid.size())
	{
		if(cnt == 0)
		{
			M_reg_0 = M[i];
		}
		++cnt;
	}

	cnt = 0;
	for (int i=grid.thread_rank(); i < N; i+= grid.size())    
	{
		int row_elem = I[i];
		int next_row_elem = I[i+1];
		float output = 0.0;
		for (int j=row_elem-1; j < next_row_elem-1; j++){
			if(cnt == 0)
			{
				val_reg_0 = val[j];
				J_reg_0 = J[j];
			}
			else if(cnt == 1)
			{
				val_reg_1 = val[j];
				J_reg_1 = J[j];
			}
			else if(cnt == 2)
			{
				val_reg_2 = val[j];
				J_reg_2 = J[j];
			}
			else if(cnt == 3)
			{
				val_reg_3 = val[j];
				J_reg_3 = J[j];
			}
			else if(cnt == 4)
			{
				val_reg_4 = val[j];
				J_reg_4 = J[j];
			}
			++cnt;
		}
	}

    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);
	
	Mat_x_Vec(I, J, val, nnz, N, a, p, Ax, cta, grid);
	iter_local++;
	Vec_Dot_Sum_F(p, Ax, d_result + iter_local * 2 - 1, N, cta, grid);
	cg::sync(grid);

    while(rsold > tol2 && iter_local <= max_iter)
	{
        Temp = d_result[iter_local * 2 - 1];
        alpha = (Temp)? rsold/Temp:0.0f;
		nalpha = -alpha;

        //A_x_X_p_Y(nalpha, Ax, r, N, grid);
        //Vec_x_Vec(M, r, Z, N, cta, grid);
		//Vec_Dot_Sum_F(r, Z, d_result + 2 * iter_local, N, cta, grid);
		{
			float temp_sum = 0.0f;
			cnt = 0;
			for (int i=grid.thread_rank(); i < N; i+= grid.size())
			{
				r_local = nalpha*Ax[i] + r[i];
				if(cnt == 0)
				{
					z_local = M_reg_0 * r_local;
				}
				else
				{
					z_local = M[i] * r_local;
				}
				temp_sum += z_local * r_local;
				r[i] = r_local;
				Z[i] = z_local;
			}
			for (int offset = 16; offset > 0; offset /= 2)
				temp_sum += __shfl_down_sync(0xffffffff, temp_sum, offset);
			if (tile32.thread_rank() == 0) 
			{
				atomicAdd(d_result + 2 * iter_local, temp_sum);
			}
		}

		cg::sync(grid);
        rnew = d_result[2 * iter_local];
        beta = (rsold) ? rnew/rsold: 0.0f;
		//A_x_X_p_Y(alpha, p, x, N, grid);
		//A_x_Y_p_X(beta, Z, p, N, grid);
        //Mat_x_Vec(I, J, val, nnz, N, a, p, Ax, cta, grid);
		{
			for (int i=grid.thread_rank(); i < N; i+= grid.size()) 
			{
				p_local = p[i];
				x[i] = alpha*p_local + x[i];
				p[i] = beta*p_local + Z[i];
			}
		}
		rsold = rnew;
		{
			cnt = 0;
			for (int i=grid.thread_rank(); i < N; i+= grid.size())
			{
				int row_elem = I[i];
				int next_row_elem = I[i+1];
				float output = 0.0;
				for (int j=row_elem-1; j < next_row_elem-1; j++){
					if(cnt == 0)
					{
						output +=  val_reg_0 * p[J_reg_0-1];
					}
					else if(cnt == 1)
					{
						output +=  val_reg_1 * p[J_reg_1-1];
					}
					else if(cnt == 2)
					{
						output +=  val_reg_2 * p[J_reg_2-1];
					}
					else if(cnt == 3)
					{
						output +=  val_reg_3 * p[J_reg_3-1];
					}
					else if(cnt == 4)
					{
						output +=  val_reg_4 * p[J_reg_4-1];
					}
					else
					{
						output += val[j] * p[J[j]-1];
					}
					cnt++;
				}
				Ax[i] = a * output;
			}
		}
		iter_local++;
        //Vec_Dot_Sum_F(p, Ax, d_result + iter_local * 2 - 1, N, cta, grid);
		{
			float temp_sum = 0.0f;
			cnt = 0;
			for (int i=grid.thread_rank(); i < N; i+= grid.size())
			{
				temp_sum += Ax[i] * p[i];
			}
			for (int offset = 16; offset > 0; offset /= 2)
				temp_sum += __shfl_down_sync(0xffffffff, temp_sum, offset);
			if (tile32.thread_rank() == 0) 
			{
				atomicAdd(d_result + 2 * iter_local - 1, temp_sum);
			}
		}
        cg::sync(grid);
    }
	if (threadIdx.x == 0 && blockIdx.x == 0){
		d_result[0] = d_result[2 * iter_local - 2];
		*Iter = iter_local;
	}
}
__device__ void Vec_Dot_Sum_F(float *vecA, float *vecB, float *result, int size, const cg::thread_block &cta, const cg::grid_group &grid)
{
   __shared__ float tmp[THREADS_PER_BLOCK];
    float temp_sum = 0.0;
    for (int i=grid.thread_rank(); i < size; i+=grid.size()){
        temp_sum += (float) (vecA[i] * vecB[i]);
    }
    tmp[cta.thread_rank()] = temp_sum;
    cg::sync(cta);
    cg::thread_block_tile<32> tile32 = cg::tiled_partition<32>(cta);
    float beta  = temp_sum;
    float temp;
    for (int i = tile32.size() / 2; i > 0; i >>= 1) {
        if (tile32.thread_rank() < i) {
            temp       = tmp[cta.thread_rank() + i];
            beta       += temp;
            tmp[cta.thread_rank()] = beta;
        }
        cg::sync(tile32);
    }
    cg::sync(cta);
    if (cta.thread_rank() == 0) {
        beta  = 0.0;
        for (int i = 0; i < cta.size(); i += tile32.size()) {
            beta  += tmp[i];
        }
        atomicAdd(result, beta);
    }
}