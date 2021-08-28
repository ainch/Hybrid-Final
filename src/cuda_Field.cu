#include "cuda_Field.cuh"
// FOR Field method 4
#define THREADS_PER_BLOCK 512   
//
void PCG_SOLVER_Laplace(){
    // OUTPUT
    // Lap_TEMP_Sol[Gsize] : Temperature Profile
    // Lap_PHI_Sol[CondNUMR][Gsize] : Each of conductor Phi Profile, This is Device value
    // Lap_SIG_Sol[CondNUMR][CondNUMR] : Each of conductor Sigma Profile for external circuit

    //Make a Field DATA set  
    //Host_PCG_DATA = (DPS_Data*)malloc(A_size*sizeof(DPS_Data));
    //checkCudaErrors(cudaMalloc((void**)&dev_PCG_DATA, A_size*sizeof(DPS_Data)));
    //Make_PCG_DATA_Init<<<A_size/4,4>>>(dev_PCG_DATA,A_size,dev_M);
    //checkCudaErrors(cudaMemcpy(Host_PCG_DATA, dev_PCG_DATA, A_size*sizeof(DPS_Data), cudaMemcpyDeviceToHost));
    // Laplace Solution
    //cudaMallocPitch(&Lap_PHI_Sol, &pitch, Gsize * sizeof(float), CondNUMR); // for Laplace Solution
    //cudaMalloc((void**) &Lap_TEMP_Sol, Gsize * sizeof(int));
    // cudaMemset((void *) array, 0, Gsize * sizeof(int));
    

           // Find good grid and block size
        //cudaOccupancyMaxPotentialBlockSize(&mingrid,&FIELD_BLOCK,(void*)cusparseSpMV,0,Gsize); 
        //FIELD_GRID = (Gsize + FIELD_BLOCK - 1) / FIELD_BLOCK;
        //printf("blockSize = %d\n",FIELD_BLOCK);
        //printf("gridSize = %d\n",FIELD_GRID);
    // For test

    //if (matA       ) { checkCudaErrors(cusparseDestroySpMat(matA)); }
    //if (vecx       ) { checkCudaErrors(cusparseDestroyDnVec(vecx)); }
    //if (vecAP      ) { checkCudaErrors(cusparseDestroyDnVec(vecAP)); }
    //if (vecp       ) { checkCudaErrors(cusparseDestroyDnVec(vecp)); }
}
void Set_MatrixPCG_cuda(){

}
void Field_Method0_Initial(){
    // CPU Conjugate Gradient
    printf(" Field Solver : [CPU] Conjugate Gradient\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    X = VFMalloc(A_size);
    B = VFMalloc(A_size);
    R0 = VFMalloc(A_size);
    AX = VFMalloc(A_size);
    P0 = VFMalloc(A_size);
    AP = VFMalloc(A_size);
    PAP = VFMalloc(A_size);
    VFInit(B,0.0,A_size);
    VFInit(AX,0.0,A_size);
    VFInit(R0,0.0,A_size);
    VFInit(P0,0.0,A_size);
    VFInit(AP,0.0,A_size);
    VFInit(PAP,0.0,A_size);
}
void Field_Method1_Initial(){
    // CPU Preconditioned Conjugate Gradient [Jacovi]
    printf(" Field Solver : [CPU] Preconditioned Conjugate Gradient\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Preconditioner[Jacovi]\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    X = VFMalloc(A_size);
    B = VFMalloc(A_size);
    R0 = VFMalloc(A_size);
    AX = VFMalloc(A_size);
    Z0 = VFMalloc(A_size);
    P0 = VFMalloc(A_size);
    AP = VFMalloc(A_size);
    PAP = VFMalloc(A_size);
    VFInit(B,0.0,A_size);
    VFInit(AX,0.0,A_size);
    VFInit(R0,0.0,A_size);
    VFInit(Z0,0.0,A_size);
    VFInit(P0,0.0,A_size);
    VFInit(AP,0.0,A_size);
    VFInit(PAP,0.0,A_size);
}
void Field_Method2_Initial(){
    printf(" Field Solver : [GPU] Conjugate Gradient\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    // Cuda Handle setting
    cublasHandle = 0;
    cublasStatus = cublasCreate(&cublasHandle);
    checkCudaErrors(cublasStatus);
    cusparseHandle = 0;
    checkCudaErrors(cusparseCreate(&cusparseHandle));
    // Data cpu > gpu
    checkCudaErrors(cudaMalloc((void**) &dev_A, 5 * A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_Aj, 5 * A_size * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_Ai, (A_size + 1) * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_b,  A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_X,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_AP,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_R,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_P,  A_size * sizeof(float)));
    checkCudaErrors(cudaMemcpy(dev_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Aj, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Ai, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice));
    /* Wrap raw data into cuSPARSE generic API objects */
    matA = NULL;
    vecx = NULL;
    vecp = NULL;
    vecAP = NULL;
    checkCudaErrors(cusparseCreateCsr(
        &matA, A_size, A_size, 5*A_size, dev_Ai, dev_Aj, dev_A, CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ONE, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecx, A_size, dev_X, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecp, A_size, dev_P, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecAP, A_size, dev_AP, CUDA_R_32F));
}
void Field_Method3_Initial(){
    printf(" Field Solver : [GPU] Conjugate Gradient\n"); 
    printf(" Cuda Function : Cuda Graphs launch\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);

    // Data cpu > gpu
    checkCudaErrors(cudaMalloc((void**) &dev_A, 5 * A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_Aj, 5 * A_size * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_Ai, (A_size + 1) * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_b,  A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_X,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_AP,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_R,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_P,  A_size * sizeof(float)));
    checkCudaErrors(cudaMemcpy(dev_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Aj, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Ai, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice));
    /* Wrap raw data into cuSPARSE generic API objects */
    matA = NULL;
    vecx = NULL;
    vecp = NULL;
    vecAP = NULL;
    checkCudaErrors(cusparseCreateCsr(
        &matA, A_size, A_size, 5*A_size, dev_Ai, dev_Aj, dev_A, CUSPARSE_INDEX_32I,
        CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ONE, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecx, A_size, dev_X, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecp, A_size, dev_P, CUDA_R_32F));
    checkCudaErrors(cusparseCreateDnVec(&vecAP, A_size, dev_AP, CUDA_R_32F));
}
void Field_Method4_Initial(){
    printf(" Field Solver : [GPU] Conjugate Gradient\n"); 
    printf(" Cuda Function : Multi Block\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    // Data cpu > gpu
    checkCudaErrors(cudaMalloc((void**) &dev_A, 5 * A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_Aj, 5 * A_size * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_Ai, (A_size + 1) * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_b,  A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_X,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_AP,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_R,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_P,  A_size * sizeof(float)));
    checkCudaErrors(cudaMemcpy(dev_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Aj, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Ai, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice));
}
void Field_Method5_Initial(){
    printf(" Field Solver : [GPU] Preconditioned Conjugate Gradient\n"); 
    printf(" Cuda Function : Multi Block\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Preconditioner[Jacovi]\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    // Data cpu > gpu
    checkCudaErrors(cudaMalloc((void**) &dev_A, 5 * A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_Aj, 5 * A_size * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_Ai, (A_size + 1) * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_b,  A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_X,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_AP,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_R,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_P,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_M,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_Z,  A_size * sizeof(float)));
    checkCudaErrors(cudaMemcpy(dev_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Aj, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Ai, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(dev_M, MatM, A_size * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemset((void *) dev_Z, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_P, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_AP, 0, A_size * sizeof(float)));
}
void Field_Method6_Initial(){
    printf(" Field Solver : [GPU] Preconditioned Conjugate Gradient\n"); 
    printf(" Cuda Function : Multi GPU\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Preconditioner[Jacovi]\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    // Data cpu > gpu
    checkCudaErrors(cudaMallocManaged((void **)&man_I, sizeof(int) * (A_size + 1)));
    checkCudaErrors(cudaMallocManaged((void **)&man_J, sizeof(int) * 5 * A_size));
    checkCudaErrors(cudaMallocManaged((void **)&man_A, sizeof(float) * 5 * A_size));
    checkCudaErrors(cudaMemcpy(man_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(man_J, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(man_I, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemAdvise(man_I, sizeof(int) * (A_size + 1), cudaMemAdviseSetReadMostly, 0));
    checkCudaErrors(cudaMemAdvise(man_J, sizeof(int) * 5 * A_size, cudaMemAdviseSetReadMostly, 0));
    checkCudaErrors(cudaMemAdvise(man_A, sizeof(float) * 5 * A_size, cudaMemAdviseSetReadMostly, 0));
    // temp memory for ConjugateGradient
    checkCudaErrors(cudaMallocManaged((void **)&man_R, A_size * sizeof(float)));
    checkCudaErrors(cudaMallocManaged((void **)&man_P, A_size * sizeof(float)));
    checkCudaErrors(cudaMallocManaged((void **)&man_AP, A_size * sizeof(float)));
    checkCudaErrors(cudaMallocManaged((void **)&man_X, A_size * sizeof(float)));
    checkCudaErrors(cudaMallocManaged((void **)&man_Z, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) man_Z, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) man_P, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) man_AP, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) man_X, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) man_R, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMallocManaged((void **)&man_M, A_size * sizeof(float)));
    checkCudaErrors(cudaMemcpy(man_M, MatM, A_size * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemAdvise(man_M, A_size * sizeof(float), cudaMemAdviseSetReadMostly, 0));
}
void Field_Method7_Initial(){
    printf(" Field Solver : [GPU] Preconditioned Conjugate Gradient\n"); 
    printf(" Cuda Function : Multi Block\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Preconditioner[IChol]\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    // Data cpu > gpu
    checkCudaErrors(cudaMalloc((void**) &dev_A, 5 * A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_Aj, 5 * A_size * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_Ai, (A_size + 1) * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_b,  A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_X,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_AP,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_R,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_P,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_M,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_Z,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_Y,  A_size * sizeof(float)));
    checkCudaErrors(cudaMemcpy(dev_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Aj, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Ai, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(dev_M, MatM, A_size * sizeof(float), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemset((void *) dev_Z, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_P, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_AP, 0, A_size * sizeof(float)));
    checkCudaErrors(cudaMemset((void *) dev_Y, 0, A_size * sizeof(float)));
}
void PCG_Laplace_TEST(){
    // Solve Laplace Equation. 
    // INPUT
    // Lap_Field_Solver_Flag 0 - [CPU] Conjugate Gradient 
    // Lap_Field_Solver_Flag 1 - [CPU] Preconditioned Conjugate Gradient 
    // Lap_Field_Solver_Flag 2 - [GPU] Conjugate Gradient 
    // Lap_Field_Solver_Flag 3 - [GPU] Conjugate Gradient + Cuda Graphs launch
    // Lap_Field_Solver_Flag 4 - [GPU] Conjugate Gradient + Multi Block
    // Lap_Field_Solver_Flag 5 - [GPU] [Jacovi] Preconditioned Conjugate Gradient + Multi Block
    // Lap_Field_Solver_Flag 6 - [GPU] [Jacovi] Preconditioned Conjugate Gradient + Multi GPU 
    // Lap_Field_Solver_Flag 7 - [GPU] [IChol or ILU] Preconditioned Conjugate Gradient + Multi Block 
    // OUTPUT
    // Lap_TEMP_Sol[Gsize] : Temperature Profile
    // Lap_PHI_Sol[CondNUMR][Gsize] : Each of conductor Phi Profile, This is Device value
    // Lap_SIG_Sol[CondNUMR][CondNUMR] : Each of conductor Sigma Profile for external circuit
    int i,j,k,TID; 
    int mingrid;
    int IIter;
    char Namebuf[256];
    float *buf;
    float **CPUsol;
    buf = VFMalloc(A_size);
    CPUsol = MFMalloc(CondNUMR,Gsize);
    vec_A_idx = (int *) malloc(ngx * ngy * sizeof(int));
    for (i = 0; i < ngx; i++) {
		for (j = 0; j < ngy; j++) {
			vec_A_idx[j + i * ngy] = A_idx[i][j];
		} 
	}
    // SPEED TEST
    cudaEvent_t start, stop;
    float gputime;
    //
    if(Lap_Field_Solver_Flag == 0){// [CPU] Conjugate Gradient 
        Field_Method0_Initial(); // Initial Setting
        for (k = 0; k < CondNUMR; k++) {
            VFCopy(B,cond_b[k],A_size);
            VFInit(X,0.0,A_size);
            cudaEventCreate(&start); cudaEventCreate(&stop);
	        cudaEventRecord( start, 0 );
            FieldIter = CG_CPU();
            cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	        cudaEventElapsedTime( &gputime, start, stop );
	        cudaEventDestroy( start );cudaEventDestroy( stop );
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d, time = %2.8f (ms)\n",FieldIter,gputime);
            // Make a Solution
            VFInit(CPUsol[k],0.0,Gsize);
            for(j=ngy-1;j>=0;j--){
                for(i=0;i<ngx;i++){
                    TID = i*ngy+j;
                    if((vec_G[TID].CondID-1)==k){
                        CPUsol[k][TID] = 1.0;
                    }
                    if(vec_A_idx[TID]){
                        CPUsol[k][TID] = X[vec_A_idx[TID]-1];
                    }
                }
            }
        }
        sprintf(Namebuf,"CPU_CG");
        if(Lap_Field_Solver_Save) Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Lap_Field_Solver_Flag == 1){// [CPU] Preconditioned Conjugate Gradient 
        Field_Method1_Initial(); // Initial Setting
        for (k = 0; k < CondNUMR; k++) {
            VFCopy(B,cond_b[k],A_size);
            VFInit(X,0.0,A_size);
            cudaEventCreate(&start); cudaEventCreate(&stop);
	        cudaEventRecord( start, 0 );
            FieldIter = PCG_CPU();
            cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	        cudaEventElapsedTime( &gputime, start, stop );
	        cudaEventDestroy( start );cudaEventDestroy( stop );
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d, time = %2.8f (ms)\n",FieldIter,gputime);
            // Make a Solution
            VFInit(CPUsol[k],0.0,Gsize);
            for(j=ngy-1;j>=0;j--){
                for(i=0;i<ngx;i++){
                    TID = i*ngy+j;
                    if((vec_G[TID].CondID-1)==k){
                        CPUsol[k][TID] = 1.0;
                    }
                    if(vec_A_idx[TID]){
                        CPUsol[k][TID] = X[vec_A_idx[TID]-1];
                    }
                }
            }
        }
        sprintf(Namebuf,"CPU_PCG");
        if(Lap_Field_Solver_Save) Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Lap_Field_Solver_Flag == 2){// [GPU] Conjugate Gradient 
		Field_Method2_Initial(); // Initial Setting 
        for (k = 0; k < CondNUMR; k++) {
            checkCudaErrors(cudaMemcpy(dev_R, cond_b[k], A_size * sizeof(float),cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
            cudaEventCreate(&start); cudaEventCreate(&stop);
	        cudaEventRecord( start, 0 );
            FieldIter = CG_GPU();
            cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	        cudaEventElapsedTime( &gputime, start, stop );
	        cudaEventDestroy( start );cudaEventDestroy( stop );
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d, time = %2.8f (ms)\n",FieldIter,gputime);
            // Make a Solution
            VFInit(CPUsol[k],0.0,Gsize);
            checkCudaErrors(cudaMemcpy(buf, dev_X, A_size * sizeof(float),cudaMemcpyDeviceToHost));
            for(j=ngy-1;j>=0;j--){
                for(i=0;i<ngx;i++){
                    TID = i*ngy+j;
                    if((vec_G[TID].CondID-1)==k){
                        CPUsol[k][TID] = 1.0;
                    }
                    if(vec_A_idx[TID]){
                        CPUsol[k][TID] = buf[vec_A_idx[TID]-1];
                    }
                }
            }
        }
        sprintf(Namebuf,"GPU_CG");
        if(Lap_Field_Solver_Save) Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Lap_Field_Solver_Flag == 3){// [GPU] Conjugate Gradient + Cuda Graphs launch
		Field_Method3_Initial(); // Initial Setting
        // stream
        for (k = 0; k < CondNUMR; k++) {
            // Cuda Handle setting
            checkCudaErrors(cudaStreamCreate(&stream1));
            cublasHandle = 0;
            cublasStatus = cublasCreate(&cublasHandle);
            checkCudaErrors(cublasStatus);
            cusparseHandle = 0;
            checkCudaErrors(cusparseCreate(&cusparseHandle));
            cusparseStatus = cusparseCreate(&cusparseHandle);
            checkCudaErrors(cusparseStatus);
            checkCudaErrors(cudaMemcpy(dev_R, cond_b[k], A_size * sizeof(float),cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
            cudaEventCreate(&start); cudaEventCreate(&stop);
	        cudaEventRecord( start, 0 );
            FieldIter = CG_GPU_CudaGraphs();
            cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	        cudaEventElapsedTime( &gputime, start, stop );
	        cudaEventDestroy( start );cudaEventDestroy( stop );
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d, time = %2.8f (ms)\n",FieldIter,gputime);
            // Make a Solution
            VFInit(CPUsol[k],0.0,Gsize);
            checkCudaErrors(cudaMemcpyAsync(buf, dev_X, A_size * sizeof(float),cudaMemcpyDeviceToHost, streamForGraph));
            checkCudaErrors(cudaStreamSynchronize(streamForGraph));
            for(j=ngy-1;j>=0;j--){
                for(i=0;i<ngx;i++){
                    TID = i*ngy+j;
                    if((vec_G[TID].CondID-1)==k){
                        CPUsol[k][TID] = 1.0;
                    }
                    if(vec_A_idx[TID]){
                        CPUsol[k][TID] = buf[vec_A_idx[TID]-1];
                    }
                }
            }
            checkCudaErrors(cudaGraphExecDestroy(graphExec));
            checkCudaErrors(cudaGraphDestroy(initGraph));
            checkCudaErrors(cudaStreamDestroy(streamForGraph));
            checkCudaErrors(cudaStreamDestroy(stream1));
            checkCudaErrors(cusparseDestroy(cusparseHandle));
            checkCudaErrors(cublasDestroy(cublasHandle));
        }
        sprintf(Namebuf,"GPU_CG_Graph");
        if(Lap_Field_Solver_Save) Field_Laplace_Solution_Save(Namebuf,CPUsol);        
    }else if(Lap_Field_Solver_Flag == 4){// [GPU] Conjugate Gradient + Multi Block
        Field_Method4_Initial(); // Initial Setting
        //Make a Field constant set  
        Host_PCG_const = (DPS_Const*)malloc(sizeof(DPS_Const));
        checkCudaErrors(cudaMalloc((void**)&dev_PCG_const,sizeof(DPS_Const)));
        Make_PCG_Const_Init<<<1,1>>>(dev_PCG_const,A_size,PCGtol);
        checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));
        //checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));

        cudaDeviceProp deviceProp;
        int sMemSize = sizeof(double) * THREADS_PER_BLOCK;
        int numBlocksPerSm = 0;
        int numThreads = THREADS_PER_BLOCK;
        checkCudaErrors(cudaGetDeviceProperties(&deviceProp, device_num));
        if (!deviceProp.managedMemory) {
            // This sample requires being run on a device that supports Unified Memory
            fprintf(stderr, "Unified Memory not supported on this device\n");
            exit(EXIT_WAIVED);
        }
        // This sample requires being run on a device that supports Cooperative Kernel Launch
        if (!deviceProp.cooperativeLaunch)
        {
            printf("\nSelected GPU (%d) does not support Cooperative Kernel Launch, Waiving the run\n", device_num);
            exit(EXIT_WAIVED);
        }
        // Statistics about the GPU device
        printf("> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n\n",
           deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);
        
        checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, gpuConjugateGradient, numThreads, sMemSize));
        int numSms = deviceProp.multiProcessorCount;
        dim3 dimGrid(numSms*numBlocksPerSm, 1, 1), dimBlock(THREADS_PER_BLOCK, 1, 1);
        float nz = 5*A_size;
        //
        double *dot_result;
        cudaMallocManaged((void **)&dot_result, sizeof(double));
        *dot_result = 0.0;
        //
        checkCudaErrors(cudaMalloc((void**)&dot_result,sizeof(float)));
        void *kernelArgs[] = {
            (void*)&dev_Ai,
            (void*)&dev_Aj,
            (void*)&dev_A,
            (void*)&dev_X,
            (void*)&dev_AP,
            (void*)&dev_P,
            (void*)&dev_R,
            (void*)&dev_PCG_const,
            (void*)&dot_result,
        };
        for (k = 0; k < CondNUMR; k++) {
            checkCudaErrors(cudaMemcpy(dev_PCG_const, Host_PCG_const,sizeof(DPS_Const), cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemcpy(dev_b, cond_b[k], A_size * sizeof(float),cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemcpy(dev_R, dev_b, A_size * sizeof(float),cudaMemcpyDeviceToDevice));
            checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
            checkCudaErrors(cudaMemset((void *) dev_AP, 0, A_size * sizeof(float)));
            checkCudaErrors(cudaMemset((void *) dev_P, 0, A_size * sizeof(float)));
            cudaEventCreate(&start); cudaEventCreate(&stop);
	        cudaEventRecord( start, 0 );
            checkCudaErrors(cudaLaunchCooperativeKernel((void *)gpuConjugateGradient, dimGrid, dimBlock, kernelArgs, sMemSize, NULL));
            checkCudaErrors(cudaDeviceSynchronize());
            cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	        cudaEventElapsedTime( &gputime, start, stop );
	        cudaEventDestroy( start );cudaEventDestroy( stop );
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf(" time = %2.8f (ms)\n",gputime);
            // Make a Solution
            VFInit(CPUsol[k],0.0,Gsize);
            checkCudaErrors(cudaMemcpy(buf, dev_X, A_size * sizeof(float),cudaMemcpyDeviceToHost));
            for(j=ngy-1;j>=0;j--){
                for(i=0;i<ngx;i++){
                    TID = i*ngy+j;
                    if((vec_G[TID].CondID-1)==k){
                        CPUsol[k][TID] = 1.0;
                    }
                    if(vec_A_idx[TID]){
                        CPUsol[k][TID] = buf[vec_A_idx[TID]-1];
                    }
                }
            }
        }
        sprintf(Namebuf,"GPU_CG_MultiBlock");
        if(Lap_Field_Solver_Save) Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Lap_Field_Solver_Flag == 5){// [GPU] [Jacovi] Preconditioned Conjugate Gradient + Multi Block
		Field_Method5_Initial(); // Initial Setting
        //Make a Field constant set  
        Host_PCG_const = (DPS_Const*)malloc(sizeof(DPS_Const));
        checkCudaErrors(cudaMalloc((void**)&dev_PCG_const,sizeof(DPS_Const)));
        Make_PCG_Const_Init<<<1,1>>>(dev_PCG_const,A_size,PCGtol);
        checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));
        //checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));

        cudaDeviceProp deviceProp;
        int sMemSize = sizeof(double) * THREADS_PER_BLOCK;
        int numBlocksPerSm = 0;
        int numThreads = THREADS_PER_BLOCK;
        checkCudaErrors(cudaGetDeviceProperties(&deviceProp, device_num));
        if (!deviceProp.managedMemory) {
            // This sample requires being run on a device that supports Unified Memory
            fprintf(stderr, "Unified Memory not supported on this device\n");
            exit(EXIT_WAIVED);
        }
        // This sample requires being run on a device that supports Cooperative Kernel Launch
        if (!deviceProp.cooperativeLaunch)
        {
            printf("\nSelected GPU (%d) does not support Cooperative Kernel Launch, Waiving the run\n", device_num);
            exit(EXIT_WAIVED);
        }
        // Statistics about the GPU device
        printf("> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n\n",
           deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);
        
        checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, gpuPreConjugateGradient, numThreads, sMemSize));
        int numSms = deviceProp.multiProcessorCount;
        dim3 dimGrid(numSms*numBlocksPerSm, 1, 1), dimBlock(THREADS_PER_BLOCK, 1, 1);
        float nz = 5*A_size;
        //
        double *dot_result;
        cudaMallocManaged((void **)&dot_result, sizeof(double));
        *dot_result = 0.0;
        //
        void *kernelArgs[] = {
            (void*)&dev_Ai,
            (void*)&dev_Aj,
            (void*)&dev_A,
            (void*)&dev_M,
            (void*)&dev_X,
            (void*)&dev_AP,
            (void*)&dev_P,
            (void*)&dev_R,
            (void*)&dev_Z,
            (void*)&dev_PCG_const,
            (void*)&dot_result,
        };
        for (k = 0; k < CondNUMR; k++) {
            checkCudaErrors(cudaMemcpy(dev_PCG_const, Host_PCG_const,sizeof(DPS_Const), cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemcpy(dev_b, cond_b[k], A_size * sizeof(float),cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemcpy(dev_R, dev_b, A_size * sizeof(float),cudaMemcpyDeviceToDevice));
            checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
            checkCudaErrors(cudaMemset((void *) dev_AP, 0, A_size * sizeof(float)));
            checkCudaErrors(cudaMemset((void *) dev_P, 0, A_size * sizeof(float)));
            cudaEventCreate(&start); cudaEventCreate(&stop);
	        cudaEventRecord( start, 0 );
            checkCudaErrors(cudaLaunchCooperativeKernel((void *)gpuPreConjugateGradient, dimGrid, dimBlock, kernelArgs, sMemSize, NULL));
            checkCudaErrors(cudaDeviceSynchronize());
            cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	        cudaEventElapsedTime( &gputime, start, stop );
	        cudaEventDestroy( start );cudaEventDestroy( stop );
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf(" time = %2.8f (ms)\n",gputime);
            // Make a Solution
            VFInit(CPUsol[k],0.0,Gsize);
            checkCudaErrors(cudaMemcpy(buf, dev_X, A_size * sizeof(float),cudaMemcpyDeviceToHost));
            for(j=ngy-1;j>=0;j--){
                for(i=0;i<ngx;i++){
                    TID = i*ngy+j;
                    if((vec_G[TID].CondID-1)==k){
                        CPUsol[k][TID] = 1.0;
                    }
                    if(vec_A_idx[TID]){
                        CPUsol[k][TID] = buf[vec_A_idx[TID]-1];
                    }
                }
            }
        }
        sprintf(Namebuf,"GPU_PCG_Jacobi_MB");
        if(Lap_Field_Solver_Save) Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Lap_Field_Solver_Flag == 6){// [GPU] [Jacovi] Preconditioned Conjugate Gradient + Multi GPU
		Field_Method6_Initial(); // Initial Setting
        //Make a Field constant set  
        Host_PCG_const = (DPS_Const*)malloc(sizeof(DPS_Const));
        checkCudaErrors(cudaMalloc((void**)&dev_PCG_const,sizeof(DPS_Const)));
        Make_PCG_Const_Init<<<1,1>>>(dev_PCG_const,A_size,PCGtol);
        checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));
        //checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));

        cudaDeviceProp deviceProp;
        int num_of_gpus = 0;
        int num_buf = device_num;
        GPUn = 2; 
        deviceN = VIMalloc(GPUn);
        checkCudaErrors(cudaGetDeviceCount(&num_of_gpus));
        if (num_of_gpus <= 1 || num_of_gpus < device_num + GPUn) {
            printf("No. of GPU on node %d\n", num_of_gpus);
            printf("Minimum Two or more GPUs are required to run this code\n");
            exit(EXIT_WAIVED);
        }
        printf("Using GPU list : %d\n",GPUn);
        for(i=0;i<GPUn;i++){
            deviceN[i] = num_buf;
            num_buf++;
            cudaGetDeviceProperties(&deviceProp, deviceN[i]); 
            printf("Name %d : %s \n",deviceN[i], deviceProp.name);
            if (!deviceProp.managedMemory) {
            // This sample requires being run on a device that supports Unified Memory
                fprintf(stderr, "Unified Memory not supported on this device\n");
                exit(EXIT_WAIVED);
            }
            // This sample requires being run on a device that supports Cooperative Kernel Launch
            if (!deviceProp.cooperativeLaunch)
            {
                printf("\nSelected GPU (%d) does not support Cooperative Kernel Launch, Waiving the run\n", device_num);
                exit(EXIT_WAIVED);
            }
        }
        //
        double *dot_result;
        cudaMallocManaged((void **)&dot_result, sizeof(double));
        checkCudaErrors(cudaMemset(dot_result, 0.0, sizeof(double)));
        //
        cudaStream_t *nStreams = (cudaStream_t *)malloc(GPUn * sizeof(cudaStream_t));
        int NNZ = 5*A_size;
        void *kernelArgs[] = {
            (void*)&man_I,
            (void*)&man_J,
            (void*)&man_A,
            (void*)&man_M,
            (void*)&man_X,
            (void*)&man_AP,
            (void*)&man_P,
            (void*)&man_R,
            (void*)&man_Z,
            (void*)&NNZ,
            (void*)&A_size,
            (void*)&PCGtol,
            (void*)&dot_result,
        };
        int sMemSize = sizeof(double) * THREADS_PER_BLOCK;
        int numBlocksPerSm = 0;
        int numThreads = THREADS_PER_BLOCK;
        num_buf = device_num;
        checkCudaErrors(cudaSetDevice(num_buf));
        checkCudaErrors(cudaGetDeviceProperties(&deviceProp, num_buf));                    
        checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(
            &numBlocksPerSm, multiGpuPreConjugateGradient, numThreads, sMemSize));
        int numSms = deviceProp.multiProcessorCount;
        dim3 dimGrid(numSms * numBlocksPerSm, 1, 1), dimBlock(THREADS_PER_BLOCK, 1, 1);
        
        int device_count = 0;
        int totalThreadsPerGPU = numSms * numBlocksPerSm * THREADS_PER_BLOCK;
        num_buf = device_num;
        // Data Divide
        for(device_count = 0;device_count<GPUn;device_count++){
            num_buf = deviceN[device_count];
            checkCudaErrors(cudaSetDevice(num_buf));
            checkCudaErrors(cudaGetDeviceProperties(&deviceProp, num_buf));
            checkCudaErrors(cudaStreamCreate(&nStreams[device_count]));
            if (deviceProp.concurrentManagedAccess) {
                int perGPUIter = A_size / (totalThreadsPerGPU * GPUn);
                int offset_Ax = device_count * totalThreadsPerGPU;
                int offset_r = device_count * totalThreadsPerGPU;
                int offset_p = device_count * totalThreadsPerGPU;
                int offset_x = device_count * totalThreadsPerGPU;
                checkCudaErrors(cudaMemPrefetchAsync(man_I, sizeof(int) * (A_size+1), num_buf,nStreams[device_count]));
                checkCudaErrors(cudaMemPrefetchAsync(man_A, sizeof(float) * 5*A_size, num_buf,nStreams[device_count]));
                checkCudaErrors(cudaMemPrefetchAsync(man_J, sizeof(int) * 5*A_size, num_buf,nStreams[device_count]));
                if (offset_Ax <= A_size) {
                    for (i = 0; i < perGPUIter; i++) {
                        cudaMemAdvise(man_AP + offset_Ax, sizeof(float) * totalThreadsPerGPU,cudaMemAdviseSetPreferredLocation, num_buf);
                        cudaMemAdvise(man_R + offset_r, sizeof(float) * totalThreadsPerGPU,cudaMemAdviseSetPreferredLocation, num_buf);
                        cudaMemAdvise(man_X + offset_x, sizeof(float) * totalThreadsPerGPU,cudaMemAdviseSetPreferredLocation, num_buf);
                        cudaMemAdvise(man_P + offset_p, sizeof(float) * totalThreadsPerGPU,cudaMemAdviseSetPreferredLocation, num_buf);
                        cudaMemAdvise(man_AP + offset_Ax, sizeof(float) * totalThreadsPerGPU,cudaMemAdviseSetAccessedBy, num_buf);
                        cudaMemAdvise(man_R + offset_r, sizeof(float) * totalThreadsPerGPU,cudaMemAdviseSetAccessedBy, num_buf);
                        cudaMemAdvise(man_P + offset_p, sizeof(float) * totalThreadsPerGPU,cudaMemAdviseSetAccessedBy, num_buf);
                        cudaMemAdvise(man_X + offset_x, sizeof(float) * totalThreadsPerGPU,cudaMemAdviseSetAccessedBy, num_buf);
                        offset_Ax += totalThreadsPerGPU * GPUn;
                        offset_r += totalThreadsPerGPU * GPUn;
                        offset_p += totalThreadsPerGPU * GPUn;
                        offset_x += totalThreadsPerGPU * GPUn;
                        if (offset_Ax >= A_size) {
                            break;
                        }
                    }
                }
            }
        }
        printf("Total threads per GPU = %d numBlocksPerSm  = %d\n",numSms * numBlocksPerSm * THREADS_PER_BLOCK, numBlocksPerSm);
        launchParamsList = (cudaLaunchParams *)malloc(GPUn * sizeof(cudaLaunchParams));
        for (i = 0; i < GPUn; i++) {
            launchParamsList[i].func = (void *)multiGpuPreConjugateGradient;
            launchParamsList[i].gridDim = dimGrid;
            launchParamsList[i].blockDim = dimBlock;
            launchParamsList[i].sharedMem = sMemSize;
            launchParamsList[i].stream = nStreams[i];
            launchParamsList[i].args = kernelArgs;
        }
        for (k = 0; k < CondNUMR; k++) {
            checkCudaErrors(cudaMemcpy(dev_PCG_const, Host_PCG_const,sizeof(DPS_Const), cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemcpy(man_R, cond_b[k], A_size * sizeof(float),cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemset((void *) man_X, 0.0, A_size * sizeof(float)));
            checkCudaErrors(cudaMemset((void *) man_AP, 0.0, A_size * sizeof(float)));
            checkCudaErrors(cudaMemset((void *) man_P, 0.0, A_size * sizeof(float)));
            cudaEventCreate(&start); cudaEventCreate(&stop);
	        cudaEventRecord(start,0);
            checkCudaErrors(cudaLaunchCooperativeKernelMultiDevice(
                launchParamsList, GPUn,
                cudaCooperativeLaunchMultiDeviceNoPreSync |
                cudaCooperativeLaunchMultiDeviceNoPostSync));
            checkCudaErrors(cudaMemPrefetchAsync(man_X, sizeof(float) * A_size, cudaCpuDeviceId));
            checkCudaErrors(cudaMemPrefetchAsync(dot_result, sizeof(double), cudaCpuDeviceId));
            for(device_count = 0;device_count<GPUn;device_count++){
                num_buf = deviceN[device_count];
                checkCudaErrors(cudaSetDevice(num_buf));
                checkCudaErrors(cudaStreamSynchronize(nStreams[device_count]));
            }
            cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	        cudaEventElapsedTime( &gputime, start, stop );
	        cudaEventDestroy( start );cudaEventDestroy( stop );
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf(" time = %2.8f (ms)",gputime);
            float r1 = *dot_result;
            printf(", residual = %e \n", r1);

            // Make a Solution
            VFInit(CPUsol[k],0.0,Gsize);
            for(j=ngy-1;j>=0;j--){
                for(i=0;i<ngx;i++){
                    TID = i*ngy+j;
                    if((vec_G[TID].CondID-1)==k){
                        CPUsol[k][TID] = 1.0;
                    }
                    if(vec_A_idx[TID]){
                        CPUsol[k][TID] = man_X[vec_A_idx[TID]-1];
                        //if(CPUsol[k][TID]!=0) printf("CHECK %g \n  ", CPUsol[k][TID]);
                    }
                }
            }

        }
        sprintf(Namebuf,"GPU_PCG_MultiGPU");
        if(Lap_Field_Solver_Save) Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Lap_Field_Solver_Flag == 7){// [GPU] [IChol] Preconditioned Conjugate Gradient + Multi Block
		Field_Method7_Initial(); // Initial Setting
        //Make a Field constant set  
        Host_PCG_const = (DPS_Const*)malloc(sizeof(DPS_Const));
        checkCudaErrors(cudaMalloc((void**)&dev_PCG_const,sizeof(DPS_Const)));
        Make_PCG_Const_Init<<<1,1>>>(dev_PCG_const,A_size,PCGtol);
        checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));
        //checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));

        cudaDeviceProp deviceProp;
        int sMemSize = sizeof(double) * THREADS_PER_BLOCK;
        int numBlocksPerSm = 0;
        int numThreads = THREADS_PER_BLOCK;
        checkCudaErrors(cudaGetDeviceProperties(&deviceProp, device_num));
        if (!deviceProp.managedMemory) {
            // This sample requires being run on a device that supports Unified Memory
            fprintf(stderr, "Unified Memory not supported on this device\n");
            exit(EXIT_WAIVED);
        }
        // This sample requires being run on a device that supports Cooperative Kernel Launch
        if (!deviceProp.cooperativeLaunch)
        {
            printf("\nSelected GPU (%d) does not support Cooperative Kernel Launch, Waiving the run\n", device_num);
            exit(EXIT_WAIVED);
        }
        // Statistics about the GPU device
        printf("> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n",
           deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);
        
        checkCudaErrors(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&numBlocksPerSm, gpuPreConjugateGradient, numThreads, sMemSize));
        int numSms = deviceProp.multiProcessorCount;
        dim3 dimGrid(numSms*numBlocksPerSm, 1, 1), dimBlock(THREADS_PER_BLOCK, 1, 1);
        float nz = 5*A_size;
        //
        double *dot_result;
        cudaMallocManaged((void **)&dot_result, sizeof(double));
        *dot_result = 0.0;
        // Make a Preconditioner matrix
        checkCudaErrors(cudaMalloc((void**) &dev_L, 5 * A_size * sizeof(float)));
	    checkCudaErrors(cudaMalloc((void**) &dev_Lj, 5 * A_size * sizeof(int)));
	    checkCudaErrors(cudaMalloc((void**) &dev_Li, (A_size + 1) * sizeof(int)));
        checkCudaErrors(cudaMemset((void *) dev_L, 0, 5 * A_size * sizeof(float)));
        checkCudaErrors(cudaMemset((void *) dev_Lj, 0, 5 * A_size * sizeof(int)));
        checkCudaErrors(cudaMemset((void *) dev_Li, 0, (A_size + 1) * sizeof(int)));
        checkCudaErrors(cudaMalloc((void**) &dev_U, 5 * A_size * sizeof(float)));
	    checkCudaErrors(cudaMalloc((void**) &dev_Uj, 5 * A_size * sizeof(int)));
	    checkCudaErrors(cudaMalloc((void**) &dev_Ui, (A_size + 1) * sizeof(int)));
        checkCudaErrors(cudaMemset((void *) dev_U, 0, 5 * A_size * sizeof(float)));
        checkCudaErrors(cudaMemset((void *) dev_Uj, 0, 5 * A_size * sizeof(int)));
        checkCudaErrors(cudaMemset((void *) dev_Ui, 0, (A_size + 1) * sizeof(int)));
        if(Preconditioner_Flag==0){
            // Incomplete Cholesky Preconditioner
            printf("Make a preconditioner : [I Cholesky]\n");
            float *L_val;
            int *Li,Lj;
            int row_elem,next_row_elem;
            L_val = VFMalloc(5*A_size);VFInit(L_val,0.0,5*A_size);
            Li = VIMalloc(A_size+1);VIInit(Li,0,A_size+1);
            Lj = VIMalloc(5*A_size);VIInit(Lj,0,5*A_size);
            k = 0;
            for (i=0; i < A_size; i++){ //ROW
                row_elem = Ai[i];
                next_row_elem = Ai[i+1];
                for (j=row_elem-1; j < next_row_elem-1; j++){ // Column
                    if(A_val[j] != 0){
                        if(j=i){
                            L_val[k] = sqrt(A_val[j]);
                            k++;
                        }else if(j<i){
                            L_val[k] = sqrt(A_val[j]); // start
                            k++;
                        }
                    }
                }
            }
 
	        exit(1);
        }else if(Preconditioner_Flag==1){
            // Incomplete LU Preconditioner
            printf("Make a preconditioner : [I LU]\n");

        }else{
            //Jacovi

        }
        printf("Complete!!\n");
        //
        void *kernelArgs[] = {
            (void*)&dev_Ai,
            (void*)&dev_Aj,
            (void*)&dev_A,
            (void*)&dev_M,
            (void*)&dev_X,
            (void*)&dev_AP,
            (void*)&dev_P,
            (void*)&dev_R,
            (void*)&dev_Z,
            (void*)&dev_PCG_const,
            (void*)&dot_result,
        };
        for (k = 0; k < CondNUMR; k++) {
            checkCudaErrors(cudaMemcpy(dev_PCG_const, Host_PCG_const,sizeof(DPS_Const), cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemcpy(dev_b, cond_b[k], A_size * sizeof(float),cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemcpy(dev_R, dev_b, A_size * sizeof(float),cudaMemcpyDeviceToDevice));
            checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
            checkCudaErrors(cudaMemset((void *) dev_AP, 0, A_size * sizeof(float)));
            checkCudaErrors(cudaMemset((void *) dev_P, 0, A_size * sizeof(float)));
            cudaEventCreate(&start); cudaEventCreate(&stop);
	        cudaEventRecord( start, 0 );
            checkCudaErrors(cudaLaunchCooperativeKernel((void *)gpuPreConjugateGradient, dimGrid, dimBlock, kernelArgs, sMemSize, NULL));
            checkCudaErrors(cudaDeviceSynchronize());
            cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	        cudaEventElapsedTime( &gputime, start, stop );
	        cudaEventDestroy( start );cudaEventDestroy( stop );
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf(" time = %2.8f (ms)\n",gputime);
            // Make a Solution
            VFInit(CPUsol[k],0.0,Gsize);
            checkCudaErrors(cudaMemcpy(buf, dev_X, A_size * sizeof(float),cudaMemcpyDeviceToHost));
            for(j=ngy-1;j>=0;j--){
                for(i=0;i<ngx;i++){
                    TID = i*ngy+j;
                    if((vec_G[TID].CondID-1)==k){
                        CPUsol[k][TID] = 1.0;
                    }
                    if(vec_A_idx[TID]){
                        CPUsol[k][TID] = buf[vec_A_idx[TID]-1];
                    }
                }
            }
        }
        sprintf(Namebuf,"GPU_PCG_IChol_MB");
        if(Lap_Field_Solver_Save) Field_Laplace_Solution_Save(Namebuf,CPUsol);
    
    }else if(Lap_Field_Solver_Flag >= 9){
        printf("Empty Test Laplace Field Solver\n");
        exit(1);
    }
}
int CG_CPU(){
    int TID,i,Iter=0;
    float tol2;
    float rsold,rnew,Temp;
    float alpha,beta;

    rsold = 0;
    for(TID=0;TID<A_size;TID++){
        AX[TID] = 0;
        for(i=Ai[TID]-Ai[0];i<Ai[TID+1]-Ai[0];i++){
            AX[TID] += A_val[i]*X[Aj[i]-Ai[0]];
        }
        R0[TID] = B[TID] - AX[TID];
        P0[TID] = R0[TID];
        rsold += R0[TID]*R0[TID]; //AtomicAdd!!
    }
    printf("First:result[0].rsold = %g\n",rsold);
    tol2 = PCGtol*PCGtol;
    while(rsold>tol2){
        Iter++;
        Temp = 0.0;
        for(TID=0;TID<A_size;TID++){
            AP[TID] = 0;
            for(i=Ai[TID]-Ai[0];i<Ai[TID+1]-Ai[0];i++){
                AP[TID] += A_val[i]*P0[Aj[i]-Ai[0]];
            }
            //if(P0[TID]!=0) printf("p[%d] = %g\n",TID,P0[TID]);
            //printf("AP[%d] = %g\n",TID,P0[TID]);
            PAP[TID] = P0[TID] * AP[TID];
            Temp += PAP[TID]; //AtomicAdd!!
            
        }
        //printf("Temp = %g\n",Temp);
        //exit(1);
        alpha = (Temp)? rsold/Temp:0.0f ;
        for(TID=0;TID<A_size;TID++){
            X[TID] = X[TID] + alpha * P0[TID];
            R0[TID] = R0[TID] - alpha * AP[TID];
            rnew += R0[TID]*R0[TID];  //AtomicAdd!!
        }
        beta = (rsold) ? rnew/rsold: 0.0f;
        for(TID=0;TID<A_size;TID++){ 
            P0[TID] = R0[TID] + beta*P0[TID];
        }
        rsold = rnew;
        rnew = 0.0;
        //if(Iter<10) printf("Iter = %d, Temp = %g, alpha = %g, beta = %g, rsold = %g\n",Iter,Temp,alpha,beta,rsold);
    }
    return Iter;
}
int PCG_CPU(){
    int TID,i,Iter=0;
    float tol2;
    float rsold,rnew,Temp;
    float alpha,beta;

    rsold = 0;
    for(TID=0;TID<A_size;TID++){
        AX[TID] = 0;
        for(i=Ai[TID]-Ai[0];i<Ai[TID+1]-Ai[0];i++){
            AX[TID] += A_val[i]*X[Aj[i]-Ai[0]];
        }
        R0[TID] = B[TID] - AX[TID];
        Z0[TID] = MatM[TID] * R0[TID];
        P0[TID] = Z0[TID];
        rsold += R0[TID]*Z0[TID]; //AtomicAdd!!
    }
    while(rsold>PCGtol2){
        Iter++;
        Temp = 0.0;
        for(TID=0;TID<A_size;TID++){
            AP[TID] = 0;
            for(i=Ai[TID]-Ai[0];i<Ai[TID+1]-Ai[0];i++){
                AP[TID] += A_val[i]*P0[Aj[i]-Ai[0]];
            }
            //printf("AP[%d] = %g\n",TID,P0[TID]);
            PAP[TID] = P0[TID] * AP[TID];
            Temp += PAP[TID]; //AtomicAdd!!
        }
        alpha = (Temp)? rsold/Temp:0.0f ;
        for(TID=0;TID<A_size;TID++){
            X[TID] = X[TID] + alpha * P0[TID];
            R0[TID] = R0[TID] - alpha * AP[TID];
            Z0[TID] = MatM[TID] * R0[TID];
            rnew += R0[TID]*Z0[TID];  //AtomicAdd!!
        }
        beta = (rsold) ? rnew/rsold: 0.0f;
        for(TID=0;TID<A_size;TID++){ 
            P0[TID] = Z0[TID] + beta*P0[TID];
        }
        rsold = rnew;
        rnew = 0.0;
        //if(Iter<10) printf("Iter = %d, Temp = %g, alpha = %g, beta = %g, rsold = %g\n",Iter,Temp,alpha,beta,rsold);
    }
    return Iter;
}
int CG_GPU(){
    int iter;
    int max_iter = 1000000;
    float a = 1.0;
    float b = 0.0;
    float na = -1.0;
    float rnew = 0.;
    float rsold,Temp;
    float alpha,beta;
    float nalpha;
    /* Allocate workspace for cuSPARSE */
    size_t bufferSize = 0;
    checkCudaErrors(cusparseSpMV_bufferSize(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &a, matA, vecx,
        &b, vecAP, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, &bufferSize));
    void *buffer = NULL;
    checkCudaErrors(cudaMalloc(&buffer, bufferSize));
    /* Begin CG */
    checkCudaErrors(cusparseSpMV(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &a, matA, vecx,
        &b, vecAP, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, buffer));
    checkCudaErrors(cublasSaxpy(cublasHandle, A_size, &na, dev_AP, 1, dev_R, 1));
    checkCudaErrors(cublasSdot(cublasHandle, A_size, dev_R, 1, dev_R, 1, &rsold));
    iter = 1;
    while (rsold > PCGtol2 && iter <= max_iter){
         if (iter > 1)
        {
            beta = rnew / rsold;
            rsold = rnew;
            cublasSscal(cublasHandle, A_size, &beta, dev_P, 1);
            cublasSaxpy(cublasHandle, A_size, &a, dev_R, 1, dev_P, 1);
        }
        else
        {
            cublasScopy(cublasHandle, A_size, dev_R, 1, dev_P, 1);
        }
        checkCudaErrors(cusparseSpMV(
            cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &a, matA,
            vecp, &b, vecAP, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, buffer));
        cublasSdot(cublasHandle, A_size, dev_P, 1, dev_AP, 1, &Temp);
        alpha = (Temp) ? rsold / Temp: 0.0f;
        cublasSaxpy(cublasHandle, A_size, &alpha, dev_P, 1, dev_X, 1);
        nalpha = -alpha;
        cublasSaxpy(cublasHandle, A_size, &nalpha, dev_AP, 1, dev_R, 1);
        cublasSdot(cublasHandle, A_size, dev_R, 1, dev_R, 1, &rnew);
        cudaDeviceSynchronize();
        //if(iter<10) printf("Iter = %d, Temp = %g, alpha = %g, beta = %g, rsold = %g\n",iter,Temp,a,b,rsold);
        iter++;
    }
    return iter;

}
int CG_GPU_CudaGraphs(){
    int iter;
    int max_iter = 1000000;
    static int init_Flag = 0;
    float a = 1.0;
    float b = 0.0;
    float na = -1.0;
    float r1;
    float *d_r1, *d_r0, *d_dot, *d_a, *d_na, *d_b;
    // SPEED TEST
    cudaEvent_t start, stop;
    float gputime;
    //
    /* Allocate workspace for cuSPARSE */
    size_t bufferSize = 0;
    checkCudaErrors(cusparseSpMV_bufferSize(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &a, matA, vecx,
        &b, vecAP, CUDA_R_32F, CUSPARSE_CSRMV_ALG1, &bufferSize));
    void *buffer = NULL;
    checkCudaErrors(cudaMalloc(&buffer, bufferSize));
    // DEVICE variable
    checkCudaErrors(cudaMalloc((void **)&d_r1, sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_r0, sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_dot, sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_a, sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_na, sizeof(float)));
    checkCudaErrors(cudaMalloc((void **)&d_b, sizeof(float)));
    //
    checkCudaErrors(cusparseSetStream(cusparseHandle, stream1));
    checkCudaErrors(cusparseSpMV(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &a, matA, vecx,
        &b, vecAP, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, buffer));
    checkCudaErrors(cublasSetStream(cublasHandle, stream1));
    checkCudaErrors(cublasSaxpy(cublasHandle, A_size, &na, dev_AP, 1, dev_R, 1));
    checkCudaErrors(cublasSetPointerMode(cublasHandle, CUBLAS_POINTER_MODE_DEVICE));
    checkCudaErrors(cublasSdot(cublasHandle, A_size, dev_R, 1, dev_R, 1, d_r1));

    iter = 1;
    // First Iteration when iter=1 starts
    checkCudaErrors(cublasScopy(cublasHandle, A_size, dev_R, 1, dev_P, 1));
    checkCudaErrors(cusparseSpMV(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &a, matA, vecp,
        &b, vecAP, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, buffer));
    checkCudaErrors(cublasSdot(cublasHandle, A_size, dev_P, 1, dev_AP, 1, d_dot));
    r1_div_x<<<1, 1, 0, stream1>>>(d_r1, d_dot, d_a);
    checkCudaErrors(cublasSaxpy(cublasHandle, A_size, d_a, dev_P, 1, dev_X, 1));
    a_minus<<<1, 1, 0, stream1>>>(d_a, d_na);
    checkCudaErrors(cublasSaxpy(cublasHandle, A_size, d_na, dev_AP, 1, dev_R, 1));
    checkCudaErrors(cudaMemcpyAsync(d_r0, d_r1, sizeof(float),cudaMemcpyDeviceToDevice, stream1));
    checkCudaErrors(cublasSdot(cublasHandle, A_size, dev_R, 1, dev_R, 1, d_r1));
    checkCudaErrors(cudaMemcpyAsync(&r1, d_r1, sizeof(float),cudaMemcpyDeviceToHost, stream1));
    checkCudaErrors(cudaStreamSynchronize(stream1));
    printf("iteration = %3d, residual = %e\n", iter, sqrt(r1));
    // First Iteration when k=1 ends
    iter++;
    checkCudaErrors(cudaStreamCreate(&streamForGraph));
    checkCudaErrors(cublasSetStream(cublasHandle, stream1));
    checkCudaErrors(cusparseSetStream(cusparseHandle, stream1));
   // if(init_Flag == 0){
    // Capture start
    checkCudaErrors(cudaStreamBeginCapture(stream1, cudaStreamCaptureModeGlobal));
    r1_div_x<<<1, 1, 0, stream1>>>(d_r1, d_r0, d_b);
    cublasSetPointerMode(cublasHandle, CUBLAS_POINTER_MODE_DEVICE);
    checkCudaErrors(cublasSscal(cublasHandle, A_size, d_b, dev_P, 1));
    cublasSetPointerMode(cublasHandle, CUBLAS_POINTER_MODE_HOST);
    checkCudaErrors(cublasSaxpy(cublasHandle, A_size, &a, dev_R, 1, dev_P, 1));
    cublasSetPointerMode(cublasHandle, CUBLAS_POINTER_MODE_DEVICE);
    checkCudaErrors(cusparseSetPointerMode(cusparseHandle, CUSPARSE_POINTER_MODE_HOST));
    checkCudaErrors(cusparseSpMV(
        cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, &a, matA, vecp,
        &b, vecAP, CUDA_R_32F, CUSPARSE_MV_ALG_DEFAULT, buffer));
    checkCudaErrors(cudaMemsetAsync(d_dot, 0, sizeof(float), stream1));
    checkCudaErrors(cublasSdot(cublasHandle, A_size, dev_P, 1, dev_AP, 1, d_dot));
    r1_div_x<<<1, 1, 0, stream1>>>(d_r1, d_dot, d_a);
    checkCudaErrors(cublasSaxpy(cublasHandle, A_size, d_a, dev_P, 1, dev_X, 1));
    a_minus<<<1, 1, 0, stream1>>>(d_a, d_na);
    checkCudaErrors(cublasSaxpy(cublasHandle, A_size, d_na, dev_AP, 1, dev_R, 1));
    checkCudaErrors(cudaMemcpyAsync(d_r0, d_r1, sizeof(float),cudaMemcpyDeviceToDevice, stream1));
    checkCudaErrors(cudaMemsetAsync(d_r1, 0, sizeof(float), stream1));
    checkCudaErrors(cublasSdot(cublasHandle, A_size, dev_R, 1, dev_R, 1, d_r1));
    checkCudaErrors(cudaMemcpyAsync((float *)&r1, d_r1, sizeof(float),cudaMemcpyDeviceToHost, stream1));
    checkCudaErrors(cudaStreamEndCapture(stream1, &initGraph));
    // Capture End
    checkCudaErrors(cudaGraphInstantiate(&graphExec, initGraph, NULL, NULL, 0));
    checkCudaErrors(cublasSetStream(cublasHandle, stream1));
    checkCudaErrors(cusparseSetStream(cusparseHandle, stream1));
        init_Flag++;
    //}

    cudaEventCreate(&start); cudaEventCreate(&stop);
	cudaEventRecord( start, 0 );
    while (r1 > PCGtol2 && iter <= max_iter) {
        checkCudaErrors(cudaGraphLaunch(graphExec, streamForGraph));
        checkCudaErrors(cudaStreamSynchronize(streamForGraph));
        //printf("iteration = %3d, residual = %e\n", iter, sqrt(r1));
    iter++;
    }
    cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
    cudaEventElapsedTime( &gputime, start, stop );
	cudaEventDestroy( start );cudaEventDestroy( stop );
    printf("time = %2.8f (ms)\n",gputime);

    return iter;
}
__global__ void PCG_LAP(float *A,int *Ai,int *Aj,DPS_Const *PCG_C,DPS_Data *PCG_D,float *X,float *b){
    int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    if(TID>=PCG_C[0].A_size) return;
    int i;
    float sum;
    int MAXITER = 20;
    //Initial
    if(TID==0){
        PCG_C[0].Iter = 0;
        PCG_C[0].rsold = 0;
        PCG_C[0].Temp = 0;
        PCG_C[0].rnew = 0;
    } 
    __syncthreads();
    // cal  AP = A * P
    for(i=Ai[TID]-1;i<Ai[TID+1]-1;i++){
        PCG_D[TID].AP += A[i] * X[Aj[i]-1];
    }
    PCG_D[TID].R = b[TID] - PCG_D[TID].AP;
    PCG_D[TID].Z = PCG_D[TID].M * PCG_D[TID].R;
    PCG_D[TID].P = PCG_D[TID].Z;
    sum = PCG_D[TID].R*PCG_D[TID].Z;
    atomicAdd(&PCG_C[0].rsold,sum);
    //if(TID==0) printf("maxNorm [%d]= %g\n",TID,PCG_C[0].rsold);  
    //if(TID==10) printf("maxNorm [%d]= %g\n",TID,PCG_C[0].rsold);  
    //if(TID==0) printf(" [%d]Initial rsold = %g, tol2 = %g\n",TID,PCG_C[0].rsold,PCG_C[0].tol2);
    while(PCG_C[0].rsold > PCG_C[0].tol2){
        if(TID==0){
            PCG_C[0].Iter++;
            PCG_C[0].Temp = 0;
            PCG_C[0].rnew = 0;
        } 
        __syncthreads();
        //if(TID==0) printf(" [%d]Iter %d start!\n",TID,PCG_C[0].Iter);
        PCG_D[TID].AP = 0;
        for(i=Ai[TID]-1;i<Ai[TID+1]-1;i++){
            PCG_D[TID].AP += A[i] * PCG_D[Aj[i]-1].P;
        }
        //printf("[%d] AP = %g\n",TID,PCG_D[TID].P);
        sum = PCG_D[TID].P * PCG_D[TID].AP;
        atomicAdd(&PCG_C[0].Temp,sum);
        //if(TID==0) printf(" [%d]Temp = %g\n",TID,PCG_C[0].Temp);       
        if(PCG_C[0].Iter>MAXITER){
            //if(TID==0) printf("Iteration[%d] = %g\n",TID,PCG_C[0].Iter);    
            break;
        } 
        if(TID==0){
            PCG_C[0].alpha = (PCG_C[0].Temp)? PCG_C[0].rsold/PCG_C[0].Temp : 0.0f; 
        } 
        __syncthreads();     
        X[TID] = X[TID] + PCG_C[0].alpha * PCG_D[TID].P;
        PCG_D[TID].R = PCG_D[TID].R - PCG_C[0].alpha * PCG_D[TID].AP;
        PCG_D[TID].Z = PCG_D[TID].M * PCG_D[TID].R;
        sum = PCG_D[TID].R * PCG_D[TID].Z;
        atomicAdd(&PCG_C[0].rnew,sum);
       //if(TID==0) printf(" [%d]rnew = %g\n",TID,PCG_C[0].rnew);     
        if(TID==0){
             PCG_C[0].beta = (PCG_C[0].rsold)? PCG_C[0].rnew/PCG_C[0].rsold : 0.0f;
        }
        __syncthreads();
        PCG_D[TID].P = PCG_D[TID].Z + PCG_C[0].beta*PCG_D[TID].P;    
        if(TID==0) PCG_C[0].rsold = PCG_C[0].rnew;
        __syncthreads();
        if(TID==0) printf(" [%d]Iter %d,Temp = %g,alpha = %g,beat = %g, rsold = %g\n",TID,PCG_C[0].Iter,PCG_C[0].Temp,PCG_C[0].alpha,PCG_C[0].beta,PCG_C[0].rsold);
    }
    if(TID==0) printf("Iter [%d]= %d, ",TID,PCG_C[0].Iter);  
    if(TID==0) printf("maxNorm [%d]= %g\n",TID,PCG_C[0].rsold);  
}
__global__ void Make_PCG_DATA_Init(DPS_Data *p, int size,float *MatrixM){
    int TID = blockIdx.x * blockDim.x + threadIdx.x;
    if(TID>=size) return;
    //printf("TID = %d, M = %g\n",TID,MatrixM[TID]);
    p[TID].R = 0.0;
    p[TID].Z = 0.0;
    p[TID].P = 0.0;
    p[TID].AP = 0.0;
    p[TID].M = MatrixM[TID];
}
__global__ void Make_PCG_Const_Init(DPS_Const *p,int Asize, float tol){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    p[i].A_size = Asize;
    p[i].Iter = 0;
    p[i].tol = tol;
    p[i].tol2 = tol*tol;
    p[i].rsold = 0.0;
    p[i].Temp = 0.0;
    p[i].rnew = 0.0;
    p[i].alpha = 0.0;
    p[i].beta = 0.0;    
}
__global__ void initVectors(float *rhs, float *x, int N) {
  size_t gid = blockIdx.x * blockDim.x + threadIdx.x;
  for (size_t i = gid; i < N; i += gridDim.x * blockDim.x) {
    rhs[i] = 1.0;
    x[i] = 0.0;
  }
}
__global__ void r1_div_x(float *r1, float *r0, float *b) {
  int gid = blockIdx.x * blockDim.x + threadIdx.x;
  if (gid == 0) {
    b[0] = r1[0] / r0[0];
  }
}
__global__ void a_minus(float *a, float *na) {
  int gid = blockIdx.x * blockDim.x + threadIdx.x;
  if (gid == 0) {
    na[0] = -(a[0]);
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

	float *row=(float *)((char *)A+height*pitch);

	PHI[TID]=row[TID];
}
__device__ void gpuSpMV(int *I, int *J, float *val, int nnz, int num_rows, float alpha, float *inputVecX, 
                        float *outputVecY, cg::thread_block &cta, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < num_rows; i+= grid.size())    {
        // i = 0 ~ A_size-1; 
        //printf("val[%d][]\n",i);
        //for(i=Ai[TID]-Ai[0];i<Ai[TID+1]-Ai[0];i++){
        //    AX[TID] += A_val[i]*X[Aj[i]-Ai[0]];
        //}
        int row_elem = I[i];
        int next_row_elem = I[i+1];
        int num_elems_this_row = next_row_elem - row_elem;
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
__device__ void gpuSaxpy(float *x, float *y, float a, int size, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < size; i+= grid.size()){        
        y[i] = a*x[i] + y[i];
    }
}
__device__ void gpuRSaxpy(float *x, float *y, float a, int size, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < size; i+= grid.size()){        
        y[i] = a*y[i] + x[i];
    }
}
__device__ void gpuDotProduct(float *vecA, float *vecB, double *result, int size, const cg::thread_block &cta, const cg::grid_group &grid)
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
__device__ void gpuCopyVector(float *srcA, float *destB, int size, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < size; i+= grid.size()){
        destB[i] = srcA[i];
    }
}
__device__ void gpuScaleVector(float *vec, float alpha, int size, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < size; i+= grid.size()){
        vec[i] = alpha*vec[i];
    }
}
__global__ void gpuConjugateGradient(int *I, int *J, float *val, float *x,  float *Ax, float *p, float *r, 
            DPS_Const *result,double *d_result)
{
    cg::thread_block cta = cg::this_thread_block();
    cg::grid_group grid = cg::this_grid();
    //int TID = blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    int k = 0;
    int max_iter = 10000;
    float a = 1.0;
    float na = -1.0;
    int nnz = 5 * result[0].A_size;
    int N = result[0].A_size;
    float rsold,rnew,Temp;
    float nalpha,alpha,beta;

    rsold = 0.0;
    if (threadIdx.x == 0 && blockIdx.x == 0){
        *d_result = 0.0;  
    } 
    gpuSpMV(I, J, val, nnz, N, a, x, Ax, cta, grid); 
    gpuSaxpy(Ax, r, na, N, grid); 
    gpuCopyVector(r, p, N, grid);
    //if(r[TID] !=0) printf("r[%d] = %g\n",TID,r[TID]);
    cg::sync(grid);
    gpuDotProduct(r, r, d_result, N, cta, grid); 
    cg::sync(grid);
    rsold = *d_result;
    //if(threadIdx.x == 0 && blockIdx.x == 0) printf("First:result[0].rsold = %g\n",rsold);
    //return;
    while (rsold > result[0].tol2 && k <= max_iter){
        k++;
        gpuSpMV(I, J, val, nnz, N, a, p, Ax, cta, grid);
        if (threadIdx.x == 0 && blockIdx.x == 0){
            *d_result = 0.0;  
        } 
        cg::sync(grid);
        //if(Ax[TID] !=0) printf("Ax[%d] = %g\n",TID,Ax[TID]);
        gpuDotProduct(p, Ax, d_result, N, cta, grid);
        cg::sync(grid);
        Temp = *d_result;
        //if(threadIdx.x == 0 && blockIdx.x == 0) printf("Temp = %g\n",Temp);
        //return;
        alpha = (Temp)? rsold/Temp:0.0f;
        gpuSaxpy(p, x, alpha, N, grid);
        nalpha = -alpha;
        gpuSaxpy(Ax, r, nalpha, N, grid);
        if (threadIdx.x == 0 && blockIdx.x == 0){
            *d_result = 0.0;  
        } 
        cg::sync(grid);
        gpuDotProduct(r, r, d_result, N, cta, grid);
        cg::sync(grid);
        rnew = *d_result;
        beta = (rsold) ? rnew/rsold: 0.0f;
        gpuRSaxpy(r, p, beta, N, grid);
        rsold = rnew;
        rnew = 0.0;
        //if(threadIdx.x == 0 && blockIdx.x == 0 && k<20) printf("Iter = %d, temp = %g,  AL = %g, BE = %g Res = %g\n",k,Temp,alpha,beta,rsold);
    }
    if(threadIdx.x == 0 && blockIdx.x == 0 ) printf("End Iter = %d, Res = %g, b = %g, a = %g\n",k,Temp,alpha,beta,rsold);
}
__device__ void gpuProductVector(float *vecA, float *vecB, float *vecC, int size, const cg::thread_block &cta, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < size; i+=grid.size()){
        vecC[i] = (vecA[i] * vecB[i]);
    }
}
__global__ void gpuPreConjugateGradient(int *I, int *J, float *val, float *M, float *x,  float *Ax, float *p, float *r, float *Z, 
            DPS_Const *result,double *d_result){
    cg::thread_block cta = cg::this_thread_block();
    cg::grid_group grid = cg::this_grid();
    //int TID = blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    int k = 0;
    int max_iter = 100000;
    float a = 1.0;
    float na = -1.0;
    int nnz = 5 * result[0].A_size;
    int N = result[0].A_size;
    float rsold,rnew,Temp;
    float nalpha,alpha,beta;

    rsold = 0.0;
    if (threadIdx.x == 0 && blockIdx.x == 0){
        *d_result = 0.0;  
    } 
    gpuSpMV(I, J, val, nnz, N, a, x, Ax, cta, grid); 
    gpuSaxpy(Ax, r, na, N, grid); 
    gpuProductVector(M, r, Z, N, cta, grid);
    gpuCopyVector(Z, p, N, grid);
    //if(r[TID] !=0) printf("r[%d] = %g\n",TID,r[TID]);
    cg::sync(grid);
    gpuDotProduct(r, Z, d_result, N, cta, grid); 
    cg::sync(grid);
    rsold = *d_result;
    //if(threadIdx.x == 0 && blockIdx.x == 0) printf("First:result[0].rsold = %g\n",rsold);
    //return;
    while (rsold > result[0].tol2 && k <= max_iter){
        k++;
        gpuSpMV(I, J, val, nnz, N, a, p, Ax, cta, grid);
        if (threadIdx.x == 0 && blockIdx.x == 0){
            *d_result = 0.0;  
        } 
        cg::sync(grid);
        //if(Ax[TID] !=0) printf("Ax[%d] = %g\n",TID,Ax[TID]);
        gpuDotProduct(p, Ax, d_result, N, cta, grid);
        cg::sync(grid);
        Temp = *d_result;
        //if(threadIdx.x == 0 && blockIdx.x == 0) printf("Temp = %g\n",Temp);
        //return;
        alpha = (Temp)? rsold/Temp:0.0f;
        gpuSaxpy(p, x, alpha, N, grid);
        nalpha = -alpha;
        gpuSaxpy(Ax, r, nalpha, N, grid);
        gpuProductVector(M, r, Z, N, cta, grid);
        if (threadIdx.x == 0 && blockIdx.x == 0){
            *d_result = 0.0;  
        } 
        cg::sync(grid);
        gpuDotProduct(r, Z, d_result, N, cta, grid);
        cg::sync(grid);
        rnew = *d_result;
        beta = (rsold) ? rnew/rsold: 0.0f;
        gpuRSaxpy(Z, p, beta, N, grid);
        rsold = rnew;
        rnew = 0.0;
        //if(threadIdx.x == 0 && blockIdx.x == 0 && k<20) printf("Iter = %d, temp = %g,  AL = %g, BE = %g Res = %g\n",k,Temp,alpha,beta,rsold);
    }
    if(threadIdx.x == 0 && blockIdx.x == 0 ) printf("End Iter = %d, Res = %g, b = %g, a = %g\n",k,Temp,alpha,beta,rsold);
}
__global__ void multiGpuPreConjugateGradient(int *I, int *J, float *val, float *M, float *x,  float *Ax, float *p, float *r, float *Z, 
            int nnz, int N, float tol, double *d_result)
{
    cg::thread_block cta = cg::this_thread_block();
    cg::grid_group grid = cg::this_grid();
    cg::multi_grid_group multi_grid = cg::this_multi_grid();

    //int TID = blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    int k = 0;
    int max_iter = 100000;
    float a = 1.0;
    float na = -1.0;
    float rsold,rnew,Temp;
    float nalpha,alpha,beta;

    rsold = 0.0; 
    cg::sync(grid);
    MultigpuSpMV(I, J, val, nnz, N, a, x, Ax, cta, multi_grid); 
    cg::sync(grid);
    MultigpuSaxpy(Ax, r, na, N, multi_grid); 
    cg::sync(grid);
    MultigpuProductVector(M, r, Z, N, cta, multi_grid);
    cg::sync(grid);
    MultigpuCopyVector(Z, p, N, multi_grid);
    cg::sync(grid);
    MultigpuDotProduct(r, Z, N, cta, multi_grid); 
    cg::sync(grid);
    if (grid.thread_rank() == 0) {
        atomicAdd_system(d_result, grid_dot_result);
        grid_dot_result = 0.0;
    }
    cg::sync(multi_grid);
    rsold = *d_result;
    //if (threadIdx.x == 0 && grid.thread_rank() == 0) printf("start : rsold1 = %g\n",rsold);
    while (rsold > tol*tol && k <= max_iter){
        k++;
        cg::sync(multi_grid);
        MultigpuSpMV(I, J, val, nnz, N, a, p, Ax, cta, multi_grid);
        if (multi_grid.thread_rank() == 0) {
            setDotResultToZero(d_result);
        }   
        cg::sync(multi_grid);
        MultigpuDotProduct(p, Ax, N, cta, multi_grid);
        cg::sync(grid);
        if (grid.thread_rank() == 0) {
            atomicAdd_system(d_result, grid_dot_result);
            grid_dot_result = 0.0;
        }
        cg::sync(multi_grid);
        Temp = *d_result;
        //if (threadIdx.x == 0 && grid.thread_rank() == 0) printf("Iter = %d, Temp = %g,",k,Temp);
        alpha = (Temp)? rsold/Temp:0.0f;
        MultigpuSaxpy(p, x, alpha, N, multi_grid);
        nalpha = -alpha;
        MultigpuSaxpy(Ax, r, nalpha, N, multi_grid);
        MultigpuProductVector(M, r, Z, N, cta, multi_grid);
        cg::sync(multi_grid);
        if (multi_grid.thread_rank() == 0) {
            setDotResultToZero(d_result);
        }   
        cg::sync(multi_grid);
        MultigpuDotProduct(r, Z, N, cta, multi_grid);
        cg::sync(grid);
        if (grid.thread_rank() == 0) {
            atomicAdd_system(d_result, grid_dot_result);
            grid_dot_result = 0.0;
        }
        cg::sync(multi_grid);
        rnew = *d_result;
        //if (threadIdx.x == 0 && grid.thread_rank() == 0) printf("rnew = %g\n",rnew);
        beta = (rsold) ? rnew/rsold: 0.0f;
        MultigpuRSaxpy(Z, p, beta, N, multi_grid);
        rsold = rnew;
        rnew = 0.0;
    }
    //if(threadIdx.x == 0 && blockIdx.x == 0 ) printf("End Iter = %d, Res = %g, b = %g, a = %g\n",k,Temp,alpha,beta,rsold);
}
__device__ void MultigpuSpMV(int *I, int *J, float *val, int nnz, int num_rows, float alpha, float *inputVecX, 
                        float *outputVecY, cg::thread_block &cta, const cg::multi_grid_group &multi_grid) {
    for (int i = multi_grid.thread_rank(); i < num_rows; i += multi_grid.size()) {
        // i = 0 ~ A_size-1; 
        int row_elem = I[i];
        int next_row_elem = I[i+1];
        int num_elems_this_row = next_row_elem - row_elem;
        float output = 0.0;
        for (int j=row_elem-1; j < next_row_elem-1; j++){
            output +=  alpha*val[j] * inputVecX[J[j]-1];
            //if(i==num_rows-1) printf("val[%d][%d] = %g, %g, %g\n",j,J[j]-1,val[j],inputVecX[J[j]-1],output);
        }
        //printf("output[%d] = %g\n",i,output);
        outputVecY[i] = output;
    }
}
__device__ void MultigpuSaxpy(float *x, float *y, float a, int size, const cg::multi_grid_group &multi_grid) {
    for (int i = multi_grid.thread_rank(); i < size; i += multi_grid.size()) {  
        y[i] = a*x[i] + y[i];
    }
}
__device__ void MultigpuRSaxpy(float *x, float *y, float a, int size, const cg::multi_grid_group &multi_grid) {
    for (int i = multi_grid.thread_rank(); i < size; i += multi_grid.size()) {  
        y[i] = a*y[i] + x[i];
    }
}
__device__ void MultigpuDotProduct(float *vecA, float *vecB, int size, const cg::thread_block &cta, const cg::multi_grid_group &multi_grid) {
   __shared__ double tmp[THREADS_PER_BLOCK];
    double temp_sum = 0.0;
    for (int i = multi_grid.thread_rank(); i < size; i += multi_grid.size()) {
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
        atomicAdd(&grid_dot_result, beta);
    }
}
__device__ void MultigpuCopyVector(float *srcA, float *destB, int size, const cg::multi_grid_group &multi_grid) {
    for (int i = multi_grid.thread_rank(); i < size; i += multi_grid.size()) {
        destB[i] = srcA[i];
    }
}
__device__ void MultigpuScaleVector(float *vec, float alpha, int size, const cg::multi_grid_group &multi_grid) {
    for (int i = multi_grid.thread_rank(); i < size; i += multi_grid.size()) {
        vec[i] = alpha*vec[i];
    }
}
__device__ void MultigpuProductVector(float *vecA, float *vecB, float *vecC, int size, const cg::thread_block &cta, const cg::multi_grid_group &multi_grid) {
    for (int i = multi_grid.thread_rank(); i < size; i += multi_grid.size()) {
        vecC[i] = (vecA[i] * vecB[i]);
    }
}
__device__ void setDotResultToZero(double *dot_result) {
  unsigned long long int *address_as_ull = (unsigned long long int *)dot_result;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS_system(address_as_ull, assumed, 0);

  } while (assumed != old);
}