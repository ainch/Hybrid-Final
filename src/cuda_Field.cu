#include "cuda_Field.cuh"
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
    // Cuda Handle setting
    cublasHandle = 0;
    cublasStatus = cublasCreate(&cublasHandle);
    checkCudaErrors(cublasStatus);
    cusparseHandle = 0;
    checkCudaErrors(cusparseCreate(&cusparseHandle));
    cusparseStatus = cusparseCreate(&cusparseHandle);
    checkCudaErrors(cusparseStatus);
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
    // stream
    checkCudaErrors(cudaStreamCreate(&stream1));
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
    checkCudaErrors(cudaMalloc((void**) &dev_A, 5 * A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_Aj, 5 * A_size * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_Ai, (A_size + 1) * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_b,  A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_X,  A_size * sizeof(float)));
    // Initialize
    checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
    //Copy
    checkCudaErrors(cudaMemcpy(dev_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Aj, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Ai, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice));
}
void Field_Method5_Initial(){
    printf(" Field Solver : [GPU] Preconditioned Conjugate Gradient\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Preconditioner[Jacovi]\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    checkCudaErrors(cudaMalloc((void**) &dev_A, 5 * A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_Aj, 5 * A_size * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_Ai, (A_size + 1) * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_b,  A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_X,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_M,  A_size * sizeof(float)));
    // Initialize
    checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
    //Copy
    checkCudaErrors(cudaMemcpy(dev_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Aj, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Ai, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(dev_M, MatM, A_size * sizeof(int), cudaMemcpyHostToDevice));
}
void Field_Method6_Initial(){
    printf(" Field Solver : [GPU] Preconditioned Conjugate Gradient\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Preconditioner[ILU]\n"); 
    printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    checkCudaErrors(cudaMalloc((void**) &dev_A, 5 * A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_Aj, 5 * A_size * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_Ai, (A_size + 1) * sizeof(int)));
	checkCudaErrors(cudaMalloc((void**) &dev_b,  A_size * sizeof(float)));
	checkCudaErrors(cudaMalloc((void**) &dev_X,  A_size * sizeof(float)));
    checkCudaErrors(cudaMalloc((void**) &dev_M,  A_size * sizeof(float)));
    // Initialize
    checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
    //Copy
    checkCudaErrors(cudaMemcpy(dev_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Aj, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice));
	checkCudaErrors(cudaMemcpy(dev_Ai, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(dev_M, MatM, A_size * sizeof(int), cudaMemcpyHostToDevice));
}
void PCG_SOLVER_Laplace(){
    // Solve Laplace Equation. 
    // INPUT
    // Field_Solver_Flag 0 - [CPU] Conjugate Gradient 
    // Field_Solver_Flag 1 - [CPU] Preconditioned Conjugate Gradient 
    // Field_Solver_Flag 2 - [GPU] Conjugate Gradient 
    // Field_Solver_Flag 3 - [GPU] Conjugate Gradient + Cuda Graphs launch
    // Field_Solver_Flag 4 - [GPU] Conjugate Gradient + Multi Block
    // Field_Solver_Flag 5 - [GPU] [Jacovi] Preconditioned Conjugate Gradient 
    // Field_Solver_Flag 6 - [GPU] [ILU] Preconditioned Conjugate Gradient 
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
    if(Field_Solver_Flag == 0){// [CPU] Conjugate Gradient 
        Field_Method0_Initial(); // Initial Setting
        for (k = 0; k < CondNUMR; k++) {
            VFCopy(B,cond_b[k],A_size);
            VFInit(X,0.0,A_size);
            FieldIter = CG_CPU();
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d\n",FieldIter);
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
        Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Field_Solver_Flag == 1){// [CPU] Preconditioned Conjugate Gradient 
        Field_Method1_Initial(); // Initial Setting
        for (k = 0; k < CondNUMR; k++) {
            VFCopy(B,cond_b[k],A_size);
            VFInit(X,0.0,A_size);
            FieldIter = PCG_CPU();
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d\n",FieldIter);
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
        Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Field_Solver_Flag == 2){// [GPU] Conjugate Gradient 
		Field_Method2_Initial(); // Initial Setting 
        for (k = 0; k < CondNUMR; k++) {
            checkCudaErrors(cudaMemcpy(dev_R, cond_b[k], A_size * sizeof(float),cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
            FieldIter = CG_GPU();
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d\n",FieldIter);
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
        Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Field_Solver_Flag == 3){// [GPU] Conjugate Gradient + Cuda Graphs launch
		Field_Method3_Initial(); // Initial Setting
        for (k = 0; k < CondNUMR; k++) {
            checkCudaErrors(cudaMemcpy(dev_R, cond_b[k], A_size * sizeof(float),cudaMemcpyHostToDevice));
            checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
            FieldIter = CG_GPU_CudaGraphs();
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d\n",FieldIter);
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
        }
        sprintf(Namebuf,"GPU_CG_Graph");
        Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Field_Solver_Flag == 4){// [GPU] Conjugate Gradient + Multi Block
		Field_Method4_Initial(); // Initial Setting
        for (k = 0; k < CondNUMR; k++) {

            
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d\n",FieldIter);
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
        Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Field_Solver_Flag == 5){// [GPU] [Jacovi] Preconditioned Conjugate Gradient
		Field_Method5_Initial(); // Initial Setting
        for (k = 0; k < CondNUMR; k++) {

            
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d\n",FieldIter);
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
        sprintf(Namebuf,"GPU_PCG_jacovi");
        Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else if(Field_Solver_Flag == 6){// [GPU] [ILU] Preconditioned Conjugate Gradient
		Field_Method6_Initial(); // Initial Setting
        for (k = 0; k < CondNUMR; k++) {

            
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d\n",FieldIter);
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
        sprintf(Namebuf,"GPU_PCG_ILU");
        Field_Laplace_Solution_Save(Namebuf,CPUsol);
    }else{

    }
    
    // Laplace Solution
    cudaMallocPitch(&Lap_PHI_Sol, &pitch, Gsize * sizeof(float), CondNUMR); // for Laplace Solution
    //cudaMalloc((void**) &Lap_TEMP_Sol, Gsize * sizeof(int));
    // cudaMemset((void *) array, 0, Gsize * sizeof(int));
    
    //Make a Field constant set  
    Host_PCG_const = (DPS_Const*)malloc(sizeof(DPS_Const));
    checkCudaErrors(cudaMalloc((void**)&dev_PCG_const,sizeof(DPS_Const)));
    Make_PCG_Const_Init<<<1,1>>>(dev_PCG_const,A_size,PCGtol);
    checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));
    //checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));
    //Make a Field DATA set  
    Host_PCG_DATA = (DPS_Data*)malloc(A_size*sizeof(DPS_Data));
    checkCudaErrors(cudaMalloc((void**)&dev_PCG_DATA, A_size*sizeof(DPS_Data)));
    Make_PCG_DATA_Init<<<A_size/4,4>>>(dev_PCG_DATA,A_size,dev_M);
    //checkCudaErrors(cudaMemcpy(Host_PCG_DATA, dev_PCG_DATA, A_size*sizeof(DPS_Data), cudaMemcpyDeviceToHost));
           // Find good grid and block size
        //cudaOccupancyMaxPotentialBlockSize(&mingrid,&FIELD_BLOCK,(void*)cusparseSpMV,0,Gsize); 
        //FIELD_GRID = (Gsize + FIELD_BLOCK - 1) / FIELD_BLOCK;
        //printf("blockSize = %d\n",FIELD_BLOCK);
        //printf("gridSize = %d\n",FIELD_GRID);
    // For test
    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);
    if (matA       ) { checkCudaErrors(cusparseDestroySpMat(matA)); }
    if (vecx       ) { checkCudaErrors(cusparseDestroyDnVec(vecx)); }
    if (vecAP      ) { checkCudaErrors(cusparseDestroyDnVec(vecAP)); }
    if (vecp       ) { checkCudaErrors(cusparseDestroyDnVec(vecp)); }
    exit(1);

}
void Set_MatrixPCG_cuda(){

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
    tol2 = PCGtol*PCGtol;
    while(rsold>tol2){
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
    int max_iter = 10000;
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
        if(iter<10) printf("Iter = %d, Temp = %g, alpha = %g, beta = %g, rsold = %g\n",iter,Temp,a,b,rsold);
        iter++;
    }
    return iter;

}
int CG_GPU_CudaGraphs(){
    int iter;
    int max_iter = 10000;
    float a = 1.0;
    float b = 0.0;
    float na = -1.0;
    float r1;
    float *d_r1, *d_r0, *d_dot, *d_a, *d_na, *d_b;
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
    while (r1 > PCGtol2 && iter <= max_iter) {
        checkCudaErrors(cudaGraphLaunch(graphExec, streamForGraph));
        checkCudaErrors(cudaStreamSynchronize(streamForGraph));
        if(iter<10) printf("iteration = %3d, residual = %e\n", iter, sqrt(r1));
    iter++;
    }
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
