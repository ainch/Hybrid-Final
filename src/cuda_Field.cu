#include "cuda_Field.cuh"

void PCG_SOLVER_Laplace(){
    // Solve Laplace Equation. (To use every time step.)
    // Goal
    // Lap_TEMP_Sol[Gsize] : Temperature Profile
    // Lap_PHI_Sol[CondNUMR][Gsize] : Each of conductor Phi Profile, This is Device value
    // Lap_SIG_Sol[CondNUMR][CondNUMR] : Each of conductor Sigma Profile for external circuit
    int i,j; 

    //////////////////////////////////////////////////////////////////////////////
    int grid,block,mingrid;
    float *buf;
    int IIter;
    buf = VFMalloc(A_size);
    printf("<FIELD SOVER>\n");
	printf(" Laplace eq. using PCG\n");
	printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    // Find good grid and block size
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)PCG_LAP,0,Gsize); 
    grid = (Gsize + block - 1) / block;
    printf("minGridSize = %d\n",mingrid);
    printf("blockSize = %d\n",block);
    printf("gridSize = %d\n",grid);
  
    for (i = 0; i < CondNUMR; i++) {
        cudaMemcpy(dev_b, cond_b[i], A_size * sizeof(float),cudaMemcpyHostToDevice);
        for(j=0;j<A_size;j++){
			//if(cond_b[i][j] !=0 )
				//printf("cond_b[%d][%d] = %g\n",i,j,cond_b[i][j]);
		} 
        IIter = 10;
        //printf("IIter = %d\n",IIter);
        //PCG_LAP<<<grid,block>>>(&IIter,Gsize,A_size,dev_A,dev_Ai,dev_Aj,dev_M,dev_AP,dev_R,dev_Z,dev_P,dev_X,dev_b);
        //printf("IIter = %d\n",IIter);
        //cudaMemcpy(buf, dev_R, A_size * sizeof(float),cudaMemcpyDeviceToHost);
        //printf("buf = %d\n",buf[1]);
    }
    //exit(1);
}
void Set_MatrixPCG_cuda(){
    int TID;
    int i,j,k; 
    int PCG_Laplace_SINGLECPU_Flag=0;
    float **CPUsol;
	//vec_cond_Garray = (int *) malloc(ngx * ngy * sizeof(int));
	//vec_boundary_Garray = (int *) malloc(ngx * ngy * sizeof(int));
	//vec_face_Garray = (int *) malloc(ngx * ngy * sizeof(int));
	//vec_area_Garray = (float *) malloc(ngx * ngy * sizeof(float));
	//vec_eps_Carray = (float *) malloc(ncx * ncy * sizeof(float));

    vec_A_idx = (int *) malloc(ngx * ngy * sizeof(int));
    for (i = 0; i < ngx; i++) {
		for (j = 0; j < ngy; j++) {
			vec_A_idx[j + i * ngy] = A_idx[i][j];
			//vec_cond_Garray[j + i * ngy] = cond_Garray[i][j];
			//vec_boundary_Garray[j + i * ngy] = boundary_Garray[i][j];
			//vec_face_Garray[j + i * ngy] = face_Garray[i][j];
			//vec_area_Garray[j + i * ngy] = area_Garray[i][j];
		} // matrix save direction ^ >
	}
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
    // Laplace Solution
    cudaMallocPitch(&Lap_PHI_Sol, &pitch, Gsize * sizeof(float), CondNUMR); // for Laplace Solution
    //cudaMalloc((void**) &Lap_TEMP_Sol, Gsize * sizeof(int));
    // cudaMemset((void *) array, 0, Gsize * sizeof(int));
    
    //Make a Field constant set  
    checkCudaErrors(cudaMalloc((void**)&dev_PCG_const,sizeof(DPS_Const)));
    Make_PCG_Const_Init<<<1,1>>>(dev_PCG_const,A_size,PCGtol);
    Host_PCG_const = (DPS_Const*)malloc(sizeof(DPS_Const));
    checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));
   
    //Make a Field DATA set  
    checkCudaErrors(cudaMalloc((void**)&dev_PCG_DATA, A_size*sizeof(DPS_Data)));
    Make_PCG_DATA_Init<<<A_size/4,4>>>(dev_PCG_DATA,A_size,dev_M);

    Host_PCG_DATA = (DPS_Data*)malloc(A_size*sizeof(DPS_Data));
    checkCudaErrors(cudaMemcpy(Host_PCG_DATA, dev_PCG_DATA, A_size*sizeof(DPS_Data), cudaMemcpyDeviceToHost));
    for (i=0;i<A_size;i++){
        printf("[%d]%g R=%g,Z=%g,P=%g,AP=%g,M=%g\n",i,MatM[i],Host_PCG_DATA[i].vecR,Host_PCG_DATA[i].vecZ,Host_PCG_DATA[i].vecP,Host_PCG_DATA[i].vecAP,Host_PCG_DATA[i].vecM);
    }
    exit(1);
    //printf("testKernel results:\n");
    //printf("point.a: %g, point.b: %g\n",CPU_BUF[0].tol,CPU_BUF[0].tol2);
    //printf("point.a: %g, point.b: %g\n",CPU_BUF[0].tol2,CPU_BUF[0].tol2);
    //printf("point.a: %d, point.b: %g\n",CPU_BUF[0].A_size,CPU_BUF[0].tol);
    //free(CPU_BUF);
    // retrieve the results
    

    
        // deallocate memory
    
    if(PCG_Laplace_SINGLECPU_Flag==1){
		printf(" Preconditioner[Jacovi]\n"); 
        printf(" Main Library set[Single CPU PCG]\n"); 
        printf(" Laplace Equation TEST\n"); 
        X = VFMalloc(A_size);
        B = VFMalloc(A_size);
        R0 = VFMalloc(A_size);
        AX = VFMalloc(A_size);
        Z0 = VFMalloc(A_size);
        P0 = VFMalloc(A_size);
        AP = VFMalloc(A_size);
        PAP = VFMalloc(A_size);
        VFInit(X,0.0,A_size);
        VFInit(B,0.0,A_size);
        VFInit(AX,0.0,A_size);
        VFInit(R0,0.0,A_size);
        VFInit(Z0,0.0,A_size);
        VFInit(P0,0.0,A_size);
        VFInit(AP,0.0,A_size);
        VFInit(PAP,0.0,A_size);
        CPUsol = MFMalloc(CondNUMR,Gsize);
        for (k = 0; k < CondNUMR; k++) {
            VFCopy(B,cond_b[k],A_size);
            FieldIter = PCG_SINGLECPU();
            printf("Solution %d",k);
            printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
            printf("FieldIter = %d\n",FieldIter);
            // Save
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
            /*
            for(j=ngy-1;j>=0;j--){
                for(i=0;i<ngx;i++){
                    TID = i*ngy+j;
                    printf("%6.2g",CPUsol[k][TID]);
                }printf("\n");
            }printf("\n");
            */
            CPU_PCG_Laplace_Solution_Save(CPUsol);
        }
        exit(1);
	}
}
__global__ void Make_PCG_DATA_Init(DPS_Data *p, int size,float *MatrixM){
    int TID = blockIdx.x * blockDim.x + threadIdx.x;
    if(TID>=size) return;
    printf("TID = %d, M = %g\n",TID,MatrixM[TID]);
    p[TID].vecR = 0.0;
    p[TID].vecZ = 0.0;
    p[TID].vecP = 0.0;
    p[TID].vecAP = 0.0;
    p[TID].vecM = MatrixM[TID];
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
int PCG_SINGLECPU(){
    int TID,i,Iter=0;
    float tol2;
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
    tol2 = PCGtol*PCGtol;
    while(rsold>tol2){
        Iter++;
        Temp = 0.0;
        for(TID=0;TID<A_size;TID++){
            AP[TID] = 0;
            for(i=Ai[TID]-Ai[0];i<Ai[TID+1]-Ai[0];i++){
                AP[TID] += A_val[i]*P0[Aj[i]-Ai[0]];
            }
            PAP[TID] = P0[TID] * AP[TID];
            Temp += PAP[TID]; //AtomicAdd!!
        }
        alpha = rsold/Temp;
        for(TID=0;TID<A_size;TID++){
            X[TID] = X[TID] + alpha * P0[TID];
            R0[TID] = R0[TID] - alpha * AP[TID];
            Z0[TID] = MatM[TID] * R0[TID];
            rnew += R0[TID]*Z0[TID];  //AtomicAdd!!
        }
        beta = rnew/rsold;
        for(TID=0;TID<A_size;TID++){ 
            P0[TID] = Z0[TID] + beta*P0[TID];
        }
        rsold = rnew;
        rnew = 0.0;
    }
    return Iter;
}
__global__ void PCG_LAP(int *Iter,int Gsize,int Asize,float *A,int *Ai,int *Aj,float *M,float *AP,float *R,float *Z,float *P,float *X,float *b){
    int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    if(TID>=Asize) return;
    int i;
    //float tol2;
    float a;

    for(i=Ai[TID]-Ai[0];i<Ai[TID+1]-Ai[0];i++){
        AP[TID] += A[i]*X[Aj[i]-Ai[0]];
    }
    R[TID] = b[TID] - AP[TID];
    Z[TID] = M[TID] * R[TID];
    P[TID] = Z[TID];
    a = R[TID]*Z[TID];
    *Iter = 100;
    
    //__syncthreads();
    //printf("TID = %d, a = %g\n",TID,a);  
    //atomicAdd(&rsold,a);
    //if(TID==0) printf("rsold %d= %g\n",TID,rsold);       
    //if(TID==10) printf("rsold %d= %g\n",TID,rsold);  

    /*
    rsold = 0;
    for(TID=0;TID<A_size;TID++){
        AX[TID] = 0;
        for(i=Ai[TID]-Ai[0];i<Ai[TID+1]-Ai[0];i++){
            AX[TID] += A_val[i]*X[Aj[i]-Ai[0]];
        }
        R0[TID] = B[TID] - AX[TID];
        Z0[TID] = MatM[TID] * R0[TID];
        P0[TID] = Z0[TID];
        rsold += R0[TID]*Z0[TID];
    }
    tol2 = PCGtol*PCGtol;
    //printf("rsold=%g\n",rsold);
    while(rsold>tol2){
        Iter++;
        Temp = 0.0;
        for(TID=0;TID<A_size;TID++){
            AP[TID] = 0;
            for(i=Ai[TID]-Ai[0];i<Ai[TID+1]-Ai[0];i++){
                AP[TID] += A_val[i]*P0[Aj[i]-Ai[0]];
                //printf("[%d],[%d] %g, %g\n",TID,Aj[i],A_val[i],P0[TID]);
            }
            //printf("PAP[%d] = %g\n",TID,AP[TID]);
            PAP[TID] = P0[TID] * AP[TID];
            Temp += PAP[TID]; //Reduction
        }
        alpha = rsold/Temp;
        for(TID=0;TID<A_size;TID++){
            X[TID] = X[TID] + alpha * P0[TID];
            R0[TID] = R0[TID] - alpha * AP[TID];
            Z0[TID] = MatM[TID] * R0[TID];
            rnew += R0[TID]*Z0[TID];  //Reduction
        }
        beta = rnew/rsold;
        for(TID=0;TID<A_size;TID++){ 
            P0[TID] = Z0[TID] + beta*P0[TID];
        }
        rsold = rnew;
        rnew = 0.0;
        //printf("alpha %d = %2g, Temp = %g, Beta %d = %2g, rsold=%2g\n",FieldIter-1,alpha,Temp,FieldIter-1,beta,rsold);
    }
    */

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
