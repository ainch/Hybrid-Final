#include "cuda_Field.cuh"

void PCG_SOLVER_Laplace(){
    // Solve Laplace Equation. (To use every time step.)
    // Goal
    // Lap_TEMP_Sol[Gsize] : Temperature Profile
    // Lap_PHI_Sol[CondNUMR][Gsize] : Each of conductor Phi Profile, This is Device value
    // Lap_SIG_Sol[CondNUMR][CondNUMR] : Each of conductor Sigma Profile for external circuit
    int i,j,k,TID; 

    //////////////////////////////////////////////////////////////////////////////
    int grid,block,mingrid;
    int IIter;
    float *buf;
    float **CPUsol;
    CPUsol = MFMalloc(CondNUMR,Gsize);
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
  
    for (k = 0; k < CondNUMR-1; k++) {
        checkCudaErrors(cudaMemcpy(dev_b, cond_b[k], A_size * sizeof(float),cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemset((void *) dev_X, 0, A_size * sizeof(float)));
        //  Method 1
        //PCG_LAP<<<grid,block>>>(dev_A,dev_Ai,dev_Aj,dev_PCG_const,dev_PCG_DATA,dev_X,dev_b);
        //  Method 2
        FieldIter = PCG_LAP_Divide(grid,block);
        //
        checkCudaErrors(cudaMemcpy(buf, dev_X, A_size * sizeof(float),cudaMemcpyDeviceToHost));
        //checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));
        //checkCudaErrors(cudaMemcpy(Host_PCG_DATA, dev_PCG_DATA, A_size*sizeof(DPS_Data), cudaMemcpyDeviceToHost));
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
        /*
        for(j=ngy-1;j>=0;j--){
            for(i=0;i<ngx;i++){
                TID = i*ngy+j;
                printf("%6.2g",CPUsol[k][TID]);
            }printf("\n");
        }printf("\n");
        */
    }
    exit(1);
}
void Set_MatrixPCG_cuda(){
    int TID;
    int i,j,k; 
    int PCG_Laplace_SINGLECPU_Flag=1;
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
    Host_PCG_const = (DPS_Const*)malloc(sizeof(DPS_Const));
    checkCudaErrors(cudaMalloc((void**)&dev_PCG_const,sizeof(DPS_Const)));
    Make_PCG_Const_Init<<<1,1>>>(dev_PCG_const,A_size,PCGtol);
    checkCudaErrors(cudaMemcpy(Host_PCG_const, dev_PCG_const, sizeof(DPS_Const), cudaMemcpyDeviceToHost));

    //Make a Field DATA set  
    Host_PCG_DATA = (DPS_Data*)malloc(A_size*sizeof(DPS_Data));
    checkCudaErrors(cudaMalloc((void**)&dev_PCG_DATA, A_size*sizeof(DPS_Data)));
    Make_PCG_DATA_Init<<<A_size/4,4>>>(dev_PCG_DATA,A_size,dev_M);
    
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
            VFInit(X,0.0,A_size);
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
        //exit(1);
	}
}
int PCG_LAP_Divide(int grid,int block){
    int Iter = 0;
    float *dev_rsold;
    float *dev_Temp;
    float *dev_rnew;
    float rsold,rnew,Temp;
    float alpha,beta;
    
    checkCudaErrors(cudaMalloc((void**) &dev_rsold,sizeof(dev_rsold)));
    checkCudaErrors(cudaMemset((void *) dev_rsold, 0,  sizeof(dev_rsold)));
    checkCudaErrors(cudaMalloc((void**) &dev_Temp,sizeof(dev_Temp)));
    checkCudaErrors(cudaMemset((void *) dev_Temp, 0,  sizeof(dev_Temp)));
    checkCudaErrors(cudaMalloc((void**) &dev_rnew,sizeof(dev_rnew)));
    checkCudaErrors(cudaMemset((void *) dev_rnew, 0,  sizeof(dev_rnew)));

    Make_PCG_DATA_Init<<<A_size/4,4>>>(dev_PCG_DATA,A_size,dev_M);
    PCG_LAP_Part0<<<grid,block>>>(dev_A,dev_Ai,dev_Aj,dev_PCG_const,dev_PCG_DATA,dev_X,dev_b,dev_rsold);
    cudaMemcpy(&rsold, dev_rsold, sizeof(dev_rsold), cudaMemcpyDeviceToHost);
    printf("Iter = %d, rsold = %g\n",Iter,rsold);
    while(rsold>Host_PCG_const[0].tol2){
        Iter++;
        PCG_LAP_Part1<<<grid,block>>>(dev_A,dev_Ai,dev_Aj,dev_PCG_const,dev_PCG_DATA,dev_Temp);
        cudaMemcpy(&Temp, dev_Temp, sizeof(dev_Temp), cudaMemcpyDeviceToHost);
        alpha = (Temp) ? rsold/Temp : 0.0f;
        printf("Iter = %d, alpha = %g\n",Iter,alpha);
        PCG_LAP_Part2<<<grid,block>>>(dev_PCG_const,dev_PCG_DATA,dev_X,alpha, dev_rnew);
        cudaMemcpy(&rnew, dev_rnew, sizeof(dev_rnew), cudaMemcpyDeviceToHost);
        beta = (rsold) ? rnew/rsold : 0.0f;
        printf("Iter = %d, beta = %g\n",Iter,beta);
        PCG_LAP_Part3<<<grid,block>>>(dev_PCG_const,dev_PCG_DATA,beta);
        rsold = rnew;
        cudaMemset((void *) dev_Temp, 0,  sizeof(dev_Temp));
        cudaMemset((void *) dev_rnew, 0,  sizeof(dev_rnew));
        if(Iter<10) printf("Iter = %d, Temp = %g, alpha = %g, beta = %g, rsold = %g\n",Iter,Temp,alpha,beta,rsold);
    }
    return Iter;
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
int PCG_SINGLECPU(){
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
            Z0[TID] = MatM[TID] * R0[TID];
            rnew += R0[TID]*Z0[TID];  //AtomicAdd!!
        }
        beta = (rsold) ? rnew/rsold: 0.0f;
        for(TID=0;TID<A_size;TID++){ 
            P0[TID] = Z0[TID] + beta*P0[TID];
        }
        rsold = rnew;
        rnew = 0.0;
        if(Iter<10) printf("Iter = %d, Temp = %g, alpha = %g, beta = %g, rsold = %g\n",Iter,Temp,alpha,beta,rsold);
    }
    return Iter;
}
__global__ void PCG_LAP(float *A,int *Ai,int *Aj,DPS_Const *PCG_C,DPS_Data *PCG_D,float *X,float *b){
    int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    if(TID>=PCG_C[0].A_size) return;
    int i;
    float sum;
    float *rs1,*rs2,*rs3;
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
__global__ void PCG_LAP_Part0(float *A,int *Ai,int *Aj,DPS_Const *PCG_C,DPS_Data *PCG_D,float *X,float *b,float *rsold){
    int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    if(TID>=PCG_C[0].A_size) return;
    int i;
    float sum;
    // cal  AP = A * P
    for(i=Ai[TID]-1;i<Ai[TID+1]-1;i++){
        PCG_D[TID].AP += A[i] * X[Aj[i]-1];
    }
    PCG_D[TID].R = b[TID] - PCG_D[TID].AP;
    PCG_D[TID].Z = PCG_D[TID].M * PCG_D[TID].R;
    PCG_D[TID].P = PCG_D[TID].Z;
    sum = PCG_D[TID].R*PCG_D[TID].Z;
    atomicAdd(rsold,sum); 
}
__global__ void PCG_LAP_Part1(float *A,int *Ai,int *Aj,DPS_Const *PCG_C,DPS_Data *PCG_D,float *Temp){
    int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    if(TID>=PCG_C[0].A_size) return;
    int i;
    float sum;
    PCG_D[TID].AP = 0;
    for(i=Ai[TID]-1;i<Ai[TID+1]-1;i++){
        PCG_D[TID].AP += A[i] * PCG_D[Aj[i]-1].P;
    }
    sum = PCG_D[TID].P * PCG_D[TID].AP;
    atomicAdd(Temp,sum);
}
__global__ void PCG_LAP_Part2(DPS_Const *PCG_C,DPS_Data *PCG_D,float *X,float alpha,float *rnew){
    int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    if(TID>=PCG_C[0].A_size) return;
    int i;
    float sum;
    printf("alpha = %g\n",alpha);
    X[TID] = X[TID] + alpha * PCG_D[TID].P;
    PCG_D[TID].R = PCG_D[TID].R - alpha * PCG_D[TID].AP;
    PCG_D[TID].Z = PCG_D[TID].M * PCG_D[TID].R;
    sum = PCG_D[TID].R * PCG_D[TID].Z;
    atomicAdd(rnew,sum);
}
__global__ void PCG_LAP_Part3(DPS_Const *PCG_C,DPS_Data *PCG_D,float beta){
    int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    if(TID>=PCG_C[0].A_size) return;
    PCG_D[TID].P = PCG_D[TID].Z + beta * PCG_D[TID].P;
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
