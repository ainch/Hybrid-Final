#include "cuda_Field.cuh"

void PCG_SOLVER_Laplace(){
    // Solve Laplace Equation. (To use every time step.)
    // Goal
    // Lap_TEMP_Sol[Gsize] : Temperature Profile
    // Lap_PHI_Sol[CondNUMR][Gsize] : Each of conductor Phi Profile, This is Device value
    // Lap_SIG_Sol[CondNUMR][CondNUMR] : Each of conductor Sigma Profile for external circuit
    int i,j;
    float *dev_A, *dev_b, *dev_R, *dev_P;	// PCG device parameter
    int *dev_Aj,*dev_Ai;			
    float *dev_AP, *dev_M, *dev_Z, *dev_X, *dev_Tmp;	// PCG device parameter
    int   *vec_A_idx;
    int   *vec_cond_Garray;
    int   *vec_boundary_Garray;
    int   *vec_face_Garray;
    float *vec_area_Garray;
    float *vec_eps_Carray;
    float *dev_Sigma;
    int   *dev_face_Garray;
    float *dev_area_Garray;
    float *dev_eps_Carray;
    cudaMalloc((void**) &dev_A, 5 * A_size * sizeof(float));
	cudaMalloc((void**) &dev_Aj, 5 * A_size * sizeof(int));
	cudaMalloc((void**) &dev_Ai, (A_size + 1) * sizeof(int));
	cudaMalloc((void**) &dev_b, A_size * sizeof(float));
	cudaMalloc((void**) &dev_R, A_size * sizeof(float));
	cudaMalloc((void**) &dev_Tmp, A_size * sizeof(float));
	cudaMalloc((void**) &dev_Z, A_size * sizeof(float));
	cudaMalloc((void**) &dev_P, A_size * sizeof(float));
	cudaMalloc((void**) &dev_AP, A_size * sizeof(float));
	cudaMalloc((void**) &dev_X, A_size * sizeof(float));
	cudaMalloc((void**) &dev_M, A_size * sizeof(float));
    cudaMalloc((void**) &dev_b, A_size * sizeof(float));
    // Initialize
    cudaMemset((void *) dev_X, 0, A_size * sizeof(float));
	cudaMemset((void *) dev_AP, 0, A_size * sizeof(float));
    //Copy
    cudaMemcpy(dev_A, A_val, 5 * A_size * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Aj, Aj, 5 * A_size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_Ai, Ai, (A_size + 1) * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_M, MatM, A_size * sizeof(float), cudaMemcpyHostToDevice);
	vec_A_idx = (int *) malloc(ngx * ngy * sizeof(int));
	//vec_cond_Garray = (int *) malloc(ngx * ngy * sizeof(int));
	//vec_boundary_Garray = (int *) malloc(ngx * ngy * sizeof(int));
	//vec_face_Garray = (int *) malloc(ngx * ngy * sizeof(int));
	//vec_area_Garray = (float *) malloc(ngx * ngy * sizeof(float));
	//vec_eps_Carray = (float *) malloc(ncx * ncy * sizeof(float));
	for (i = 0; i < ngx; i++) {
		for (j = 0; j < ngy; j++) {
			vec_A_idx[j + i * ngy] = A_idx[i][j];
			//vec_cond_Garray[j + i * ngy] = cond_Garray[i][j];
			//vec_boundary_Garray[j + i * ngy] = boundary_Garray[i][j];
			//vec_face_Garray[j + i * ngy] = face_Garray[i][j];
			//vec_area_Garray[j + i * ngy] = area_Garray[i][j];
		} // matrix save direction ^ >
	}
    //////////////////////////////////////////////////////////////////////////////
	
    int grid,block,mingrid;
    // Find good grid and block size
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)PCG,0,Gsize); 
    grid = (Gsize + block - 1) / block;
    printf("minGridSize = %d\n",mingrid);
    printf("blockSize = %d\n",block);
    printf("gridSize = %d\n",grid);

    PCGtol *= 1e-3;
    for (i = 0; i < CondNUMR; i++) {
        cudaMemcpy(dev_b, cond_b[i], A_size * sizeof(float),cudaMemcpyHostToDevice);
        for(j=0;j<A_size;j++){
			//if(cond_b[i][j] !=0 )
				//printf("cond_b[%d][%d] = %g\n",i,j,cond_b[i][j]);
		} 
        PCG<<<grid,block>>>(FieldIter,Gsize,A_size,dev_A,dev_Ai,dev_Aj,dev_X,dev_b);
        printf("FieldIter = %d\n",FieldIter);
    }
    PCGtol *= 1e3;
    exit(1);
}
void Set_MatrixPCG_cuda(){
    

    float *dev_phi_dw;
    float *dev_phi_u;
    printf("<FIELD SOVER>\n");
	printf(" Laplace eq. using PCG\n");
	printf(" Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);
    if(PCG_Method==0){
		printf(" Preconditioner[Jacovi]\n"); 
        printf(" Main Library set[cusparseScsrmv]\n"); 
	}else if(PCG_Method==1){
        printf(" Preconditioner[Jacovi]\n");
        printf(" Main Library set[cusparseSpMV]\n"); 
	}else if(PCG_Method==2){
        printf(" Preconditioner[Jacovi]\n");
        printf(" Main Library set[cusparseSpMV+GraphLaunch]\n"); 
	}else if(PCG_Method==3){

	}

    // Laplace Solution
    cudaMallocPitch(&Lap_PHI_Sol, &pitch, Gsize * sizeof(float), CondNUMR); // for Laplace Solution
    //cudaMalloc((void**) &Lap_TEMP_Sol, Gsize * sizeof(int));
   // cudaMemset((void *) array, 0, Gsize * sizeof(int));
}
__global__ void PCG(int Iter,int Gsize,int Asize,float *A,int *Ai,int *Aj,float *X,float *b){
    int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;
    if(TID>=Asize) return;
    int i;
    float AP;
    // r0 = b-AX
    // Data access
    // TID = 0~(Asize-1)
    // St_ID = Ai[0];
    // ID = Ai[TID] - 1;
    // Dn = Ai[TID]-Ai[TID-1]

    
    Iter = TID;
    //printf("Iter = %d\n",Iter);
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