#include "cuda_Field.cuh"
//
void PCG_SOLVER_Laplace(){
    int i,j,k,TID;
    // OUTPUT
    // Lap_TEMP_Sol[Gsize] : Temperature Profile
    // Lap_PHI_Sol[CondNUMR][Gsize] : Each of conductor Phi Profile, This is Device value
    // Lap_SIG_Sol[CondNUMR][CondNUMR] : Each of conductor Sigma Profile for external circuit
    void *kernelArgs[] = {
        (void*)&dev_Ai,
        (void*)&dev_Aj,
        (void*)&dev_A,
        (void*)&dev_X,
        (void*)&dev_M,
        (void*)&dev_AP,
        (void*)&dev_P,
        (void*)&dev_R,
        (void*)&dev_Z,
        (void*)&N,
        (void*)&nz,
        (void*)&PCGtol2,
        (void*)&FIter,
        (void*)&dot_result,
    };
    for (k = 0; k < CondNUMR; k++) {
        checkCudaErrors(cudaMemcpy(dev_R, cond_b[k], N * sizeof(float),cudaMemcpyHostToDevice));
        checkCudaErrors(cudaMemset((void *) dev_X, 0, N * sizeof(float)));
        checkCudaErrors(cudaLaunchCooperativeKernel((void *)PCG,FIELD_GRID,FIELD_BLOCK, kernelArgs, sMemSize, NULL));
        checkCudaErrors(cudaDeviceSynchronize());
        printf(" Laplace Solution %d",k);
        printf(" : Conductor %d = 1 V, Other CondUCTOR = 0 V\n",k);
        printf(" - Iter = %d (ms), rsold^2 = %g\n",*FIter,*dot_result);

    }

    // Laplace Solution
    //cudaMallocPitch(&Lap_PHI_Sol, &pitch, Gsize * sizeof(float), CondNUMR); // for Laplace Solution
    //cudaMalloc((void**) &Lap_TEMP_Sol, Gsize * sizeof(int));
    // cudaMemset((void *) array, 0, Gsize * sizeof(int));

    exit(1);
}
void Set_MatrixPCG_cuda(){
    int i,j;
    N = A_size;
    nz = 5 * A_size;
    printf(" Field Solver : [GPU] Preconditioned Conjugate Gradient\n"); 
    printf(" Laplace Equation\n"); 
    printf(" Matrix Size = %d, ngx x ngy = %d X %d = %d\n", N, ngx, ngy, Gsize);
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
    //Unified memory value for Field residual
    cudaMallocManaged((void **)&dot_result, sizeof(double));
    *dot_result = 0.0;
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
    checkCudaErrors(cudaMemcpy(dev_A_idx, A_idx, Gsize * sizeof(int),cudaMemcpyHostToDevice));
}
__device__ void Mat_x_Vec(int *I, int *J, float *val, int nnz, int num_rows, float alpha, float *inputVecX, 
                        float *outputVecY, cg::thread_block &cta, const cg::grid_group &grid){
    for (int i=grid.thread_rank(); i < num_rows; i+= grid.size())    {
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
    for (int i=grid.thread_rank(); i < size; i+=grid.size()) vecC[i] = (vecA[i] * vecB[i]);
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
    //if(r[TID] !=0) printf("r[%d] = %g\n",TID,r[TID]);
    cg::sync(grid);
    Vec_Dot_Sum(r, Z, d_result, N, cta, grid); 
    cg::sync(grid);
    rsold = *d_result;
    //if(threadIdx.x == 0 && blockIdx.x == 0) printf("First:rsold = %g N = %d, nnz = %d\n",rsold,N,nnz);
    //return;
    while (rsold > tol2 && *Iter <= max_iter){
        Mat_x_Vec(I, J, val, nnz, N, a, p, Ax, cta, grid);
        if (threadIdx.x == 0 && blockIdx.x == 0){
            *Iter = *Iter + 1;
            *d_result = 0.0;  
        } 
        cg::sync(grid);
        //if(Ax[TID] !=0) printf("Ax[%d] = %g\n",TID,Ax[TID]);
        Vec_Dot_Sum(p, Ax, d_result, N, cta, grid);
        cg::sync(grid);
        Temp = *d_result;
        //if(threadIdx.x == 0 && blockIdx.x == 0) printf("Temp = %g\n",Temp);
        //return;
        alpha = (Temp)? rsold/Temp:0.0f;
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
        rnew = 0.0;
        //if(threadIdx.x == 0 && blockIdx.x == 0 && k<20) printf("Iter = %d, temp = %g,  AL = %g, BE = %g Res = %g\n",k,Temp,alpha,beta,rsold);
    }
    //if(threadIdx.x == 0 && blockIdx.x == 0 ) printf("End Iter = %d, Res = %g, b = %g, a = %g\n",k,Temp,alpha,beta,rsold);
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