#include "cuda_Field.cuh"


void PCG_SOLVER_Laplace(){
    // Solve Laplace Equation. (To use every time step.)
    // Goal
    // Lap_TEMP_Sol[Gsize] : Temperature Profile
    // Lap_PHI_Sol[CondNUMR][Gsize] : Each of conductor Phi Profile, This is Device value
    // Lap_SIG_Sol[CondNUMR][CondNUMR] : Each of conductor Sigma Profile for external circuit

    float *dev_A, *dev_b, *dev_R, *dev_P;				// PCG device parameter
    float *dev_AP, *dev_M, *dev_Z, *dev_X, *dev_Tmp;	// PCG device parameter
    float *dev_phi_dw;
    float *dev_phi_u;

    cudaMalloc((void**) &dev_b, A_size * sizeof(float));
    cudaMalloc((void**) &dev_phi_dw, ngx * sizeof(float));
	cudaMalloc((void**) &dev_phi_u, ngx * sizeof(float));

	int i;
    fprintf(stderr, "<FIELD SOVER>\n");
	fprintf(stderr, " Laplace eq. using PCG\n");
	fprintf(stderr, " Matrix Size = %d X %d = %d\n", A_size, A_size, A_size*A_size);

    PCGtol *= 1e-3;
    for (i = 0; i < CondNUMR; i++) {
        cudaMemcpy(dev_b, cond_b[i], A_size * sizeof(float),cudaMemcpyHostToDevice);
		//cudaMemcpy(dev_phi_dw, phi_dw[i], ngx * sizeof(float),cudaMemcpyHostToDevice);
		//cudaMemcpy(dev_phi_u, phi_u[i], ngx * sizeof(float),cudaMemcpyHostToDevice);
    }
    PCGtol *= 1e3;
}
void Set_MatrixPCG_cuda(){
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
    float *dev_phi_dw;
    float *dev_phi_u;

    // Laplace Solution
    cudaMallocPitch(&Lap_PHI_Sol, &pitch, Gsize * sizeof(float), CondNUMR); // for Laplace Solution
    //cudaMalloc((void**) &Lap_TEMP_Sol, Gsize * sizeof(int));
   // cudaMemset((void *) array, 0, Gsize * sizeof(int));
}
__global__ void SaveAT2D(float *A, size_t pitch, int height, float *PHI, int n)
{
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;

	if(TID>=n) return;

	float *row=(float *)((char *)A+height*pitch);

	row[TID]=PHI[TID];
}
__global__ void LoadAT2D(float *A, size_t pitch, int height, float *PHI, int n)
{
	int TID=blockDim.x*(gridDim.x*blockIdx.y+blockIdx.x)+threadIdx.x;

	if(TID>=n) return;

	float *row=(float *)((char *)A+height*pitch);

	PHI[TID]=row[TID];
}