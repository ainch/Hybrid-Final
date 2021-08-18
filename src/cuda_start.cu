#include "cuda_start.cuh"

void info_Device()
{
	int count,i;
	cudaDeviceProp prop;
	// Determine the number of CUDA capable GPUs
	checkCudaErrors(cudaGetDeviceCount(&count));
	if(count<1){
		printf("No CUDA Capable GPU(s) Detected \n");
		return;
	}
    // Display the GPU processor specification
	printf("<CUDA Parallel Progaraming>\n");
	printf("- Number of CUDA devices:\t%d\n", count);
	for(i=0;i<count;i++){
		cudaDeviceProp dprop;
		checkCudaErrors(cudaGetDeviceProperties(&dprop, i));
		printf(">>>>>>>> Device %d is a %s\n", i, dprop.name);
	}
    printf("GPU selected Device ID = %d \n", device_num);
	checkCudaErrors(cudaGetDeviceProperties (&prop, device_num));
	checkCudaErrors(cudaSetDevice(device_num));
    // This will pick the best possible CUDA capable device
	printf("-------------------- Device%d --------------------\n", device_num);
	printf("Name: %s \n", prop.name);
	printf("Computer capability: %d.%d \n", prop.major, prop.minor);
	//printf("Clock rate: %d\n", prop.clockRate);
	printf("Total Global Mem.: %u Mbytes\n", prop.totalGlobalMem/1024/1024);
	//printf("Total constant Mem.: %u Kbytes\n", prop.totalConstMem/1024);
	printf("Shared Mem. per block: %u bytes \n", prop.sharedMemPerBlock);
	//printf("Registers available per block.: %u #\n", prop.regsPerBlock);
	//printf("Max Mem. pitch: %ld\n", prop.memPitch);
	printf("Multiprocessor count: %d\n", prop.multiProcessorCount);
	printf("Threads in warp: %d \n", prop.warpSize);
	printf("Max thread dimensions: (%d, %d, %d) \n",prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
	printf("Max grid dimensions: (%d, %d, %d) \n",prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
	printf("---------------------------------------------------\n");	
}
void start_cuda(){

	/*** Field solver ***/
	if(PCG_Method==0){
		
	}else if(PCG_Method==1){

	}else if(PCG_Method==2){

	}else if(PCG_Method==3){

	}
	//FieldSolver = PCG_SOLVER;
	/*** Move ***/
	if(ConstB_Flag){
		//MOVE = MoveB_cuda;
	}else {
		//MOVE = Move_cuda;
	}
	switch(MainGas){
	case ARGON:
		//SORT_BOUNDARY = AR_SortAndBoundary_cuda;
		//MCC		= ArMccDiag_cuda;
		//CONTIEQ = Ar_solve_continuity_eqn;
		//DIAG    = Diagnostic;
		break;
	case OXYGEN:
		//SORT_BOUNDARY = Oxy_SortAndBoundary_cuda;
		//MCC		= OxyMccDiag_cuda;
		//CONTIEQ = Oxy_solve_continuity_eqn;
		//DIAG    = Diagnostic;
		break;
	case ARO2:
		//SORT_BOUNDARY = ARO2_SortAndBoundary_cuda;
		//MCC		= ARO2_MccDiag_cuda;
		//CONTIEQ = ARO2_solve_continuity_eqn;
		//DIAG    = Diagnostic;
		break;
	}
	/*** Deposit ***/
	//DEPOSIT = Deposit_cuda;

}
