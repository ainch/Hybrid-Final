#include "cuda_start.cuh"
void start_cuda(){
	if(External_Flag){
		EFIELD = Efield_cuda;
	}else{
		EFIELD = Efield_cuda_Basic;
	}
	/*** Move ***/
	if(ConstB_Flag){
		//MOVE = MoveB_cuda;
	}else {
		MOVE = Move_cuda;
	}
	SORT_BOUNDARY = SortBounndary_cuda;

	switch(MainGas){
	case ARGON:
		MCC	= MCC_Ar_cuda;
		MCC_Basic = MCC_Ar_cuda;
		break;
	case OXYGEN:
		MCC		= MCC_O2_cuda;
		MCC_Basic = MCC_O2_cuda;
		break;
	case ARO2:
		MCC		= MCC_ArO2_cuda;
		MCC_Basic = MCC_ArO2_cuda;
		break;
	}
	if(CSS_Flag) CONTIEQ = Solve_Continuity_eqn_check;
	else CONTIEQ = Solve_Continuity_eqn;
	DEPOSIT = Deposit_cuda;
	DIAG    = Diagnostic;
}

void info_Device()
{
	int count,i,isp;
	int NeedMemory = 0;
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
	printf("Total constant Mem.: %u Kbytes\n", prop.totalConstMem/1024);
	printf("Shared Mem. per block: %u bytes \n", prop.sharedMemPerBlock);
	printf("Registers available per block.: %u #\n", prop.regsPerBlock);
	printf("Max Mem. pitch: %ld\n", prop.memPitch);
	printf("Multiprocessor count: %d\n", prop.multiProcessorCount);
	printf("Threads in warp: %d \n", prop.warpSize);
	printf("Max thread dimensions: (%d, %d, %d) \n",prop.maxThreadsDim[0], prop.maxThreadsDim[1], prop.maxThreadsDim[2]);
	printf("Max grid dimensions: (%d, %d, %d) \n",prop.maxGridSize[0], prop.maxGridSize[1], prop.maxGridSize[2]);
	printf("---------------------------------------------------\n");	
	checkCudaErrors(cudaGetDeviceProperties(&prop, device_num));
    if (!prop.managedMemory){
        // This sample requires being run on a device that supports Unified Memory
        fprintf(stderr, "Unified Memory not supported on this device\n");
        exit(EXIT_WAIVED);
    }
    // This sample requires being run on a device that supports Cooperative Kernel Launch
    if (!prop.cooperativeLaunch){
        printf("\nSelected GPU (%d) does not support Cooperative Kernel Launch, Waiving the run\n", device_num);
        exit(EXIT_WAIVED);
    }
	// Aprroximation
    for(isp=0;isp<nsp;isp++){
        NeedMemory += 6 * SP[isp].MAXNP*sizeof(float)/1024/1024;
    }
	//eedMemory += 1.1 * Gsize * sizeof(float);
	if(prop.totalGlobalMem/1024/1024<NeedMemory){
		printf("Error : Insufficient GPU memory.\n");
		printf(" GPU Global: %u Mbytes\n", prop.totalGlobalMem/1024/1024);
		printf(" Need Memory: %u Mbytes (Approx.)\n", NeedMemory);
		printf(" You should reduce MaxNP.\n");
		exit(1);
	}
}
