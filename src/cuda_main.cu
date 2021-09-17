#include "cuda_main.cuh"

extern "C" void main_cuda()
{
    int isp;
    cudaEvent_t start, stop;
    printf("-------------GPU_CUDA START-------------\n");
    info_Device();
    start_cuda();
    Set_Device_Parameter();
    Set_Particle_cuda();
    Set_NullCollisionTime_cuda();
	Set_DiagParameter_cuda();
	Set_MatrixPCG_cuda();
	if(Lap_Field_Solver_Test) PCG_Laplace_TEST();
    PCG_SOLVER_Laplace();
	Deposit_cuda();
    while(1){
		PCG_SOLVER();
		Efield_cuda();
        //START
		//(*MOVE)();
		//(*SORT_BOUNDARY)();
		//(*MCC)();
		//if(add_izrate_flag==1) AddIzRate_cuda(); // ADD IONIZATION RATE instead of ICP SOURCE
		//(*DEPOSIT)();
		//if(Meta_flag){// && t>5e-7)
		//	(*CONTIEQ)();
		//}
		//(*DIAG)();
		//SaveDumpFile();
        //
        
        t+=dt; // real time
        tstep++; // step
        if((tstep%CYCLE_NUM) == 0) cstep++; // Number of Cycle step
        printf("TIME = %2.4g (s), STEP = %d (#), CYCLE = %d\n",t,tstep,cstep);
        //if(t>1e-3) break;    
        if(tstep>=1) break;    
        /*
        if((tstep%CYCLE_NUM) == 0){
            cstep++; // Number of Cycle step
            // time calculate
	        while(totaltime > 1000){
		        TotalT_S++;
		        totaltime = totaltime - 1000;
	        }
	        while(TotalT_S >= 60){
			        TotalT_M++;
			        TotalT_S -= 60;
	        }
	        while(TotalT_M >= 60){
			        TotalT_H++;
			        TotalT_M -= 60;
	        }
	        while(TotalT_H >= 24){
			        TotalT_D++;
			        TotalT_H -= 24;
	        }
	        fprintf(stderr, "\nDump at t=%1.5e(s)\n",t);
	        for (isp = 0; isp < nsp; isp++){
		        fprintf(stderr, "%s : %d, ", SP[isp].name, SP[isp].np);
	        }fprintf(stderr, "\n");
	        fprintf(stderr, "Domain size : %d X %d =%d,  ", ngx, ngy, Gsize);
	        fprintf(stderr, "Time: %d(d), %d(h), %d(m), %d(s)\n",TotalT_D,TotalT_H,TotalT_M,TotalT_S);
	        time_sum = gputime_field+gputime_efield+gputime_move+gputime_sort+gputime_mcc+gputime_continue+gputime_deposit+gputime_diag+gputime_trace+gputime_dump;
	        fprintf(stderr, "Total : time = %2.8f	(s)\n", time_sum * 0.001);
	        fprintf(stderr, "Field	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_field * 0.001, gputime_field * 100 / time_sum);
	        fprintf(stderr, "Efield	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_efield * 0.001, gputime_efield * 100 / time_sum);
	        fprintf(stderr, "Move	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_move * 0.001, gputime_move * 100 / time_sum);
	        fprintf(stderr, "Sort	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_sort * 0.001, gputime_sort * 100 / time_sum);
	        fprintf(stderr, "Mcc	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_mcc * 0.001, gputime_mcc * 100 / time_sum);
	        fprintf(stderr, "CONTI	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_continue * 0.001, gputime_continue * 100 / time_sum);
	        fprintf(stderr, "Depo	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_deposit * 0.001, gputime_deposit * 100 / time_sum);
	        fprintf(stderr, "Diag	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_diag * 0.001, gputime_diag * 100 / time_sum);
	        fprintf(stderr, "Dump	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_dump * 0.001, gputime_dump * 100 / time_sum);
	        fprintf(stderr, "------------------------------------------------------------------------------\n");
	
        } 
        */
    }
}


__global__ void testKernel(point *p)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    //printf("B data[%d] %g, %g\n",i,p[i].a,p[i].b);
    p[i].a = 1.1;
    p[i].b = 2.2;
    //printf("A data[%d] %g, %g\n",i,p[i].a,p[i].b);
    
}
__global__ void MakeVectorForMoveKernel(int ngx,int ngy,point *p)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    p[i].a = (float)ngx;
    p[i].b = (float)ngy;
	printf("point.a: %g, point.b: %g\n",p[i].a,p[i].b);
}
int test(void)
{
    // set number of points 
    int numPoints    = 16;
    int gpuBlockSize = 4;
    int pointSize    = sizeof(point);
    int numBytes     = numPoints * pointSize;
    int gpuGridSize  = numPoints / gpuBlockSize;
    point *cpuPointArray;
    point *gpuPointArray;
        // allocate memory
    cpuPointArray = (point*)malloc(numBytes);
    cudaMalloc((void**)&gpuPointArray, numBytes);

    // launch kernel
    testKernel<<<gpuGridSize,gpuBlockSize>>>(gpuPointArray);

    // retrieve the results
    checkCudaErrors(cudaMemcpy(cpuPointArray, gpuPointArray, numBytes, cudaMemcpyDeviceToHost));
    //printf("testKernel results:\n");
    for(int i = 0; i < numPoints; ++i)
    {
        //printf("point.a: %g, point.b: %g\n",cpuPointArray[i].a,cpuPointArray[i].b);
    }
        // deallocate memory
    free(cpuPointArray);
    cudaFree(gpuPointArray);

    return 0;
}
