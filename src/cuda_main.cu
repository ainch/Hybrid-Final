#include "cuda_main.cuh"
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"\nGPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}
extern "C" void main_cuda()
{
    int isp,i,k,sum;
    float dsum,dsum2;
    int KEY0, KEY1, KEY2;
    printf("-------------GPU_CUDA START-------------\n");
    info_Device();
    start_cuda();
    Set_Device_Parameter();
    Set_Particle_cuda();
    Set_NullCollisionTime_cuda();
	Set_Diagnostic_cuda();
    Set_SortBoundary_cuda();
	Set_MatrixPCG_cuda();
	if(Lap_Field_Solver_Test) PCG_Laplace_TEST();
    PCG_SOLVER_Laplace();
    Set_Fluid_cuda();
	Deposit_cuda();
    while(cstep<Basic_Flag){
        t+=dt; // real time
        tstep++; // step
        if((tstep%CYCLE_NUM) == 0) cstep++;
        PCG_SOLVER();
        Efield_cuda();
        Move_cuda();
        SortBounndary_cuda();
        if(MainGas == ARGON) MCC_Ar_cuda();
        else if(MainGas == OXYGEN) MCC_O2_cuda();
        else if(MainGas == ARO2) MCC_ArO2_cuda();
        Deposit_cuda();
        Diagnostic();
        SaveDumpFile(0,0,0);
        printf("TIME = %1.4e (s), Iter = %3d, res = %1.3e\r",t,*FIter,*dot_result);
        if(isnan(*dot_result) || isinf(*dot_result)){
            printf("\nField solver Error!\n");
            exit(1);
        }
    }
    KEY2 = 0, KEY1 = 0, KEY0 = 0; // Save DumpFile version setting; 
    while(t<1e-3){
        t+=dt; // real time
        tstep++; // step
        if((tstep%CYCLE_NUM) == 0) cstep++;
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
	    cudaEventRecord( start, 0 );
		PCG_SOLVER();
        cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    cudaEventElapsedTime( &gputime, start, stop );
	    cudaEventDestroy( start );cudaEventDestroy( stop );
        gputime_field+=gputime;
		totaltime+=gputime;        
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        fprintf(stderr,"TIME = %1.4e (s),[%3d][%3d], Iter = %3d, res = %1.3e\r",t,tstep,cstep,*FIter,*dot_result);
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
	    cudaEventRecord( start, 0 );
		Efield_cuda();
        cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    cudaEventElapsedTime( &gputime, start, stop );
	    cudaEventDestroy( start );cudaEventDestroy( stop );
        gputime_efield+=gputime;
		totaltime+=gputime;        
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
	    cudaEventRecord( start, 0 );
        (*MOVE)();
        cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    cudaEventElapsedTime( &gputime, start, stop );
	    cudaEventDestroy( start );cudaEventDestroy( stop );
        gputime_move+=gputime;
		totaltime+=gputime;        
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
	    cudaEventRecord( start, 0 );
        (*SORT_BOUNDARY)();
        cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    cudaEventElapsedTime( &gputime, start, stop );
	    cudaEventDestroy( start );cudaEventDestroy( stop );
        gputime_sort+=gputime;
		totaltime+=gputime;        
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
	    cudaEventRecord( start, 0 );
        (*MCC)();
        cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    cudaEventElapsedTime( &gputime, start, stop );
	    cudaEventDestroy( start );cudaEventDestroy( stop );
        gputime_mcc+=gputime;
		totaltime+=gputime;        
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
	    cudaEventRecord( start, 0 );
		(*DEPOSIT)();
        cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    cudaEventElapsedTime( &gputime, start, stop );
	    cudaEventDestroy( start );cudaEventDestroy( stop );
        gputime_deposit+=gputime;
		totaltime+=gputime;        
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
		cudaEventRecord( start, 0 );
		//if(Conti_Flag) (*CONTIEQ)();
		cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
		cudaEventElapsedTime( &gputime, start, stop );
		cudaEventDestroy( start );cudaEventDestroy( stop );
		gputime_continue+=gputime;
		totaltime+=gputime;
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
	    cudaEventRecord( start, 0 );
        (*DIAG)();
        cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    cudaEventElapsedTime( &gputime, start, stop );
	    cudaEventDestroy( start );cudaEventDestroy( stop );
        gputime_diag+=gputime;
		totaltime+=gputime;        
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
	    cudaEventRecord( start, 0 );
        Tecplot_save();
        cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    cudaEventElapsedTime( &gputime, start, stop );
	    cudaEventDestroy( start );cudaEventDestroy( stop );
        gputime_Tec+=gputime;
		totaltime+=gputime;        
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        ///////////////////////////////////////////////////////////////////////////
        cudaEventCreate(&start); cudaEventCreate(&stop);
		cudaEventRecord( start, 0 );
		SaveDumpFile(KEY2,KEY1,KEY0);
		cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
		cudaEventElapsedTime( &gputime, start, stop );
		cudaEventDestroy( start );cudaEventDestroy( stop );
		gputime_dump+=gputime;
		totaltime+=gputime;        
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );
        ///////////////////////////////////////////////////////////////////////////  
        if(isnan(*dot_result) || isinf(*dot_result)){
            printf("\nField solver Error!\n");
            exit(1);
        }
        //if(tstep > 3) exit(1);
    }
}