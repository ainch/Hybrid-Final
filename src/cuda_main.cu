#include "cuda_main.cuh"

extern "C" void main_cuda()
{
    int isp,i,k,sum,dsum,dsum2;
    float dsum3,dsum4;
    cudaEvent_t start, stop;
    float gputime;
    printf("-------------GPU_CUDA START-------------\n");
    info_Device();
    start_cuda();
    Set_Device_Parameter();
    Set_Particle_cuda();
    Set_NullCollisionTime_cuda();
	Set_DiagParameter_cuda();
    Set_SortBoundary_cuda();
	Set_MatrixPCG_cuda();
	if(Lap_Field_Solver_Test) PCG_Laplace_TEST();
    PCG_SOLVER_Laplace();
    Set_Fluid_cuda();
	Deposit_cuda();
    /*
    gputime_field	=0.0;
	gputime_efield	=0.0;
	gputime_diag	=0.0;
	gputime_move	=0.0;
	gputime_mcc		=0.0;
	gputime_continue=0.0;
	gputime_deposit	=0.0;
	gputime_sort	=0.0;
	gputime_trace	=0.0;
	gputime_dump	=0.0;
	totaltime		=0.0;
	TotalT_D		=0;
	TotalT_H		=0;
	TotalT_M		=0;
	TotalT_S		=0;
    */
    int np[5],np2[5],np3[5];
    float np4[5];
    np[0] = 0;np[1] = 0;np[2] = 0;np[3] = 0;np[4] = 0;
    np2[0] = 0;np2[1] = 0;np2[2] = 0;np2[3] = 0;np2[4] = 0;
    np3[0] = 0;np3[1] = 0;np3[2] = 0;np3[3] = 0;np3[4] = 0;
    np4[0] = 0.0f;np4[1] = 0.0f;np4[2] = 0.0f;np4[3] = 0.0f;np4[4] = 0.0f;
    while(1){
        //cudaEventCreate(&start); cudaEventCreate(&stop);
	    //cudaEventRecord( start, 0 );
		PCG_SOLVER();
        //cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    //cudaEventElapsedTime( &gputime, start, stop );
	    //cudaEventDestroy( start );cudaEventDestroy( stop );
        //gputime_field+=gputime;
		//totaltime+=gputime;
        //cudaEventCreate(&start); cudaEventCreate(&stop);
	    //cudaEventRecord( start, 0 );
		Efield_cuda();
        //cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    //cudaEventElapsedTime( &gputime, start, stop );
	    //cudaEventDestroy( start );cudaEventDestroy( stop );
        //gputime_efield+=gputime;
		//totaltime+=gputime;
        //cudaEventCreate(&start); cudaEventCreate(&stop);
	    //cudaEventRecord( start, 0 );
        Move_cuda();
        //cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    //cudaEventElapsedTime( &gputime, start, stop );
	    //cudaEventDestroy( start );cudaEventDestroy( stop );
        //gputime_move+=gputime;
		//totaltime+=gputime;
        //cudaEventCreate(&start); cudaEventCreate(&stop);
	    //cudaEventRecord( start, 0 );
        SortBounndary_cuda();
        //cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    //cudaEventElapsedTime( &gputime, start, stop );
	    //cudaEventDestroy( start );cudaEventDestroy( stop );
        //gputime_sort+=gputime;
		//totaltime+=gputime;
        //cudaEventCreate(&start); cudaEventCreate(&stop);
	    //cudaEventRecord( start, 0 );
        if(MainGas == ARGON) MCC_Ar_cuda();
        if(MainGas == OXYGEN) MCC_O2_cuda();
        if(MainGas == ARO2) MCC_ArO2_cuda();
        //cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    //cudaEventElapsedTime( &gputime, start, stop );
	    //cudaEventDestroy( start );cudaEventDestroy( stop );
        //gputime_mcc+=gputime;
		//totaltime+=gputime;
        //cudaEventCreate(&start); cudaEventCreate(&stop);
	    //cudaEventRecord( start, 0 );
		Deposit_cuda();
        //cudaEventRecord( stop, 0 ); cudaEventSynchronize( stop );
	    //cudaEventElapsedTime( &gputime, start, stop );
	    //cudaEventDestroy( start );cudaEventDestroy( stop );
        //gputime_deposit+=gputime;
		//totaltime+=gputime;
        Tecplot_save();
        //
        t+=dt; // real time
        tstep++; // step
        //if(tstep == 1){
        if((tstep%CYCLE_NUM) == 0){
            cstep++; // Number of Cycle step
            cudaMemcpy(Host_G_sp, dev_G_sp, nsp * Gsize * sizeof(GPG),cudaMemcpyDeviceToHost);
            for(isp=0;isp<nsp;isp++){
                sum = 0;
                dsum = 0;
                dsum2 = 0;
                for(i=0;i<Gsize;i++){
                    if(vec_G[i].DensRegion){
                        sum +=Host_G_sp[isp*Gsize+i].PtNumInCell;
                        dsum +=Host_G_sp[isp*Gsize+i].PtNumMCCInCell;
                        dsum2 +=Host_G_sp[isp*Gsize+i].PtNullMCCInCell;   
                    }
                }
                np[isp] = sum;
                np2[isp] = dsum;
                np3[isp] = dsum2;
            }
            printf("Np = [%d],[%d],[%d],[%d],[%d]\n",np[0],np[1],np[2],np[3],np[4]);
            printf("mcc = [%d],[%d],[%d],[%d],[%d]\n",np2[0],np2[1],np2[2],np2[3],np2[4]);
            printf("null = [%d],[%d],[%d],[%d],[%d]\n",np3[0],np3[1],np3[2],np3[3],np3[4]);
            //printf("[%] = [%3.2g %],[%3.2g %],[%3.2g %],[%3.2g %],[%3.2g %]\n"
            //                        ,100*(float)np2[0]/(np[0]+np2[0]),100*(float)np2[1]/(np[1]+np2[1])
            //                        ,100*(float)np2[2]/(np[2]+np2[2]),100*(float)np2[3]/(np[3]+np2[3])
            //                        ,100*(float)np2[4]/(np[4]+np2[4]));
        }
        printf("TIME = %1.4e (s),[%3d][%3d], Iter = %3d, res = %1.3e\n",t,tstep,cstep,*FIter,*dot_result);
        //if(t>1e-3) break;    
        //if(tstep == 1){
        //if(cstep==200){
            cudaMemcpy(Host_G_sp, dev_G_sp, nsp * Gsize * sizeof(GPG),cudaMemcpyDeviceToHost);
            for(isp=0;isp<nsp;isp++){
                sum = 0;
                dsum = 0;
                dsum2 = 0;
                for(i=0;i<Gsize;i++){
                    if(vec_G[i].DensRegion){
                        sum +=Host_G_sp[isp*Gsize+i].PtNumInCell; 
                        dsum +=Host_G_sp[isp*Gsize+i].PtNumMCCInCell;
                        dsum2 +=Host_G_sp[isp*Gsize+i].PtNullMCCInCell;   
                    }
                }
                np[isp] = sum;
                np2[isp] = dsum;
                np3[isp] = dsum2;
            }
            if(MainGas == ARGON){
                printf("Np = [%d],[%d]\n",np[0],np[1]);
                printf("mcc = [%d],[%d]\n",np2[0],np2[1]);
                printf("null = [%d],[%d]\n",np3[0],np3[1]);
                printf("Coll = [%d],[%d]\n",np2[0]-np3[0],np2[1]-np3[1]);
                printf("[%] = [%3.2g %],[%3.2g %]\n",100*(float)np2[0]/(np[0]+np2[0]),100*(float)np2[1]/(np[1]+np2[1]));
                cudaMemcpy(MCC_rate, dev_MCC_rate, Msize * sizeof(float), cudaMemcpyDeviceToHost);
                for(isp=0;isp<TnRct;isp++){
                    for(i=0;i<Gsize;i++){
                        if(vec_G[i].DensRegion){
                            if(isp<=4){
                                np4[0] += MCC_rate[i*TnRct+isp];
                            }else if(isp>4 && isp<=6){
                                np4[1] += MCC_rate[i*TnRct+isp];
                            }
                        }
                    }
                    printf("R[%d] = [%g],[%g]\n",isp,np4[0],np4[1]);
                    np4[0] = 0.0f;np4[1] = 0.0f;
                }
            }
            if(MainGas == OXYGEN){
                printf("Np = [%d],[%d],[%d],[%d]\n",np[0],np[1],np[2],np[3]);
                printf("mcc = [%d],[%d],[%d],[%d]\n",np2[0],np2[1],np2[2],np2[3]);
                printf("null = [%d],[%d],[%d],[%d]\n",np3[0],np3[1],np3[2],np3[3]);
                printf("Coll = [%d],[%d],[%d],[%d]\n",np2[0]-np3[0],np2[1]-np3[1],np2[2]-np3[2],np2[3]-np3[3]);
                printf("[%] = [%3.2g %],[%3.2g %],[%3.2g %],[%3.2g %]\n"
                                    ,100*(float)np2[0]/(np[0]+np2[0]),100*(float)np2[1]/(np[1]+np2[1])
                                    ,100*(float)np2[2]/(np[2]+np2[2]),100*(float)np2[3]/(np[3]+np2[3]));
                cudaMemcpy(MCC_rate, dev_MCC_rate, Msize * sizeof(float), cudaMemcpyDeviceToHost);
                for(isp=0;isp<TnRct;isp++){
                    for(i=0;i<Gsize;i++){
                        if(vec_G[i].DensRegion){
                            if(isp<=40){
                                np4[0] += MCC_rate[i*TnRct+isp];
                            }else if(isp>40 && isp<=46){
                                np4[3] += MCC_rate[i*TnRct+isp];
                            }else if(isp>46 && isp<=52){
                                np4[1] += MCC_rate[i*TnRct+isp];
                            }else if(isp>52 && isp<=57){
                                np4[2] += MCC_rate[i*TnRct+isp];
                            }
                        }
                    }
                    printf("R[%d] = [%g],[%g],[%g],[%g]\n",isp,np4[0],np4[1],np4[2],np4[3]);
                    np4[0] = 0.0f;np4[1] = 0.0f;np4[2] = 0.0f;np4[3] = 0.0f;
                }
            }
            if(MainGas == ARO2){
                printf("Np = [%d],[%d],[%d],[%d],[%d]\n",np[0],np[1],np[2],np[3],np[4]);
                //printf("mcc = [%d],[%d],[%d],[%d],[%d]\n",np2[0],np2[1],np2[2],np2[3],np2[4]);
                //printf("null = [%d],[%d],[%d],[%d],[%d]\n",np3[0],np3[1],np3[2],np3[3],np3[4]);
                printf("   Coll = [%d],[%d],[%d],[%d],[%d]\n",np2[0]-np3[0],np2[1]-np3[1],np2[2]-np3[2],np2[3]-np3[3],np2[4]-np3[4]);
                printf("[%] = [%3.2g %],[%3.2g %],[%3.2g %],[%3.2g %],[%3.2g %]\n"
                                    ,100*(float)np2[0]/(np[0]+np2[0]),100*(float)np2[1]/(np[1]+np2[1])
                                    ,100*(float)np2[2]/(np[2]+np2[2]),100*(float)np2[3]/(np[3]+np2[3])
                                    ,100*(float)np2[4]/(np[4]+np2[4]));
                cudaMemcpy(MCC_rate, dev_MCC_rate, Msize * sizeof(float), cudaMemcpyDeviceToHost);
                for(isp=0;isp<TnRct;isp++){
                    for(i=0;i<Gsize;i++){
                        if(vec_G[i].DensRegion){
                            if(isp<=45){
                                np4[0] += MCC_rate[i*TnRct+isp];
                            }else if(isp>45 && isp<=51){
                                np4[4] += MCC_rate[i*TnRct+isp];
                            }else if(isp>51 && isp<=59){
                                np4[2] += MCC_rate[i*TnRct+isp];
                            }else if(isp>59 && isp<=64){
                                np4[3] += MCC_rate[i*TnRct+isp];
                            }else if(isp>64 && isp<=67){
                                np4[1] += MCC_rate[i*TnRct+isp];
                            }
                        }
                    }
                    //printf("R[%d] = [%g],[%g],[%g],[%g],[%g]\n",isp,np4[0],np4[1],np4[2],np4[3],np4[4]);
                    np4[0] = 0.0f;np4[1] = 0.0f;np4[2] = 0.0f;np4[3] = 0.0f;np4[4] = 0.0f;
                }
            }
            //time_sum = gputime_field+gputime_efield+gputime_move+gputime_sort+gputime_mcc+gputime_continue+gputime_deposit+gputime_diag+gputime_trace+gputime_dump;
            //fprintf(stderr, "\n");
	        //fprintf(stderr, "Total : time = %2.8f	(s)\n", time_sum * 0.001);
	        //fprintf(stderr, "Field	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_field * 0.001, gputime_field * 100 / time_sum);
	        //fprintf(stderr, "Efield	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_efield * 0.001, gputime_efield * 100 / time_sum);
	        //fprintf(stderr, "Move	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_move * 0.001, gputime_move * 100 / time_sum);
	        //fprintf(stderr, "Sort	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_sort * 0.001, gputime_sort * 100 / time_sum);
	        //fprintf(stderr, "Mcc	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_mcc * 0.001, gputime_mcc * 100 / time_sum);
	        //fprintf(stderr, "Depo	: time = %2.8f	(s)		rate = %g	(%)\n",	gputime_deposit * 0.001, gputime_deposit * 100 / time_sum);
	        //fprintf(stderr, "------------------------------------------------------------------------------\n");
            //break; 
            //exit(1);
        //}    
        if(isnan(*dot_result) || isinf(*dot_result)){
            printf("\n");
            cudaMemcpy(Host_G_sp, dev_G_sp, nsp * Gsize * sizeof(GPG),cudaMemcpyDeviceToHost);
            for(isp=0;isp<nsp;isp++){
                sum = 0;
                dsum3 = 0.0;
                dsum4 = 0.0;
                k = 0;
                for(i=0;i<Gsize;i++){
                    if(vec_G[i].DensRegion){
                        k++;
                        sum += Host_G_sp[isp*Gsize+i].PtNumInCell;
                        dsum3 += Host_G_sp[isp*Gsize+i].den * SP[isp].np2c/dx/dy;
                        dsum4 += Host_G_sp[isp*Gsize+i].sigma;
                        if(isnan(Host_G_sp[isp*Gsize+i].sigma)) printf("\tsigma[%d].[%d] = [%g]\n",isp,i,Host_G_sp[isp*Gsize+i].sigma);
                    }
                }
                printf("\tNP - %s : %d, %g, %g\n",SP[isp].name,sum,dsum3/k,dsum4/k);
                if(dsum3/k == 0) exit(1);
                //if(isp !=0 && isnan(dsum)) exit(1);
            }
            exit(1);
        }
        
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
