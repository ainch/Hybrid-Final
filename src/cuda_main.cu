#include "cuda_main.cuh"

extern "C" void main_cuda()
{
    info_Device();
    start_cuda();
    test();
}


__global__ void testKernel(point *p)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    printf("B data[%d] %g, %g\n",i,p[i].a,p[i].b);
    p[i].a = 1.1;
    p[i].b = 2.2;
    printf("A data[%d] %g, %g\n",i,p[i].a,p[i].b);
}
__global__ void MakeVectorForMoveKernel(int ngx,int ngy,point *p)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    p[i].a = ngx;
    p[i].b = ngy;
}
int test(void)
{
    // set number of points 
    int numPoints    = 1;
    int gpuBlockSize = 1;
    int pointSize    = sizeof(point);
    int numBytes     = numPoints * pointSize;
    int gpuGridSize  = numPoints / gpuBlockSize;
    point *cpuPointArray;
    point *gpuPointArray;
        // allocate memory
    cpuPointArray = (point*)malloc(numBytes);
    gpuPointArray = (point*)malloc(numBytes);
    cudaMalloc((void**)&gpuPointArray, numBytes);

    // launch kernel
    MakeVectorForMoveKernel<<<gpuGridSize,gpuBlockSize>>>(ngx,ngy,gpuPointArray);
    testKernel<<<gpuGridSize,gpuBlockSize>>>(gpuPointArray);

    // retrieve the results
    checkCudaErrors(cudaMemcpy(cpuPointArray, gpuPointArray, numBytes, cudaMemcpyDeviceToHost));
    printf("testKernel results:\n");
    for(int i = 0; i < numPoints; ++i)
    {
        printf("point.a: %g, point.b: %g\n",cpuPointArray[i].a,cpuPointArray[i].b);
    }
        // deallocate memory
    free(cpuPointArray);
    cudaFree(gpuPointArray);

    return 0;
}
