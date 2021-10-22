#include "cuda_Diagnostic.cuh"

void Diagnostic(){


}
void Set_Diagnostic_cuda(){
    // Host BUF VECTOR
    Host_G_buf = VFMalloc(Gsize);
    Host_C_buf = VFMalloc(Csize);
    VFInit(Host_G_buf,0.0,Gsize);
    VFInit(Host_C_buf,0.0,Csize);

}