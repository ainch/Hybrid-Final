#include "cuda_Fluid.cuh"

void Set_Fluid_cuda(){
    int isp,i,j,GID;
    int size,grid,block,mingrid; 
    Conti_Flag = 0;
    
    checkCudaErrors(cudaMalloc((void**)&dev_FG, nfsp * sizeof(Fluid)));
    checkCudaErrors(cudaMemcpy(dev_FG, FG, nfsp * sizeof(Fluid), cudaMemcpyHostToDevice));

    // Density, D, nfsp * Csize
    checkCudaErrors(cudaMalloc((void**)&dev_C_F, nfsp * Csize *  sizeof(GFC)));
    checkCudaErrors(cudaMemcpy(dev_C_F, Host_C_F, nfsp * Csize *  sizeof(GFC), cudaMemcpyHostToDevice));

    // Flux_x, Flux_y,  Gummel_ax,bx,ay,by : [nfsp] * [Gsize]
    Host_G_F = (GFG *) malloc(nfsp * Gsize * sizeof(GFG)); //__Global_Fluid_Gsize_Data
    for(isp=0;isp<nfsp;isp++){
        for(i=0;i<Gsize;i++){
            GID = isp*Gsize + i;
            Host_G_F[GID].Flux_x = 0.0f;
            Host_G_F[GID].Gummel_ax = 0.0f;
            Host_G_F[GID].Gummel_bx = 0.0f;
            Host_G_F[GID].Flux_y = 0.0f;
            Host_G_F[GID].Gummel_ay = 0.0f;
            Host_G_F[GID].Gummel_by = 0.0f;
        }
    }
    checkCudaErrors(cudaMalloc((void**)&dev_G_F, nfsp * Gsize *  sizeof(GFG)));
    checkCudaErrors(cudaMemcpy(dev_G_F, Host_G_F, nfsp * Gsize *  sizeof(GFG), cudaMemcpyHostToDevice));

    // Calculate Diffusion coefficient
    if(MainGas == ARGON) Cal_D_Argon<<<CONTI_GRID,CONTI_BLOCK>>>(nfsp,ncy,Csize,Total_Pressure,dev_CvecSet,dev_GvecSet,dev_C_F);
    if(MainGas == OXYGEN) Cal_D_Oxygen<<<CONTI_GRID,CONTI_BLOCK>>>(nfsp,ncy,Csize,Total_Pressure,dev_CvecSet,dev_GvecSet,dev_C_F);
    if(MainGas == ARO2) Cal_D_ArO2<<<CONTI_GRID,CONTI_BLOCK>>>(nfsp,ncy,Csize,Total_Pressure,dev_CvecSet,dev_GvecSet,dev_C_F);
    checkCudaErrors(cudaMemcpy(Host_C_F, dev_C_F, nfsp * Csize *  sizeof(GFC), cudaMemcpyDeviceToHost));

    // Region check
    Conti_xnum = Cal_XRegion_check();
    Conti_ynum = Cal_YRegion_check();
    // number of tri cal = nfsp * Conti_xnum + nfsp * Conti_ynum
    // tridiag data size = nfsp * (3 * ncx * Conti_ynum)
    // tridiag data size = nfsp * (3 * ncy * Conti_xnum)
    size = nfsp * Conti_xnum;
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)Cal_D_Argon,0,size); 
    grid = (size + block - 1) / block;
    CONTIx_GRID = dim3(grid, 1, 1);
    CONTIx_BLOCK = dim3(block, 1, 1);
    size = nfsp * Conti_ynum;
    cudaOccupancyMaxPotentialBlockSize(&mingrid,&block,(void*)Cal_D_Argon,0,size); 
    grid = (size + block - 1) / block;
    CONTIy_GRID = dim3(grid, 1, 1);
    CONTIy_BLOCK = dim3(block, 1, 1);
    //
    Conti_x = (Con_RegionX *) malloc(nfsp * Conti_xnum * sizeof(Con_RegionX));
    Conti_y = (Con_RegionY *) malloc(nfsp * Conti_ynum * sizeof(Con_RegionY));
    for(isp=0;isp<nfsp;isp++){
        Set_Con_Region(isp, Conti_x, Conti_y);
        Set_Con_Boundary(isp, Conti_x, Conti_y);
    } 
    checkCudaErrors(cudaMalloc((void**)&dev_Conti_x, nfsp * Conti_xnum *  sizeof(Con_RegionX)));
    checkCudaErrors(cudaMemcpy(dev_Conti_x, Conti_x, nfsp * Conti_xnum *  sizeof(Con_RegionX), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMalloc((void**)&dev_Conti_y, nfsp * Conti_ynum *  sizeof(Con_RegionY)));
    checkCudaErrors(cudaMemcpy(dev_Conti_y, Conti_y, nfsp * Conti_ynum *  sizeof(Con_RegionY), cudaMemcpyHostToDevice));

    // Calculate Gummel Coef 
    calculate_gummel_coef_x<<<CONTIx_GRID,CONTIx_BLOCK>>>(nfsp, Conti_xnum, ngy, Gsize, Csize, dx, dev_Conti_x, dev_FG, dev_C_F, dev_G_F);
    calculate_gummel_coef_y<<<CONTIy_GRID,CONTIy_BLOCK>>>(nfsp, Conti_ynum, ngy, Gsize, Csize, dy, dev_Conti_y, dev_FG, dev_C_F, dev_G_F);

    // Calculate Tridiagonal A matrix 
    // tridiag size 3 * ncx * Conti_xnum
    // tridiag size 3 * ncy * Conti_ynum


    cudaDeviceSynchronize();
    exit(1);
}
__global__ void calculate_gummel_coef_x(int nfsp, int size, int ngy, int Gsize,int Csize, float dx, Con_RegionX *val, Fluid *info, GFC *Cdata, GFG *Gdata){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=nfsp*size) return;
    int i,isp,ID;
    ID = TID%size;
    isp = TID/size;
    int x1,x2,yy;
    float aa=0,bb;
    
    x1=val[TID].x1;
	x2=val[TID].x2;
	yy=val[TID].yy;
    Gdata[isp*Gsize + x1*ngy+yy].Gummel_bx =- val[TID].fg1*(-0.25*info[isp].Vel);
    Gdata[isp*Gsize + x1*ngy+yy].Gummel_ax =  0;
	for(i=x1+1;i<=x2;i++){
		bb=(Cdata[isp*Csize + (i-1)*(ngy-1)+yy].D + Cdata[isp*Csize + i*(ngy-1)+yy].D)/2;
		Gdata[isp*Gsize + i*ngy+yy].Gummel_ax=bb/dx;
		Gdata[isp*Gsize + i*ngy+yy].Gummel_bx=bb/dx;
        printf("%g ",bb/dx);
	}
    Gdata[isp*Gsize + (x2+1)*ngy+yy].Gummel_bx = 0;
    Gdata[isp*Gsize + (x2+1)*ngy+yy].Gummel_ax =- val[TID].fg2*(0.25*info[isp].Vel);

}
__global__ void calculate_gummel_coef_y(int nfsp, int size, int ngy, int Gsize,int Csize, float dy, Con_RegionY *val, Fluid *info, GFC *Cdata, GFG *Gdata){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=nfsp*size) return;
    int j,isp,ID;
    ID = TID%size;
    isp = TID/size;
    int xx,y1,y2;
    float aa=0,bb;

    xx=val[TID].xx;
	y1=val[TID].y1;
	y2=val[TID].y2;
    Gdata[isp*Gsize + xx*ngy+y1].Gummel_by =- val[TID].fg1*(-0.25*info[isp].Vel);
    Gdata[isp*Gsize + xx*ngy+y1].Gummel_ay =  0;
	for(j=y1+1;j<=y2;j++){
		bb=(Cdata[isp*Csize + xx*(ngy-1)+j-1].D + Cdata[isp*Csize + xx*(ngy-1)+j].D)/2;
		Gdata[isp*Gsize + xx*ngy+j].Gummel_ay=bb/dy;
		Gdata[isp*Gsize + xx*ngy+j].Gummel_by=bb/dy;
	}
    Gdata[isp*Gsize + xx*ngy+y2+1].Gummel_by = 0;
    Gdata[isp*Gsize + xx*ngy+y2+1].Gummel_ay =- val[TID].fg2*(0.25*info[isp].Vel);

}
void Set_Con_Boundary(int isp, Con_RegionX *Cx,Con_RegionY *Cy){
    int i, j ,k,GID,GID1,GID2,GID3;
	int x1, x2, y1, y2, xx, yy;
	for(i=0;i<Conti_xnum;i++) {
		x1=Cx[i].x1;
		x2=Cx[i].x2;
		yy=Cx[i].yy;
        GID = x1*ngy + yy;
        GID1 = x1*ngy + yy+1;
        GID2 = (x2+1)*ngy + yy;
        GID3 = (x2+1)*ngy + yy + 1;
		if(vec_G[GID].Boundary==NEUMANN || vec_G[GID1].Boundary==NEUMANN) {
			if(x1==0) 
                Cx[isp*Conti_xnum+i].fg1=0;
			else 
                Cx[isp*Conti_xnum+i].fg1=FG[isp].Gamma1;
		}
		else if(vec_G[GID].Boundary==DIELECTRIC) {
			Cx[isp*Conti_xnum+i].fg1=FG[isp].Gamma1;
		}
		else if(vec_G[GID].Boundary==CONDUCTOR) {
			Cx[isp*Conti_xnum+i].fg1=FG[isp].Gamma1;
		}
		else Cx[isp*Conti_xnum+i].fg1=FG[isp].Gamma1;

		if(vec_G[GID2].Boundary==NEUMANN || vec_G[GID3].Boundary==NEUMANN) {
			if(x2==ncx-1) Cx[isp*Conti_xnum+i].fg2=0;
			else Cx[isp*Conti_xnum+i].fg2=FG[isp].Gamma1;
		}
		else if(vec_G[GID2].Boundary==DIELECTRIC) {
			Cx[isp*Conti_xnum+i].fg2=FG[isp].Gamma1;
		}
		else if(vec_G[GID2].Boundary==CONDUCTOR) {
			Cx[isp*Conti_xnum+i].fg2=FG[isp].Gamma1;
		}
		else Cx[isp*Conti_xnum+i].fg2=FG[isp].Gamma1;

	}
	for(i=0;i<Conti_ynum;i++) {
        xx=Cy[i].xx;
		y1=Cy[i].y1;
		y2=Cy[i].y2;
        GID = xx*ngy + y1;
        GID1 = (xx+1)*ngy + y1;
        GID2 = xx*ngy + y2 + 1;
        GID3 = (xx+1)*ngy + y2+1;
		if(vec_G[GID].Boundary==NEUMANN || vec_G[GID1].Boundary==NEUMANN) {
			if(y1==0) Cy[isp*Conti_ynum+i].fg1=0;
			else Cy[isp*Conti_ynum+i].fg1=FG[isp].Gamma1;
		}
		else if(vec_G[GID].Boundary==DIELECTRIC) {
			Cy[isp*Conti_ynum+i].fg1=FG[isp].Gamma1;
		}
		else if(vec_G[GID].Boundary==CONDUCTOR) {
			Cy[isp*Conti_ynum+i].fg1=FG[isp].Gamma1;
		}
		else Cy[isp*Conti_ynum+i].fg1=FG[isp].Gamma1;
		if(vec_G[GID2].Boundary==NEUMANN || vec_G[GID3].Boundary==NEUMANN) {
			if(y2==ncy-1) Cy[isp*Conti_xnum+i].fg2=0;
			else Cy[isp*Conti_ynum+i].fg2=FG[isp].Gamma1;
		}
		else if(vec_G[GID2].Boundary==DIELECTRIC) {
			Cy[isp*Conti_ynum+i].fg2=FG[isp].Gamma1;
		}
		else if(vec_G[GID2].Boundary==CONDUCTOR) {
			Cy[isp*Conti_ynum+i].fg2=FG[isp].Gamma1;
		}
		else Cy[isp*Conti_ynum+i].fg2=FG[isp].Gamma1;
	}
}
void Set_Con_Region(int isp, Con_RegionX *Cx,Con_RegionY *Cy){
    int i, j ,k,add,CID,CID1,CID2;
	int x1, x2, y1, y2, xx, yy;
	k=0;
	for(j=0;j<ncy;j++) {
		for(i=0;i<ncx;i++) {
            CID = i*ncy + j;
            CID1 = (i-1)*ncy + j;
            CID2 = (i+1)*ncy + j;
			if(vec_C[CID].PlasmaRegion == 1) { // plasma
				if(i == 0) { //left side wallCycle of check
					Cx[isp * Conti_xnum + k].x1 = i;
					Cx[isp * Conti_xnum + k].yy = j;
				}
				else if(vec_C[CID1].PlasmaRegion != 1) { // if left side is not plasma region
					Cx[isp * Conti_xnum + k].x1 = i;
					Cx[isp * Conti_xnum + k].yy = j;
				}
				if(i == ncx-1) { //right side wall
					Cx[isp * Conti_xnum + k].x2 = i;
					k++;
				}
				else if(vec_C[CID2].PlasmaRegion != 1) { // if Right side is not plasma region
					Cx[isp * Conti_xnum + k].x2 = i;
					k++;
				}
			}
		}
	}
	k=0;
	for(i=0;i<ncx;i++) {
		for(j=0;j<ncy;j++) {
            CID = i*ncy + j;
            CID1 = i*ncy + j-1;
            CID2 = i*ncy + j+1;
			if(vec_C[CID].PlasmaRegion == 1) {
				if(j == 0) {
					Cy[isp * Conti_ynum + k].xx = i;
					Cy[isp * Conti_ynum + k].y1 = j;
				}
				else if(vec_C[CID1].PlasmaRegion != 1) {
					Cy[isp * Conti_ynum + k].xx = i;
					Cy[isp * Conti_ynum + k].y1 = j;
				}

				if(j == ncy-1) {
					Cy[isp * Conti_ynum + k].y2 = j;
					k++;
				}
				else if(vec_C[CID2].PlasmaRegion != 1) {
					Cy[isp * Conti_ynum + k].y2 = j;
					k++;
				}
			}
		}
	}
}
int Cal_XRegion_check(){
    int i,j,k;
    k=0;
    for(j=0;j<ncy;j++){
        for(i=0;i<ncx;i++){
            if(vec_C[i*ncy+j].PlasmaRegion == 1){
                if(i == ncx-1) {
                    k++;
                }else if(vec_C[(i+1)*ncy+j].PlasmaRegion != 1) {
                    k++;
                }
            }
        }
    }
    return k;
}
int Cal_YRegion_check(){
    int i,j,k;
    k=0;
    for(i=0;i<ncx;i++){
        for(j=0;j<ncy;j++){
            if(vec_C[i*ncy+j].PlasmaRegion == 1){
                if(j == ncy-1) {
                    k++;
                }else if(vec_C[i*ncy+j+1].PlasmaRegion != 1) {
                    k++;
                }
            }
        }
    }
    return k;
}
__global__ void Cal_D_Argon(int nfsp, int ncy, int Csize, float press, GCA *vec_C, GGA *vec_G, GFC *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=nfsp*Csize) return;
    int ID = TID%Csize;
    int xi = ID/ncy;
    int yi = ID%ncy;
    int GID;
    float T = 0.0f;
    float epsi_k = 124.0f; 	// epsilon/K from Viscosity
	float sigmasq = 11.682724f; 	// r0, A. from Viscosity (3.418)^2
    float omega;

    GID = xi * (ncy+1) + yi;
    T = 0.25*(vec_G[GID].Temp + vec_G[GID+1].Temp + vec_G[GID+ncy+1].Temp + vec_G[GID+ncy+2].Temp);
    omega = Ar_meta_omega(T/epsi_k);
    data[TID].D = 1.997279e-4*sqrt(T/39.948)*T/(press*sigmasq*omega)*20;
    //printf("D[%d] = %g\n",TID,data[TID].D);
}
__global__ void Cal_D_Oxygen(int nfsp, int ncy, int Csize, float press, GCA *vec_C, GGA *vec_G, GFC *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=nfsp*Csize) return;
    int ID = TID%Csize;
    int isp = TID/Csize;
    int xi = ID/ncy;
    int yi = ID%ncy;
    int GID;
    float T = 0.0f;
	float sigmasq1 = 0.10892f; 	// r0, A. from Viscosity 1/(3.03)^2  for OP  OD
    float sigmasq2 = 0.08650f; 	// r0, A. from Viscosity 1/(3.40)^2  for O2A O2B

    GID = xi * (ncy+1) + yi;
    T = 0.25*(vec_G[GID].Temp + vec_G[GID+1].Temp + vec_G[GID+ncy+1].Temp + vec_G[GID+ncy+2].Temp);
    if(isp <= 1){
        // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
        data[TID].D = 2.63e-3*sqrt(T*T*T*24/32/32)/press*sigmasq2;
        //if(ID == 510) printf("D2[%d] = %g\n",TID,data[TID].D);
    }else{
        // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
        data[TID].D = 2.63e-3*sqrt(T*T*T*24/16/16)/press*sigmasq1;
        //if(ID == 510) printf("D1[%d] = %g\n",TID,data[TID].D);
    }
}
__global__ void Cal_D_ArO2(int nfsp, int ncy, int Csize, float press, GCA *vec_C, GGA *vec_G, GFC *data){
    int TID = threadIdx.x + blockIdx.x * blockDim.x;
    if(TID>=nfsp*Csize) return;
    int ID = TID%Csize;
    int isp = TID/Csize;
    int xi = ID/ncy;
    int yi = ID%ncy;
    int GID;
    float T = 0.0f;
    //float epsi_k = 124.0f; 	// epsilon/K from Viscosity
	//float sigmasq = 11.682724f; // r0, A. from Viscosity (3.418)^2  for AR meta
    float sigmasq0 = 0.08559f;  // r0, A. from Viscosity 1/(3.418)^2 for AR meta
	float sigmasq1 = 0.10892f; 	// r0, A. from Viscosity 1/(3.03)^2  for OP  OD
    float sigmasq2 = 0.08650f; 	// r0, A. from Viscosity 1/(3.40)^2  for O2A O2B
    //float omega;

    GID = xi * (ncy+1) + yi;
    T = 0.25*(vec_G[GID].Temp + vec_G[GID+1].Temp + vec_G[GID+ncy+1].Temp + vec_G[GID+ncy+2].Temp);
    if(isp == 0){
        //omega = Ar_meta_omega(T/epsi_k);
        //data[TID].D = 1.997279e-4*sqrt(T/39.948)*T/(press*sigmasq*omega)*20;
        data[TID].D = 2.63e-3*sqrt(T*T*T*24/39.948/39.948)/press*sigmasq0;
        //if(ID == 510) printf("D0[%d] = %g\n",TID,data[TID].D);
    }else if(isp == 1){
        // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
        data[TID].D = 2.63e-3*sqrt(T*T*T*24/32/32)/press*sigmasq2;
        //if(ID == 510) printf("D1[%d] = %g\n",TID,data[TID].D);
    }else if(isp == 2){
        // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
        data[TID].D = 2.63e-3*sqrt(T*T*T*24/32/32)/press*sigmasq2;
        //if(ID == 510) printf("D2[%d] = %g\n",TID,data[TID].D);
    }else if(isp == 3){
        // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
        data[TID].D = 2.63e-3*sqrt(T*T*T*24/16/16)/press*sigmasq1;
        //if(ID == 510) printf("D3[%d] = %g\n",TID,data[TID].D);
    }else if(isp == 4){
        // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
        data[TID].D = 2.63e-3*sqrt(T*T*T*24/16/16)/press*sigmasq1;
        //if(ID == 510) printf("D4[%d] = %g\n",TID,data[TID].D);
    }
}
__device__ float Ar_meta_omega(float value){
    int   i;
	float tstar[27] = {0.250, 0.267, 0.286, 0.308, 0.333, 0.364, 0.400, 0.444,
	                   0.500, 0.571, 0.667, 0.800, 1.000, 1.330, 2.000, 2.500,
	                   2.860, 3.330, 4.000, 5.000, 6.670, 10.00, 12.50, 16.67,
	                   25.00, 50.00, 100.0},
	      omega[27] = {1.8714, 1.8344, 1.7962, 1.7560, 1.7142, 1.6702, 1.6238,
	                   1.5748, 1.5228, 1.4676, 1.4090, 1.3466, 1.2802, 1.2096,
	                   1.1348, 1.1042, 1.0888, 1.0738, 1.0592, 1.0450, 1.0318,
	                   1.0198, 1.0154, 1.0112, 1.0072, 1.0034, 1.000};
	float temp1;
    float temp2 = DBL_MIN;
	if      (value<tstar[0])  temp2=omega[0];
	else if (value>tstar[26]) temp2=omega[26];
	else
		for (i=1;i<27;i++)
			if (value<tstar[i]) {
				temp1=(value-tstar[i-1])/(tstar[i]-tstar[i-1]);
				temp2=omega[i]*temp1+omega[i-1]*(1-temp1);
			}
	return temp2;
}