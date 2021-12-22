#include "cuda_Fluid.cuh"

void Set_Fluid_cuda(){
    int isp,i,j,k,GID;
    int size,grid,block,mingrid; 

    checkCudaErrors(cudaMalloc((void**)&dev_FG, nfsp * sizeof(Fluid)));
    checkCudaErrors(cudaMemcpy(dev_FG, FG, nfsp * sizeof(Fluid), cudaMemcpyHostToDevice));

    Sync_Fluid_GFCtoGFG_forDen(Fluid_sp, Fluid_Den); 
    checkCudaErrors(cudaMalloc((void**)&dev_FG_Den, nfsp * Gsize * sizeof(GFG)));
    checkCudaErrors(cudaMemcpy(dev_FG_Den, Fluid_Den, nfsp * Gsize * sizeof(GFG), cudaMemcpyHostToDevice));

    checkCudaErrors(cudaMalloc((void**)&dev_FG_Src, nfsp * Gsize * sizeof(GFG)));
    checkCudaErrors(cudaMemcpy(dev_FG_Src, Fluid_Src, nfsp * Gsize * sizeof(GFG), cudaMemcpyHostToDevice));
   
    // Calculate Diffusion coefficient
    if(MainGas == ARGON) 
        Cal_D_Argon(vec_C,vec_G,Fluid_sp);
    else if(MainGas == OXYGEN) 
        Cal_D_Oxygen(vec_C,vec_G,Fluid_sp);
    else if(MainGas == ARO2) 
        Cal_D_ArO2(vec_C,vec_G,Fluid_sp);

    // Gummel coefficient
    calculate_gummel_coef_x();
    calculate_gummel_coef_y();

    for(isp=0;isp<nfsp;isp++) {
        Calculate_Flux_x(isp);
        Calculate_Flux_y(isp);
    }
}
void Solve_Continuity_eqn(){
    int isp;
    for(isp=0;isp<nfsp;isp++) {
		Solve_Density_x(isp);
		Calculate_Flux_x(isp);
		Solve_Density_y(isp);
		Calculate_Flux_y(isp);
	}
}
void Solve_Continuity_eqn_check(){
    int isp;
    for(isp=0;isp<nfsp;isp++) {
        if(FG[isp].CSS_Flag){
		    Solve_Density_x(isp);
		    Calculate_Flux_x(isp);
		    Solve_Density_y(isp);
		    Calculate_Flux_y(isp);
        }
	}
}
void Calculate_Flux_x(int isp){
    int i, k, x1, x2, yy;
	float left_den, right_den;
    GFC *pointer;
    Con_RegionX *pointer2;
    pointer = &(Fluid_sp[isp]);
    pointer2 = &(Conti_x[isp]);

	for(k=0;k<Conti_xnum;k++){
		x1=pointer2->x1[k];
		x2=pointer2->x2[k];
		yy=pointer2->yy[k];
		pointer->flux_x[x1][yy] = -Conti_x[isp].fg1[k] * pointer->gummel_bx[x1][yy] * pointer->den[x1][yy];
		for(i=x1+1;i<=x2;i++) {
			left_den = pointer->den[i-1][yy];
			right_den = pointer->den[i][yy];
			pointer->flux_x[i][yy] = pointer->gummel_ax[i][yy]*left_den - pointer->gummel_bx[i][yy]*right_den;
		}
		pointer->flux_x[x2+1][yy] = Conti_x[isp].fg2[k]*pointer->gummel_ax[x2+1][yy]*pointer->den[x2][yy];
	}
}
void Calculate_Flux_y(int isp){
    int j, k, y1, y2, xx;
	float lower_den, upper_den;
    GFC *pointer;
    Con_RegionY * pointer2;
    pointer = &(Fluid_sp[isp]);
    pointer2 = &(Conti_y[isp]);

	for(k=0;k<Conti_ynum;k++){
		xx=pointer2->xx[k];
		y1=pointer2->y1[k];
		y2=pointer2->y2[k];
		pointer->flux_y[xx][y1] = -Conti_y[isp].fg1[k]*pointer->gummel_by[xx][y1]*pointer->den[xx][y1];
		for(j=y1+1;j<=y2;j++) {
			lower_den = pointer->den[xx][j-1];
			upper_den = pointer->den[xx][j];
			pointer->flux_y[xx][j] = pointer->gummel_ay[xx][j]*lower_den - pointer->gummel_by[xx][j]*upper_den;
		}
		pointer->flux_y[xx][y2+1] = Conti_y[isp].fg2[k]*pointer->gummel_ay[xx][y2+1]*pointer->den[xx][y2];
	}
}
void Solve_Density_x(int isp){
    int k,i,x1,x2,yy,num,xx;
    float src, ga, gb, gc;
    float subd[ngx], diag[ngx], superd[ngx], rhs[ngx], density[ngx], gam[ngx];
    float del_f;
    GFC *pointer;
    Con_RegionX *pointer2;
    pointer = &(Fluid_sp[isp]);
    pointer2 = &(Conti_x[isp]);
    del_f = dtc/dx;
    for(k=0;k<Conti_xnum;k++){
        x1=pointer2->x1[k];
		x2=pointer2->x2[k];
		yy=pointer2->yy[k];
		num = x2-x1+1;
		for(i=0;i<num;i++) {
			xx = i + x1;
            ga =  -pointer->gummel_ax[xx][yy]*del_f;
            gb =  1+(pointer->gummel_ax[xx+1][yy]+pointer->gummel_bx[xx][yy])*del_f;
            gc =  -pointer->gummel_bx[xx+1][yy]*del_f;
            src = pointer->den[xx][yy] + pointer->Source[xx][yy]*dtc - (pointer->flux_y[xx][yy+1]-pointer->flux_y[xx][yy])*dtc/dy;
            subd[i] =  ga;
			diag[i] = gb;
			superd[i] = gc;
			rhs[i] =  src;
		}
        tridiag(subd,diag,superd,rhs,density,gam,0,num-1);
        for(i=0;i<num;i++) {
			if(vec_C[(i+x1)*ncy+yy].PlasmaRegion == 1 && density[i] > 1e9){
				pointer->den[i+x1][yy] = density[i];
            }else{
				pointer->den[i+x1][yy] = 0.0f;
			}
		}
    }
}
void Solve_Density_y(int isp){
    int k,j,xx,y1,y2,num,yy;
    float src, ga, gb, gc;
    float subd[ngy], diag[ngy], superd[ngy], rhs[ngy], density[ngy], gam[ngy];
    float del_f;
    GFC *pointer;
    Con_RegionY * pointer2;
    pointer = &(Fluid_sp[isp]);
    pointer2 = &(Conti_y[isp]);
    del_f = dtc/dy;
    for(k=0;k<Conti_ynum;k++){
		xx=pointer2->xx[k];
		y1=pointer2->y1[k];
		y2=pointer2->y2[k];
		num = y2-y1+1;
		for(j=0;j<num;j++) {
			yy = j + y1;
            ga =  -pointer->gummel_ay[xx][yy]*del_f;
            gb =  1+(pointer->gummel_ay[xx][yy+1]+pointer->gummel_by[xx][yy])*del_f;
            gc =  -pointer->gummel_by[xx][yy+1]*del_f;
            src = pointer->den[xx][yy] + pointer->Source[xx][yy]*dtc - (pointer->flux_x[xx+1][yy]-pointer->flux_x[xx][yy])*dtc/dx;
            subd[j] = ga;
			diag[j] = gb;
			superd[j] = gc;
			rhs[j] = src;
		}
		tridiag(subd,diag,superd,rhs,density,gam,0,num-1);
        for(j=0;j<num;j++) {
			if(vec_C[xx*ncy+j+y1].PlasmaRegion == 1 && density[j] > 1e9){
				pointer->den[xx][j+y1]=density[j];
            }else{
				pointer->den[xx][j+y1]=0.0f;
			}
		}
	}
}
int tridiag(float *a, float *b, float *c, float *d, float *x, float *gam, int min, int max)
{
	int j;
	float bet;

	x[min]=d[min]/(bet=b[min]);
	for(j=min+1;j<=max;j++) {
		gam[j-1]=c[j-1]/bet;
		bet=b[j]-a[j]*c[j-1]/bet;
		x[j]=(d[j]-a[j]*x[j-1])/bet;
		if(bet==0) return 0;
	}

	for(j=max-1;j>=min;j--) x[j] -= gam[j]*x[j+1];

	return 1;
}
void calculate_gummel_coef_x(){
    int isp,k,i;
    int x1,x2,yy;
    for(isp=0;isp<nfsp;isp++){
        for(k=0;k<Conti_xnum;k++){
		    x1=Conti_x[isp].x1[k];
		    x2=Conti_x[isp].x2[k];
		    yy=Conti_x[isp].yy[k];
		    Fluid_sp[isp].gummel_bx[x1][yy]=-Conti_x[isp].fg1[k] * (-0.25*sqrt(vec_G[x1*ngy+yy].Temp*1.38e-23/FG[isp].mass));
		    Fluid_sp[isp].gummel_ax[x1][yy]=0;
		    for(i=x1+1;i<=x2;i++){
			    Fluid_sp[isp].gummel_ax[i][yy] = (Fluid_sp[isp].D[i-1][yy] + Fluid_sp[isp].D[i-1][yy])/2/dx;
			    Fluid_sp[isp].gummel_bx[i][yy] = Fluid_sp[isp].gummel_ax[i][yy];
		    }
		    Fluid_sp[isp].gummel_ax[x2+1][yy]=Conti_x[isp].fg2[k] * (0.25*sqrt(vec_G[x2*ngy+yy].Temp*1.38e-23/FG[isp].mass));
		    Fluid_sp[isp].gummel_bx[x2+1][yy]=0;
	    }
    }
}
void calculate_gummel_coef_y(){
    int isp,k,j;
    int xx,y1,y2;
    float aa=0,bb;
    for(isp=0;isp<nfsp;isp++){
        for(k=0;k<Conti_ynum;k++){
		    xx=Conti_y[isp].xx[k];
		    y1=Conti_y[isp].y1[k];
		    y2=Conti_y[isp].y2[k];
            Fluid_sp[isp].gummel_by[xx][y1]=-Conti_y[isp].fg1[k] * (-0.25*sqrt(vec_G[xx*ngy+y1].Temp*1.38e-23/FG[isp].mass));
		    Fluid_sp[isp].gummel_ay[xx][y1]=0;
		    for(j=y1+1;j<=y2;j++){
			    Fluid_sp[isp].gummel_ay[xx][j] = (Fluid_sp[isp].D[xx][j-1] + Fluid_sp[isp].D[xx][j])/2/dy;
			    Fluid_sp[isp].gummel_by[xx][j] = Fluid_sp[isp].gummel_ay[xx][j];
		    }
		    Fluid_sp[isp].gummel_ay[xx][y2+1]=Conti_y[isp].fg2[k] * (0.25*sqrt(vec_G[xx*ngy+y2].Temp*1.38e-23/FG[isp].mass));
		    Fluid_sp[isp].gummel_by[xx][y2+1]=0;
        }
    }
}
void Cal_D_Argon(GCA *vecC, GGA *vecG, GFC *data){
    int isp,i,j,CID,GID;
    float T = 0.0f;
    float epsi_k = 124.0f; 	// epsilon/K from Viscosity
	float sigmasq = 11.682724f; 	// r0, A. from Viscosity (3.418)^2
    float omega;
    for(isp=0;isp<nfsp;isp++){
        for(i=0;i<ncx;i++){
		    for(j=0;j<ncy;j++){
                CID = i * ncy + j;
                GID = i * ngy + j;
                T = 0.25*(vecG[GID].Temp + vecG[GID+1].Temp + vecG[GID+ngy].Temp + vecG[GID+ngy+1].Temp);
                omega = Ar_meta_omega(T/epsi_k);
                data[isp].D[i][j] = 1.997279e-4*sqrt(T/39.948)*T/(Total_Pressure*sigmasq*omega)*20;
            }
        }
    }
    //printf("D[%d] = %g\n",TID,data[TID].D);
}
void Cal_D_Oxygen(GCA *vecC, GGA *vecG, GFC *data){
    int isp,i,j,CID,GID;
    float T = 0.0f;
	float sigmasq1 = 0.10892f; 	// r0, A. from Viscosity 1/(3.03)^2  for OP  OD
    float sigmasq2 = 0.08650f; 	// r0, A. from Viscosity 1/(3.40)^2  for O2A O2B
    float omega;
    for(isp=0;isp<nfsp;isp++){
        for(i=0;i<ncx;i++){
		    for(j=0;j<ncy;j++){
                CID = i * ncy + j;
                GID = i * ngy + j;
                T = 0.25*(vecG[GID].Temp + vecG[GID+1].Temp + vecG[GID+ngy].Temp + vecG[GID+ngy+1].Temp);
                if(isp <= 1){
                    // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
                    data[isp].D[i][j] = 2.63e-3*sqrt(T*T*T*24/32/32)/Total_Pressure*sigmasq2;
                }else{
                     // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
                    data[isp].D[i][j] = 2.63e-3*sqrt(T*T*T*24/16/16)/Total_Pressure*sigmasq1;  
                }
            }
        }
    }
}
void Cal_D_ArO2(GCA *vecC, GGA *vecG, GFC *data){   
    int isp,i,j,CID,GID;
    float T = 0.0f;
	float sigmasq0 = 0.08559f;  // r0, A. from Viscosity 1/(3.418)^2 for AR meta
	float sigmasq1 = 0.10892f; 	// r0, A. from Viscosity 1/(3.03)^2  for OP  OD
    float sigmasq2 = 0.08650f; 	// r0, A. from Viscosity 1/(3.40)^2  for O2A O2B
    float omega;
    for(isp=0;isp<nfsp;isp++){
        for(i=0;i<ncx;i++){
		    for(j=0;j<ncy;j++){
                CID = i * ncy + j;
                GID = i * ngy + j;
                T = 0.25*(vecG[GID].Temp + vecG[GID+1].Temp + vecG[GID+ngy].Temp + vecG[GID+ngy+1].Temp);
                if(isp == 0){
                    //data[TID].D = 1.997279e-4*sqrt(T/39.948)*T/(press*sigmasq*omega)*20;
                    data[isp].D[i][j] = 2.63e-3*sqrt(T*T*T*24/39.948/39.948)/Total_Pressure*sigmasq0;
                }else if(isp == 1){
                    // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
                    data[isp].D[i][j] = 2.63e-3*sqrt(T*T*T*24/32/32)/Total_Pressure*sigmasq2;
                }else if(isp == 2){
                    // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
                    data[isp].D[i][j] = 2.63e-3*sqrt(T*T*T*24/32/32)/Total_Pressure*sigmasq2;
                }else if(isp == 3){
                    // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
                    data[isp].D[i][j] = 2.63e-3*sqrt(T*T*T*24/16/16)/Total_Pressure*sigmasq1;
                }else if(isp == 4){
                    // D = [(2.63X10^3)/P/(SIGMA)^2] * sqrt( T^3 *(Ma+Mb)/2/Ma/Mb);
                    data[isp].D[i][j] = 2.63e-3*sqrt(T*T*T*24/16/16)/Total_Pressure*sigmasq1;
                }
            }
        }
    }
}
float Ar_meta_omega(float value){
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
void Sync_Fluid_GFCtoGFG_forDen(GFC *A, GFG *B){ // [ncx][ncy] --> [Gsize], cpu DATA > GPU
    int isp,i,j,GID;
    float buf;
    for(isp=0;isp<nfsp;isp++){
        buf = 0.0f;
        for(i=0;i<ncx;i++){
		    for(j=0;j<ncy;j++){
                GID = isp * Gsize + i * ngy + j;
                B[GID].n = A[isp].den[i][j];
                if(i == ncx-1 || j == ncy-1){
                    GID = isp * Gsize + (i+1) * ngy + j+1;
                    B[GID].n = A[isp].den[i][j];
                }
                buf += A[isp].den[i][j];
            }
        }
        FG[isp].ave_Den = buf/Csize;
    }
}
void Sync_Fluid_GFGtoGFC_forSource(GFG *A, GFC *B){ // [nsp*Gsize+Gsize] --> [ncx][ncy], GPU data > CPU
    int isp,i,j,GID;
    for(isp=0;isp<nfsp;isp++){
        for(i=0;i<ncx;i++){
		    for(j=0;j<ncy;j++){
                GID = isp * Gsize + i * ngy + j;
                B[isp].Source[i][j] = A[GID].n;
            }
        }
    }
}