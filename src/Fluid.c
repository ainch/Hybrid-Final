#include "Fluid.h"

void Fluid_Setting(){
    int i,j,isp,CID,GID;
    float xx,yy,wv;
    printf("Fluid Setting\n");
    Fluid_sp = (GFC *) malloc(nfsp * sizeof(GFC));
    if(FG[0].Loadtype == 0)       printf("\tFluid Load Type : UNIFORM\n");
    else if(FG[0].Loadtype == 1)  printf("\tFluid Load Type : EXPONETIAL\n");
    else if(FG[0].Loadtype == 2)  printf("\tFluid Load Type : COSINE\n");
    else if(FG[0].Loadtype == 4){
        printf("\tFluid Load Type : SmartLoad\n");
        printf("\tError : not yet!\n");
        exit(1);
    }
    // Density and Source Setting for GPU
    Fluid_Den = (GFG *) malloc(nfsp * Gsize * sizeof(GFG));
    Fluid_Src = (GFG *) malloc(nfsp * Gsize * sizeof(GFG));
    // Density setting
    for(isp=0;isp<nfsp;isp++){
        FG[isp].ave_Den = 0;
        Fluid_sp[isp].D = MFMalloc(ncx,ncy);
        Fluid_sp[isp].den = MFMalloc(ncx,ncy);
        Fluid_sp[isp].Source = MFMalloc(ncx,ncy);
        Fluid_sp[isp].gummel_ax = MFMalloc(ngx,ncy); 
        MFInit(Fluid_sp[isp].gummel_ax,0.0f,ngx,ncy); 
        Fluid_sp[isp].gummel_bx = MFMalloc(ngx,ncy);
        MFInit(Fluid_sp[isp].gummel_bx,0.0f,ngx,ncy); 
        Fluid_sp[isp].flux_x = MFMalloc(ngx,ncy);
        MFInit(Fluid_sp[isp].flux_x,0.0f,ngx,ncy); 
        Fluid_sp[isp].gummel_ay = MFMalloc(ncx,ngy);
        MFInit(Fluid_sp[isp].gummel_ax,0.0f,ngx,ncy); 
        Fluid_sp[isp].gummel_by = MFMalloc(ncx,ngy);
        MFInit(Fluid_sp[isp].gummel_bx,0.0f,ngx,ncy); 
        Fluid_sp[isp].flux_y = MFMalloc(ncx,ngy);
        MFInit(Fluid_sp[isp].flux_x,0.0f,ngx,ncy); 
        for(i=0;i<ncx;i++){
		    for(j=0;j<ncy;j++){
                Fluid_sp[isp].D[i][j] = 0.0f;
                Fluid_sp[isp].Source[i][j] = 0.0f;
                //
                Fluid_sp[isp].den[i][j] = 0.0f;
                if(DumpFlag == 0){
                    xx = (float) i * dx;
			        yy = (float) j * dy;
                    CID = i * ncy + j;
                    if(vec_C[CID].PlasmaRegion != 0){
                        if(FG[0].Loadtype == 0){ // UNIFORM
                            Fluid_sp[isp].den[i][j] = FG[isp].InitDens;
                        //printf("DEN[%d]= %g\n",CID,Fluid_sp[CID].den);
                        }else if(FG[0].Loadtype == 1){// EXPONETIAL
           	                wv = exp(-1 * ((xx - FG[0].x_center)/FG[0].x_fall)*((xx - FG[0].x_center)/FG[0].x_fall))
                                *exp(-1 * ((yy - FG[0].y_center)/FG[0].y_fall)*((yy - FG[0].y_center)/FG[0].y_fall));
				            Fluid_sp[isp].den[i][j] = FG[isp].InitDens * wv;
                            //printf("DEN[%d]= %g\n",CID,Fluid_sp[CID].den);
                        }else if(FG[0].Loadtype == 2){// COSINE
                            wv = fabs(cos((xx - FG[0].x_center)*M_PI/2/FG[0].x_fall)*cos((yy - FG[0].y_center)*M_PI/2/FG[0].y_fall));
				            Fluid_sp[isp].den[i][j] = FG[isp].InitDens * wv;
                        }
                        FG[isp].ave_Den += Fluid_sp[isp].den[i][j];
                    }
                }
                GID = isp * Gsize + i * ngy + j;
                Fluid_Den[GID].n = 0.0f;
                Fluid_Src[GID].n = 0.0f;
                if(i == ncx-1 || j == ncy-1){
                    GID = isp * Gsize + i+1 * ngy + j+1;
                    Fluid_Den[GID].n = 0.0f;
                    Fluid_Src[GID].n = 0.0f;
                }
            }
        }
    }
    // Plasma Reigon
    Conti_Flag = 0;
    // Region check
    Conti_xnum = Cal_XRegion_check();
    printf("\tX direction Number of line = %d (<ngy(%d))\n",Conti_xnum,ngy);
    Conti_ynum = Cal_YRegion_check();
    printf("\tY direction Number of line = %d (<ngx(%d))\n",Conti_ynum,ngx);
    Conti_x = (Con_RegionX *) malloc(nfsp * sizeof(Con_RegionX));
    Conti_y = (Con_RegionY *) malloc(nfsp * sizeof(Con_RegionY));
    for(isp=0;isp<nfsp;isp++){
        Conti_x[isp].x1 = VIMalloc(Conti_xnum);
        Conti_x[isp].x2 = VIMalloc(Conti_xnum);
        Conti_x[isp].yy = VIMalloc(Conti_xnum);
        Conti_x[isp].fg1 = VFMalloc(Conti_xnum);
        Conti_x[isp].fg2 = VFMalloc(Conti_xnum);
        Conti_y[isp].xx = VIMalloc(Conti_ynum);
        Conti_y[isp].y1 = VIMalloc(Conti_ynum);
        Conti_y[isp].y2 = VIMalloc(Conti_ynum);
        Conti_y[isp].fg1 = VFMalloc(Conti_ynum);
        Conti_y[isp].fg2 = VFMalloc(Conti_ynum);
        Set_Con_Region(isp, Conti_x, Conti_y);
        Set_Con_Boundary(isp, Conti_x, Conti_y);
    }
}
void Set_Con_Boundary(int isp, Con_RegionX *Cx,Con_RegionY *Cy){
    int i, j ,k,GID,GID1,GID2,GID3;
	int x1, x2, y1, y2, xx, yy;
	for(i=0;i<Conti_xnum;i++) {
		x1=Cx[isp].x1[i];
		x2=Cx[isp].x2[i];
		yy=Cx[isp].yy[i];
        GID = x1*ngy + yy;
        GID1 = x1*ngy + yy+1;
        GID2 = (x2+1)*ngy + yy;
        GID3 = (x2+1)*ngy + yy + 1;
		if(vec_G[GID].Boundary==NEUMANN || vec_G[GID1].Boundary==NEUMANN) {
			if(x1==0) 
                Cx[isp].fg1[i] = 0.0f;
			else 
                Cx[isp].fg1[i] = FG[isp].Gamma1;
		}
		else if(vec_G[GID].Boundary==DIELECTRIC) {
			Cx[isp].fg1[i]=FG[isp].Gamma1;
		}
		else if(vec_G[GID].Boundary==CONDUCTOR) {
			Cx[isp].fg1[i]=FG[isp].Gamma1;
		}
		else Cx[isp].fg1[i]=FG[isp].Gamma1;

		if(vec_G[GID2].Boundary==NEUMANN || vec_G[GID3].Boundary==NEUMANN) {
			if(x2==ncx-1) Cx[isp].fg2[i]=0.0f;
			else Cx[isp].fg2[i]=FG[isp].Gamma1;
		}
		else if(vec_G[GID2].Boundary==DIELECTRIC) {
			Cx[isp].fg2[i]=FG[isp].Gamma1;
		}
		else if(vec_G[GID2].Boundary==CONDUCTOR) {
			Cx[isp].fg2[i]=FG[isp].Gamma1;
		}
		else Cx[isp].fg2[i]=FG[isp].Gamma1;

	}
	for(i=0;i<Conti_ynum;i++) {
        xx=Cy[isp].xx[i];
		y1=Cy[isp].y1[i];
		y2=Cy[isp].y2[i];
        GID = xx*ngy + y1;
        GID1 = (xx+1)*ngy + y1;
        GID2 = xx*ngy + y2 + 1;
        GID3 = (xx+1)*ngy + y2+1;
		if(vec_G[GID].Boundary==NEUMANN || vec_G[GID1].Boundary==NEUMANN) {
			if(y1==0) Cy[isp].fg1[i]=0;
			else Cy[isp].fg1[i]=FG[isp].Gamma1;
		}
		else if(vec_G[GID].Boundary==DIELECTRIC) {
			Cy[isp].fg1[i]=FG[isp].Gamma1;
		}
		else if(vec_G[GID].Boundary==CONDUCTOR) {
			Cy[isp].fg1[i]=FG[isp].Gamma1;
		}
		else Cy[isp].fg1[i]=FG[isp].Gamma1;
		if(vec_G[GID2].Boundary==NEUMANN || vec_G[GID3].Boundary==NEUMANN) {
			if(y2==ncy-1) Cy[isp].fg2[i]=0;
			else Cy[isp].fg2[i]=FG[isp].Gamma1;
		}
		else if(vec_G[GID2].Boundary==DIELECTRIC) {
			Cy[isp].fg2[i]=FG[isp].Gamma1;
		}
		else if(vec_G[GID2].Boundary==CONDUCTOR) {
			Cy[isp].fg2[i]=FG[isp].Gamma1;
		}
		else Cy[isp].fg2[i]=FG[isp].Gamma1;
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
					Cx[isp].x1[k] = i;
					Cx[isp].yy[k] = j;
				}
				else if(vec_C[CID1].PlasmaRegion != 1) { // if left side is not plasma region
					Cx[isp].x1[k] = i;
					Cx[isp].yy[k] = j;
				}
				if(i == ncx-1) { //right side wall
					Cx[isp].x2[k] = i;
					k++;
				}
				else if(vec_C[CID2].PlasmaRegion != 1) { // if Right side is not plasma region
					Cx[isp].x2[k] = i;
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
					Cy[isp].xx[k] = i;
					Cy[isp].y1[k]  = j;
				}
				else if(vec_C[CID1].PlasmaRegion != 1) {
					Cy[isp].xx[k]  = i;
					Cy[isp].y1[k]  = j;
				}

				if(j == ncy-1) {
					Cy[isp].y2[k]  = j;
					k++;
				}
				else if(vec_C[CID2].PlasmaRegion != 1) {
					Cy[isp].y2[k]  = j;
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
