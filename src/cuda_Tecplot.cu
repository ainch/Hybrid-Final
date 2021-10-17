#include "cuda_Tecplot.cuh"

void Tecplot_save(){
    int isp,i;
    int sum;
    static int Movie_Init = 1;
    static int PTMovie_Init = 1;
    if(TecplotS_2D_Flag){
        if(cstep == 1){
            Tecplot_2D();
        }else if(cstep >= 1 &&(cstep/TecplotS_2D_Ncycle == 1.0f)){
            Tecplot_2D();
        }
    }
    if(TecplotS_Movie_Flag){
        if((cstep%TecplotS_Movie_Ncycle == 0)){
            TecplotS_Movie_Count++;
            if(TecplotS_Movie_SCYCLE == TecplotS_Movie_Count){
                cudaMemcpy(Host_G_sp, dev_G_sp, nsp * Gsize * sizeof(GPG),cudaMemcpyDeviceToHost);
                cudaMemcpy(vec_G, dev_GvecSet, Gsize * sizeof(GGA),cudaMemcpyDeviceToHost);
                if(Movie_Init){
                    Tecplot_Gsize_Movie(Movie_Init);
                    Movie_Init--;
                }
                Tecplot_Gsize_Movie(Movie_Init);
                TecplotS_Movie_Count = 0;
            }
        }
    }
    if(TecplotS_PT_Movie_Flag){
        if((cstep%TecplotS_PT_Movie_Ncycle == 0)){
            TecplotS_PT_Movie_Count++;
            if(TecplotS_PT_Movie_SCYCLE == TecplotS_PT_Movie_Count){
                cudaMemcpy(Host_G_sp, dev_G_sp, nsp * Gsize * sizeof(GPG),cudaMemcpyDeviceToHost);
                for(isp=0;isp<nsp;isp++){
                    sum = 0;
                    for(i=0;i<Gsize;i++)
                        sum += Host_G_sp[isp*Gsize + i].PtNumInCell;
                    SP[isp].np = sum;
                }
                cudaMemcpy(Host_sp, dev_sp, Total_maxnp * sizeof(GCP),cudaMemcpyDeviceToHost);
                checkCudaErrors(cudaMemcpy(dev_info_sp, SP, nsp * sizeof(Species), cudaMemcpyHostToDevice));
                Copy_GCPtoHCP(SP, Host_sp, PtD, Host_G_sp);
                if(PTMovie_Init){
                    for(isp=0;isp<nsp;isp++) Tecplot_PT_Movie(PTMovie_Init,isp);
                    PTMovie_Init--;
                }
                for(isp=0;isp<nsp;isp++) Tecplot_PT_Movie(PTMovie_Init,isp);
                TecplotS_PT_Movie_Count = 0;
            }
        }  
    }
}
void Tecplot_PT_Movie(int Init,int isp){
    char filename[512];
    FILE *fp;
    int i,j,k;
    sprintf(filename,"%s_PT%d_Movie.dat",InputFile,isp);
    if(Init){
        PT_Movie_S_count = 0;
        fp = fopen(filename,"w");
        fprintf(fp, "TITLE = \"2D PIC Movie\"\n");
	    fprintf(fp, "VARIABLES = \"X(m)\",\"Y(m)\",\"Vx(m/s)\",\"Vy(m/s)\",\"Vz(m/s)\"\n");
        // GEOMETRY
        fprintf(fp, "GEOMETRY\nF=POINT\nCS=GRID\nX=0.00,Y=0.00,Z=0.00\nC=BLACK\nS=GLOBAL\nL=SOLID\nPL=4\nLT=0.1\n");
        fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\nDRAWORDER=AFTERDATA\nMFC=\"\"\n");
        fprintf(fp, "T=RECTANGLE %g %g\n",xlength,ylength);
        for(i=0;i<CondNUM;i++){
            fprintf(fp, "GEOMETRY\nF=POINT\nCS=GRID\n");
            fprintf(fp, "X=%g,Y=%g,Z=0.00\n",CondX0[i]*dx,CondY0[i]*dy);
            fprintf(fp, "C=BLACK\nS=GLOBAL\nL=SOLID\nPL=4\nLT=0.1\n");
            fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\nDRAWORDER=AFTERDATA\nMFC=\"\"\n");
            fprintf(fp, "T=RECTANGLE %g %g\n",CondX1[i]*dx-CondX0[i]*dx,CondY1[i]*dy-CondY0[i]*dy);
        }
        for(i=0;i<DielNUM;i++){
            fprintf(fp, "GEOMETRY\nF=POINT\nCS=GRID\n");
            fprintf(fp, "X=%g,Y=%g,Z=0.00\n",DielX0[i]*dx,DielY0[i]*dy);
            fprintf(fp, "C=BLACK\nS=GLOBAL\nL=SOLID\nPL=4\nLT=0.1\n");
            fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\nDRAWORDER=AFTERDATA\nMFC=\"\"\n");
            fprintf(fp, "T=RECTANGLE %g %g\n",DielX1[i]*dx-DielX0[i]*dx,DielY1[i]*dy-DielY0[i]*dy);
        }
    }else{
        PT_Movie_S_count++;
        fp = fopen(filename,"a");
    }
    fprintf(fp, "ZONE T=\"ZONE %d\"\n",PT_Movie_S_count);
	fprintf(fp, " STRANDID=0, SOLUTIONTIME=1\n");
	fprintf(fp, " I=%d, J=1, K=1, ZONETYPE=Ordered\n",SP[isp].np);
	fprintf(fp, " DATAPACKING=POINT\n");
	fprintf(fp, " DT=(SINGLE SINGLE SINGLE SINGLE SINGLE)\n");
    for (k = 0; k < SP[isp].np; k++) {
		fprintf(fp,"%3.5g %3.5g %3.5g %3.5g %3.5g\n",dx*PtD[isp].x[k],dy*PtD[isp].y[k],PtD[isp].vx[k],PtD[isp].vy[k],PtD[isp].vz[k]);
	}
	fclose(fp);				
}
void Tecplot_Gsize_Movie(int Init){
    char filename[512];
    FILE *fp;
    int isp,i,j,k;
    int Flag1 = 1;
    int Flag2 = 1;
    int Flag3 = 1;
    sprintf(filename,"%s_Movie.dat",InputFile);
    if(Init){
        fp = fopen(filename,"w");
        fprintf(fp, "TITLE = \"2D PIC Movie\"\n");
	    fprintf(fp, "VARIABLES = \"X (m)\", \"Y (m)\", ");
        if(Flag1) fprintf(fp, "\"Ex\", ");
        if(Flag2) fprintf(fp, "\"Ey\", ");
        if(Flag3){
            for(isp=0;isp<nsp;isp++){
                fprintf(fp, "\"Den%d\", ",isp);
            }
        } 
        fprintf(fp, "\n");
        // GEOMETRY
        fprintf(fp, "GEOMETRY\nF=POINT\nCS=GRID\nX=0.00,Y=0.00,Z=0.00\nC=BLACK\nS=GLOBAL\nL=SOLID\nPL=4\nLT=0.1\n");
        fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\nDRAWORDER=AFTERDATA\nMFC=\"\"\n");
        fprintf(fp, "T=RECTANGLE %g %g\n",xlength,ylength);
        for(i=0;i<CondNUM;i++){
            fprintf(fp, "GEOMETRY\nF=POINT\nCS=GRID\n");
            fprintf(fp, "X=%g,Y=%g,Z=0.00\n",CondX0[i]*dx,CondY0[i]*dy);
            fprintf(fp, "C=BLACK\nS=GLOBAL\nL=SOLID\nPL=4\nLT=0.1\n");
            fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\nDRAWORDER=AFTERDATA\nMFC=\"\"\n");
            fprintf(fp, "T=RECTANGLE %g %g\n",CondX1[i]*dx-CondX0[i]*dx,CondY1[i]*dy-CondY0[i]*dy);
        }
        for(i=0;i<DielNUM;i++){
            fprintf(fp, "GEOMETRY\nF=POINT\nCS=GRID\n");
            fprintf(fp, "X=%g,Y=%g,Z=0.00\n",DielX0[i]*dx,DielY0[i]*dy);
            fprintf(fp, "C=BLACK\nS=GLOBAL\nL=SOLID\nPL=4\nLT=0.1\n");
            fprintf(fp, "CLIPPING=CLIPTOVIEWPORT\nDRAWORDER=AFTERDATA\nMFC=\"\"\n");
            fprintf(fp, "T=RECTANGLE %g %g\n",DielX1[i]*dx-DielX0[i]*dx,DielY1[i]*dy-DielY0[i]*dy);
        }
    }else{
        fp = fopen(filename,"a");
    }
    fprintf(fp, "ZONE I = %d, J = %d\n", ngy, ngx);
	fprintf(fp, "ZONETYPE = Ordered, DATAPACKING = POINT\n");
		for(i=0;i<ngx;i++) {
			for(j=0;j<ngy;j++) {
                fprintf(fp,"%g %g ",i*dx,j*dy);
                if(Flag1) fprintf(fp,"%g ",vec_G[i*ngy+j].Ex);
                if(Flag2) fprintf(fp,"%g ",vec_G[i*ngy+j].Ey);
                if(Flag3){
                    for(isp=0;isp<nsp;isp++){
                        fprintf(fp,"%g ",Host_G_sp[isp*Gsize + i*ngy+j].den);
                    }
                } 
				fprintf(fp,"\n");
			}
		}
	fclose(fp);				
}
void Tecplot_2D(){

}