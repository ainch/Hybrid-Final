#include "load.h"
void V001_LoadDUMP(FILE *LF){

}
void V000_LoadDUMP(FILE *LF){
    int buf0,isp,i;
    float buf1;
    // Time
    fread(&t, 8, 1, LF);
    fread(&tstep, 4, 1, LF);
    fread(&cstep, 4, 1, LF);
    // Gas info
    fread(&buf0, 4, 1, LF); // MainGas
    if(buf0 != MainGas){
        printf(" - Error! : Gas species do not match. input[%d] != dump[%d]\n",MainGas,buf0);
        exit(1);
    }
    for(isp=0;isp<nsp;isp++){
        fread(&buf1, 4, 1, LF); // NP2C
        if(buf1 != SP[isp].np2c){
            printf(" - Warning! : %s np2c do not match. input[%g] != dump[%g]\n",buf1,SP[isp].np2c);
        }
        SP[isp].np2c = buf1;
        fread(&SP[isp].np, 4, 1, LF); // NP
        for (i = 0; i < SP[isp].np; i++) {
            fread(&PtD[isp].CellID[i], 4, 1, LF);
            fread(&PtD[isp].x[i], 4, 1, LF);
            fread(&PtD[isp].y[i], 4, 1, LF);
            fread(&PtD[isp].vx[i], 4, 1, LF);
            fread(&PtD[isp].vy[i], 4, 1, LF);
            fread(&PtD[isp].vz[i], 4, 1, LF);
        }
    }
    //Fluid info
    fread(&Conti_Flag, 4, 1, LF);
    for(isp=0;isp<nfsp;isp++){
        for(i=0;i<ncx;i++){
            fread(Fluid_sp[isp].den[i], 4, ncy, LF);
            fread(Fluid_sp[isp].Source[i], 4, ncy, LF);
        }
        for(i=0;i<ngx;i++)
            fread(Fluid_sp[isp].flux_x[i], 4, ncy, LF);
        for(i=0;i<ncx;i++)
            fread(Fluid_sp[isp].flux_y[i], 4, ngy, LF);
    }
    // History data
    fread(&hist_count, 4, 1, LF);
    fread(&dHIST, 4, 1, LF);
	fread(t_array, 4, hist_count, LF);
	for (isp = 0; isp < nsp; isp++) {
		fread(HistPt[isp].np, 4, hist_count, LF);
	}
    for (isp = 0; isp < nfsp; isp++) {
		fread(HistFG[isp].np, 4, hist_count, LF);
	}
	fread(iter_array, 4, hist_count, LF);
    for (i = 0; i < CondNUMR; i++) {
		fread(Current_hist[i], 4, hist_count, LF);
		for (isp = 0; isp < nsp; isp++) {
			fread(SP_current_hist[isp][i], 4, hist_count, LF);
		}
		fread(Volt_hist[i], 4, hist_count, LF);
		fread(Volt_cond_hist[i], 4, hist_count, LF);
		fread(Surf_charge_hist[i], 4, hist_count, LF);
	}
    fread(&hist_ave_count, 4, 1, LF);
    fread(t_ave_array, 4, hist_ave_count, LF);
    for (isp = 0; isp < nsp; isp++) {
		fread(Hist_ave_Pt[isp].np, 4, hist_ave_count, LF);
        fread(Hist_ave_Pt_stack[isp].np, 4, hist_ave_count, LF);
	}
    // 2D data
    for (i = 0; i < nsp*Gsize; i++) {
        fread(&Host_G_sp[i].den, 4, 1, LF);    // pt density : GPG
        fread(&Host_G_sp[i].ave_den, 4, 1, LF); // pt Ave_density : GPG
        fread(&Host_G_sp[i].sigma, 4, 1, LF); // pt sigma
    }
    for (i = 0; i < Gsize; i++) {
        if(DumpFlag){
            fread(&vec_G[i].Temp, 4, 1, LF); // BG Tem : GGA
            fread(&vec_G[i].BackVel1, 4, 1, LF); // BG vel : GGA
            fread(&vec_G[i].BackDen1, 4, 1, LF); // BG density : GGA
            fread(&vec_G[i].BackVel2, 4, 1, LF); // BG vel : GGA
            fread(&vec_G[i].BackDen2, 4, 1, LF); // BG density : GGA
        }else{
            fread(&buf1, 4, 1, LF); // BG Tem : GGA
            fread(&buf1, 4, 1, LF); // BG vel : GGA
            fread(&buf1, 4, 1, LF); // BG density : GGA
            fread(&buf1, 4, 1, LF); // BG vel : GGA
            fread(&buf1, 4, 1, LF); // BG density : GGA
        }
        fread(&vec_G[i].Lap_Pot, 4, 1, LF); // BG Lap_Pot
        fread(&vec_G[i].Pois_Pot, 4, 1, LF); // BG Pois_Pot
        fread(&vec_G[i].Ex, 4, 1, LF); // BG Ex
        fread(&vec_G[i].Ey, 4, 1, LF); // BG Ey
    } 
    // Field
    fread(vec_Potential, 4, Gsize, LF);
    fread(ave_Potential, 4, Gsize, LF);
    fread(ave_Source, 4, Gsize, LF);
    fread(ave_Sigma, 4, Gsize,LF);
    fread(ave_Ex, 4, Gsize, LF);
    fread(ave_Ey, 4, Gsize,LF);
    // MCC
    fread(ave_MCC_rate, 4, TnRct * Gsize, LF); // MCC AVE RATE : mccrate
}
void LoadDumpFile(){
    FILE *DumpDeck;
    int KEY2, KEY1, KEY0;
    printf("Read The Dump file.\n");
    DumpDeck = fopen(DumpFile, "r+b");
    fread(&KEY2, 4, 1, DumpDeck);
    fread(&KEY1, 4, 1, DumpDeck);
    fread(&KEY0, 4, 1, DumpDeck);
    printf(" - Dump Version : [%d].[%d].[%d]\n",KEY2, KEY1, KEY0);
    if(KEY2>=0 && KEY1>=0 && KEY0>=0)  V000_LoadDUMP(DumpDeck);
    if(KEY2>=0 && KEY1>=0 && KEY0>=1)  V001_LoadDUMP(DumpDeck);
    fclose(DumpDeck);
    printf("Dump file read complete.\n");
}
void CG_Matrix_Setting(float *A, int *Ai, int *Aj, float **b, float *M, float *Atemp, float *btemp) {
	int i, j, k, x1, y1, x2, y2;
	int m_count, l_count;
	int LEFT_NBC, RIGHT_NBC, UP_NBC, DOWN_NBC;
	float LEFT_COEF, RIGHT_COEF, UP_COEF, DOWN_COEF;
	float LEFT_COEF1, RIGHT_COEF1, UP_COEF1, DOWN_COEF1;
	float dydx, dxdy;
	int bd_num, left_p, right_p;
	float temp;

	dydx=0.5*dy/dx;
	dxdy=0.5*dx/dy;
	m_count=0; //element number
	l_count=0; //matrix line number

	Ai[0]=1; //one-based indexing
	for(i=0;i<ngx;i++) {
		for(j=0;j<ngy;j++) {
			if(A_idx[i][j]) {
				//Calculate coefficients
				x1=(i==0) ? i : i-1;
				y1=(j==ncy) ? j-1 : j;
				y2=(j==0) ? j : j-1;
				LEFT_COEF=dydx*(vec_C[x1*ncy+y1].eps_r+vec_C[x1*ncy+y2].eps_r);

				x1=(i==ncx) ? i-1 : i;
				y1=(j==ncy) ? j-1 : j;
				y2=(j==0) ? j : j-1;
				RIGHT_COEF=dydx*(vec_C[x1*ncy+y1].eps_r+vec_C[x1*ncy+y2].eps_r);
			
				y1=(j==ncy) ? j-1 : j;
				x1=(i==ncx) ? i-1 : i;
				x2=(i==0) ? i : i-1;
				UP_COEF=dxdy*(vec_C[x1*ncy+y1].eps_r+vec_C[x2*ncy+y1].eps_r);
				
				y1=(j==0) ? j : j-1;
				x1=(i==ncx) ? i-1 : i;
				x2=(i==0) ? i : i-1;
				DOWN_COEF=dxdy*(vec_C[x1*ncy+y1].eps_r+vec_C[x2*ncy+y1].eps_r);
			 	
				LEFT_NBC=(i==0) ? 2 : 1;
				RIGHT_NBC=(i==ncx) ? 2 : 1;
				UP_NBC=(j==ncy) ? 2 : 1;
				DOWN_NBC=(j==0) ? 2 : 1;

				LEFT_COEF1 = LEFT_COEF;
				RIGHT_COEF1 = RIGHT_COEF;
				UP_COEF1 = UP_COEF;
				DOWN_COEF1 = DOWN_COEF;
                // LEFT SUFACE CHECK 
                if(i) {
                    if(vec_G[(i-1)*ngy+(j)].CondID) { // If the left grid is the Conductor
                        b[vec_G[(i-1)*ngy+(j)].CondID-1][l_count]+=LEFT_COEF;
                        btemp[l_count]+=LEFT_COEF*vec_G[(i-1)*ngy+(j)].Temp;
                    }
                    else if(vec_G[(i-1)*ngy+(j)].Boundary==DIRICHLET) {
                        btemp[l_count]+=LEFT_COEF*vec_G[(i-1)*ngy+(j)].Temp;
                    }
                    else {
                        A[m_count]=(-LEFT_COEF)*RIGHT_NBC;
                        Aj[m_count]=A_idx[i-1][j];
                        Atemp[m_count]=(-LEFT_COEF)*RIGHT_NBC;
                        m_count++;
                    }
                }
                //DOWN
                if(j) { 
                    if(vec_G[(i)*ngy+(j-1)].CondID) {
                        b[vec_G[(i)*ngy+(j-1)].CondID-1][l_count]+=DOWN_COEF;
                        btemp[l_count]+=DOWN_COEF*vec_G[(i)*ngy+(j-1)].Temp;
                    }
                    else if(vec_G[(i)*ngy+(j-1)].Boundary==DIRICHLET) {
                        btemp[l_count]+=DOWN_COEF*vec_G[(i)*ngy+(j-1)].Temp;
                    }
                    else {
                        A[m_count]=(-DOWN_COEF)*UP_NBC;
                        Aj[m_count]=A_idx[i][j-1];
                        Atemp[m_count]=(-DOWN_COEF);
                        m_count++;
                    }
                }
                //MIDDLE
                A[m_count]=LEFT_COEF1+DOWN_COEF1+UP_COEF1+RIGHT_COEF1;
                M[l_count]=1/A[m_count];
                Aj[m_count]=A_idx[i][j];
                Atemp[m_count]=LEFT_COEF+DOWN_COEF+UP_COEF+RIGHT_COEF;
                m_count++;
                //UP
                if(j!=ncy) {
                    if(vec_G[(i)*ngy+(j+1)].CondID) {
                        b[vec_G[(i)*ngy+(j+1)].CondID-1][l_count]+=UP_COEF;
                        btemp[l_count]+=UP_COEF*vec_G[(i)*ngy+(j+1)].Temp;
                    }
                    else if(vec_G[(i)*ngy+(j+1)].Boundary==DIRICHLET) {
                        btemp[l_count]+=UP_COEF*vec_G[(i)*ngy+(j+1)].Temp;
                    }
                    else{
                        A[m_count]=(-UP_COEF)*DOWN_NBC;
                        Aj[m_count]=A_idx[i][j+1];

                        Atemp[m_count]=(-UP_COEF);
                        m_count++;
                    }
                }
                //RIGHT
                if(i!=ncx) { //for Neumann B.C
                    if(vec_G[(i+1)*ngy+(j)].CondID) {// DIRICHLET B.C
                        b[vec_G[(i+1)*ngy+(j)].CondID-1][l_count]+=RIGHT_COEF;
                        btemp[l_count]+=RIGHT_COEF*vec_G[(i+1)*ngy+(j)].Temp;
                    }
                    else if(vec_G[(i+1)*ngy+(j)].Boundary==DIRICHLET) {
                        btemp[l_count]+=RIGHT_COEF*vec_G[(i+1)*ngy+(j)].Temp;
                    }
                    else {
                        A[m_count]=(-RIGHT_COEF)*LEFT_NBC;
                        Aj[m_count]=A_idx[i+1][j];
                        Atemp[m_count]=(-RIGHT_COEF);
                        m_count++;
                    }
                }
                l_count++;
                Ai[l_count]=m_count+1;
			}
		}
	}
}