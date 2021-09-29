#include "load.h"

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