#include "Particle.h"

void SetParticleLoad(int isp, float Ninit, int load_type, float x_left,float x_right, float y_top, float y_bottom, float vti) {
    int i, j, k, n;
	float y0, y1, ylen, x0, x1, xlen, inv_den;
	int n_per_cell, Num, index;
	float xx, yy;
	float **wv; 			//weighting value

	n_per_cell = (int) (Ninit * dx * dy * zlength / SP[isp].np2c);
	SP[isp].np = 0;
    y0 = y_top / dy;
	y1 = y_bottom / dy - 1;
	ylen = y1 - y0 + 1;
	x0 = x_left / dx;
	x1 = x_right / dx - 1;
	xlen = x1 - x0 + 1;
	///// Loading the velocities and positions depending on\
	///// the given flag.  This Loader can load particles\
	///// using a quiet start method, uniform, or randomly
	switch (load_type) {
	case UNIFORM: {
		index = 0;
		for (i = 0; i < ncx; i++) {
			for (j = 0; j < ncy; j++) {
				if (vec_C[i*ncy+j].PlasmaRegion == 1) {
					if (i >= x0 && i < x1 && j >= y0 && j < y1){
						Num = n_per_cell;
                    }else{
                        Num = 0;
                    }	
					for (n = 0; n < Num; n++) {
						PtD[isp].vx[index] = PtD[isp].vy[index] = PtD[isp].vz[index] = 0.0;
						maxwellv(&PtD[isp].vx[index], &PtD[isp].vy[index],&PtD[isp].vz[index], vti);
						PtD[isp].x[index] = (float) i + frand();
						PtD[isp].y[index] = (float) j + frand();
						PtD[isp].CellID[index] = i * ngy + j;
						index++;
					}
				}
			}
		}
		SP[isp].np = index;
		printf("\tTotal Number of %s particle = %d, N_per_cell : %d\n",SP[isp].name, SP[isp].np, n_per_cell);
		if (SP[isp].np > NP_LIMIT) {
			printf("LOAD: too many particles, species %s\n", SP[isp].name);
			exit(1);
		}
		break;
	}
	case EXPONETIAL: {
		wv = MFMalloc(ncx,ncy);
		MFInit(wv,0.0,ncx,ncy);
		for (i = 0; i < ncx; i++) {
			for (j = 0; j < ncy; j++) {
				xx = (float) i * dx + dx / 2;
				yy = (float) j * dy + dy / 2;
				wv[i][j] = exp(
						-1 * ((xx - x_left) / x_right)
								* ((xx - x_left) / x_right))
						* exp(
								-1 * ((yy - y_bottom) / y_top)
										* ((yy - y_bottom) / y_top));
			}
		}
		index = 0;
		for (i = 0; i < ncx; i++) {
			for (j = 0; j < ncy; j++) {
				if (vec_C[i*ncy+j].PlasmaRegion == 1) {
					Num = n_per_cell * wv[i][j];
					if (Num < 1) {
						Num = SP[isp].Ratio;
					}
					for (n = 0; n < Num; n++) {
						PtD[isp].vx[index] = PtD[isp].vy[index] = PtD[isp].vz[index] = 0.0;
						maxwellv(&PtD[isp].vx[index], &PtD[isp].vy[index],&PtD[isp].vz[index], vti);
						PtD[isp].x[index] = (float) i + frand();
						PtD[isp].y[index] = (float) j + frand();
						PtD[isp].CellID[index] = i * ngy + j;
						index++;
					}
				}
			}
		}
		SP[isp].np = index;
		printf("\tTotal Number of %s particle = %d\n",SP[isp].name, SP[isp].np);
		if (SP[isp].np > NP_LIMIT) {
			printf("LOAD: too many particles, species %s\n", SP[isp].name);
			exit(1);
		}
		break;
	}
	case COSINE: {
		wv = MFMalloc(ncx,ncy);
		MFInit(wv,0.0,ncx,ncy);
		for (i = 0; i < ncx; i++) {
			for (j = 0; j < ncy; j++) {
				xx = (float) i * dx + dx / 2;
				yy = (float) j * dy + dy / 2;
				if(xx>=0.06){
					wv[i][j] = 0;
				}else{
					wv[i][j] = cos((xx - x_left) * M_PI / 2 / x_right)
						* cos((yy - y_bottom) * M_PI / 2 / y_top);
				}
				
			}
		}
		index = 0;
		for (i = 0; i < ncx; i++) {
			for (j = 0; j < ncy; j++) {
				if (vec_C[i*ncy+j].PlasmaRegion == 1) {
					Num = (int) n_per_cell * wv[i][j];
					if (Num < 1) {
						Num = SP[isp].Ratio;
					}
					for (n = 0; n < Num; n++) {
						PtD[isp].vx[index] = PtD[isp].vy[index] = PtD[isp].vz[index] = 0.0;
						maxwellv(&PtD[isp].vx[index], &PtD[isp].vy[index],&PtD[isp].vz[index], vti);
						PtD[isp].x[index] = (float) i + frand();
						PtD[isp].y[index] = (float) j + frand();
						PtD[isp].CellID[index] = i * ngy + j;
						index++;
					}
				}
			}
		}
		SP[isp].np = index;
		printf("\tTotal Number of %s particle = %d\n",SP[isp].name, SP[isp].np);
		if (SP[isp].np > NP_LIMIT) {
			printf("LOAD: too many particles, species %s\n", SP[isp].name);
			exit(1);
		}
		break;
	}
	case NP_RAIO: {
		index = 0;
		for (i = 0; i < ncx; i++) {
			for (j = 0; j < ncy; j++) {
				if (vec_C[i*ncy+j].PlasmaRegion == 1) {
					if (i >= x0 && i < x1 && j >= y0 && j < y1)
						Num = SP[isp].Ratio;
					else
						Num = 0;
					for (n = 0; n < Num; n++) {
						PtD[isp].vx[index] = PtD[isp].vy[index] = PtD[isp].vz[index] = 0.0;
						maxwellv(&PtD[isp].vx[index], &PtD[isp].vy[index],&PtD[isp].vz[index], vti);
						PtD[isp].x[index] = (float) i + frand();
						PtD[isp].y[index] = (float) j + frand();
						PtD[isp].CellID[index] = i * ngy + j;
						index++;
					}
				}
			}
		}
		SP[isp].np = index;
		printf("\tTotal Number of %s particle = %d\n",SP[isp].name, SP[isp].np);
		if (SP[isp].np > NP_LIMIT) {
			printf("LOAD: too many particles, species %s\n", SP[isp].name);
			exit(1);
		}
		break;
	}
    case SMARTLOAD: {
        printf("LoadType : SMARTLOAD DEVELOP YET");
        exit(1);
        break;
    }
	default: {
		puts("LOAD: Bad value for loader flag");
		exit(-1);
		break;
	}
	}
}
void maxwellv(float *vx_local, float *vy_local, float *vz_local, float vti)
{
  static int nvel, init_flag= 1;

  int i, n;
  float vmag, aphi, dv, rr, sintheta, costheta;

  if (init_flag) {
    nvel= (int)(1./(1-F((float)NVTS)));
    dv = sqrt(M_PI)/(4.0*nvel);
    if(nvel > 100000)
      puts("Warning: Your choice of NVTS has made nvel > 1e5");

    vsave= (float *) malloc(nvel*sizeof(float));

    init_flag= 0;
    i=n=0;
    for (n=0; n<nvel; n++) {
      rr=(1.0*n)/nvel;
      while (F(i*dv)< rr) i++;
      vmag=i*dv;
			vsave[n]=sqrt(2.0)*i*dv;
    }
	 h_nvel=nvel;
//	 fprintf(stderr,"h_vel = %d \n", h_nvel);
  }

  n = (int)((nvel-1)*frand());
  aphi=2*M_PI*frand();
  costheta = 1-2*frand();
  sintheta = sqrt(1-costheta*costheta);
  *vx_local = vti*vsave[n]*sintheta*cos(aphi);
  *vy_local = vti*vsave[n]*sintheta*sin(aphi);
  *vz_local = vti*vsave[n]*costheta;
}
float frand()
{
  long a = 16807, m = 2147483647, q = 127773, r = 2836;
  long hi, lo;
  float fnumb;
  /* static long seed=31207321; */

  hi = seed/q;
  lo = seed - q*hi;
  seed = a*lo - r*hi;
  /* "seed" will always be a legal integer of 32 bits (including sign). */
  if(seed <= 0) seed = seed + m;
  fnumb = seed/2147483646.0;

  return(fnumb);
}
float F(float v)
{
  return(-2*v*exp(-v*v)/sqrt(M_PI) +erf(v));
}
