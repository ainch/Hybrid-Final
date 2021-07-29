/******************************************/
/* def.h */
#ifndef __PhysicalConstant__
#define __PhysicalConstant__
#define M_PI 3.14159265358979323846
#define EPS0 8.8542e-12      /* (F/m)  */
#define KB 1.38064852e-23
#define C0 3.0e8
#define SNG_MIN 1e-200
#define DBL_MIN 1E-200
#define NperTORR   8.3221e20
#define AMU		1.66053e-27
#define NVTS    3.3
#endif
#ifndef __Simulation__
#define __Simulation__
enum BoundaryCondition {
    DIRICHLET = 1,
    NEUMANN,
    CONDUCTOR,
    DIELECTRIC,
};
enum BoundaryType {
    LEFT,
    RIGHT,
    UP,
    DOWN,
    UL_CORN,
    UR_CORN,
    LL_CORN,
    LR_CORN,
    NO_FACE,
};
enum LoadType {
    SMARTLOAD,
    UNIFORM,
    EXPONETIAL,
    COSINE,
};
enum GasType {
    ARGON,
    OXYGEN,
    ARO2,
};
#endif
#ifndef __MATH__
#define __MATH__
#define TRUE 1
#define FALSE 0
#ifndef max
#define max(x,y) (((x)>(y)) ? (x) : (y))
#endif
#ifndef min
#define min(x,y) (((x)<(y)) ? (x) : (y))
#endif
#endif

#ifndef TotalDomain
typedef struct TEST
{
    float a, b;
    float *c;
    float *d;
} point;
typedef struct _ChargedParticle
{
	//SIZE = MAXNP
	int CellID;
    float x;  
    float y;
    float vx;
	float vy;
	float vz;
} CP;
#endif

#ifndef TotalDomain
typedef struct tag_species{
	//Constant
	int  Loadtype;		// density load type
	float x_center,x_fall;		// density load position
	float y_center,y_fall;		// density load position
	int  S_ID;			// number name
	float InitDens;		//initial density
	float Temp;			//Load temperature
	float np2c;			// np2c
	int MAXNP;			// maximum np
	//char name[10];		// string name
	//int load_type;		// particle load type
	//float init_xleft, init_xright, init_ybottom, init_ytop;
	//
	//
	//
	//int   PTN;			// if loadtype = 0 or 1, it is not working
	// if loadtype = 2 or 3, minimum number of particle in cell
	// if loadtype = 4, Number of particle for load type = 4 case	
	//float q,m;			// mass and Charge q, charge q is q*np2c
	//float vti;			// Thermal velocity
	//int np;				// number of particle
	//int *MaxPtNumInCell;// Maximum Particle number in cell
	//int *PtNumInCell;	// Particle number in cell
	//float *vec_den,**den; // density
	//float qm;			// charge q times mass m
	//float ascale;		// charge q times mass m times time step dt
	//float Escale;		// half of mass m of constant Q
	//float **sigma;		// for surface charge at dielectric wall
	//int 	*CSeecFlag;
	//float 	*CSeec;
	//int 	*DSeecFlag;
	//float 	*DSeec;
//	float q_density;
//	int AvePtNumInCell;
	//int view_np;
	// diag parameter
	//float *vec_Tx, *vec_Ty, *vec_Tz;
	//float *vec_acc_den,**ave_den;
	//float *vec_acc_Tx, *vec_acc_Ty, *vec_acc_Tz;
	//float *vec_acc_fx, *vec_acc_fy, *vec_acc_fz;
	//float *vec_JdotE, *vec_Jx, *vec_Jy;
	//float *vec_acc_JdotE, *vec_acc_Jx, *vec_acc_Jy;
	//int *HistNp;
	//float **Current_hist;
	//float sum_eloss, ave_eloss;
	//float sum_wall_eloss, ave_wall_eloss;
	//float *KE_hist;
}Species;
#endif
#ifndef TotalDomain
typedef struct fluid{
	//int  Cname;					// species name
	//char name[10];				// string name
	int   Loadtype;				// density load type
	float x_center,x_fall;		// density load position
	float y_center,y_fall;		// density load position
	float InitDens;				// initial density
	float Temp;					// Temperature
	//float m;					// Mass
	//float *HistDen;				//
	//float **den;				// density of each of cell using update_source();
	//float **sum_den;			// for Average density
	//float **ave_den;			// AverageDensity
	//float **D;					// Diffusion coefficient each of cell
	//float **source;				// source each of cell
	//float **gummel_ax;			// Velocity component for x direction flux, positive direction is left->right
	//float **gummel_bx;			// Velocity component for x direction flux, positive direction is left->right
	//float **gummel_ay;			// Velocity component for y direction flux, positive direction is down->up
	//float **gummel_by;			// Velocity component for y direction flux, positive direction is down->up
	//float **flux_x;				// x direction flux
	//float **flux_y;				// y direction flux
	//float **vt;					// mean velocity
	//float GasThermalVelocity;	// thermal velosity
	//int C_E_flag;				// Continuity equation for steady state
	//int C_E_Check;				// Cycle to judge convergence
	//int C_E_num;				// Cycle to solve the continuity equation
	//float C_E_conver;           // convergence judge value
	//float **T;					// Temperature
	//float n0;					// unit change density [cm^-3]
	//float sig;					// sigma
} Fluid;
#endif
#ifndef TotalDomain
typedef struct BackGround{
	//int  Cname;					// species Number
	//char name[10];				// string name
	float Pres;				    // Pressure
	float Temp;				// Temperature
	//float GasDensity;			// density [m^-3]
	//float GasThermalVelocity;	// thermal velocity
	//float m;					// Mass
	//float **den;				// Temperature gadient density
//	float vt;
//	float sig;
} BackG;
#endif




#ifndef TotalDomain
typedef struct tag_totaldomain{
	//int	ncx,ncy,ngx,ngy;
	//float *phi;
	//int *boundary_flag;
}TotalDomain;
#endif

#ifndef TotalDomain
typedef struct tag_dev_species{
	// Constant
	//int *np;
	//int *IndexArray;
	//float *den,*sigma;
	//float *real_den;
	//float *sum_vx,*sum_vy,*sum_vz;
	//float *Tx, *Ty, *Tz;
	//int   *CellID;
	//float *Lx,*Ly;					// Localized x,y
	//float *vx, *vy, *vz;
	//int   *MaxPtNumInCell;
	//int   *PtNumInCell;
	//int   *PtNumMoveInterCell;
	//int   *CSeecFlag;
	//float *CSeec;
	//int   *DSeecFlag;
	//float *DSeec;
	//float *acc_SEED;									// Second Electron Emission Diagnostic
	//float *acc_CWB;										// Collision with wall
	//float *acc_den;
	//float *acc_Tx, *acc_Ty, *acc_Tz;
	//float *acc_fx, *acc_fy, *acc_fz;
	//float *JdotE, *Jx, *Jy;
	//float *acc_JdotE, *acc_Jx, *acc_Jy;
	//float *charge;										// 0(zero charge) or 1 (charge)
	//float *eloss;
	//float *wall_eloss;
	//float *KE;
}dev_Species;
#endif
#ifndef TotalDomain
typedef struct MCC_Diagnostic{
	//float *vec_rate;
	//float **rate;
	//float ave_eloss;
	//float eloss;
}MCCR;
#endif
#ifndef TotalDomain
typedef struct dev_MCC_Diagnostic{
	//float *num;	  // cell
	//float *rate;  // grid
	//float *eloss; // grid
}dev_MCCR;
#endif
#ifndef TotalDomain
typedef struct Function{
	//int size;
	//float engy;
	//float *x, *y;
} Func;
#endif
#ifndef TotalDomain
typedef struct plasma_region {
	//int x1, x2, y1, y2;		//
	//float fg1, fg2;			// boundary condition for flux
	//int region_num;//
} plasma_region;
#endif
#ifndef TotalDomain
typedef struct tag_EDF{
	//int num, n_bin;
	//float E_min, E_max, dE;
	//float *xl, *xr, *yb, *yt;
	//float **edf, **edf_sum, **edf_ave;
	//float *dev_edf;
	//int *xs, *xe, *ys, *ye;
}EDF;
#endif




