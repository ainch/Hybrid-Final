/******************************************/
/* def.h */
#ifndef __PhysicalConstant__
#define __PhysicalConstant__
#define M_PI 3.14159265358979323846
#define EPS0 8.8542e-12      /* (F/m)  */
#define KB 1.38064852e-23
#define C0 3.0e8
#define CQ 1.602e-19
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
    UNIFORM,
    EXPONETIAL,
    COSINE,
	NP_RAIO,
	SMARTLOAD,
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
typedef struct __Host_Gsize_Array{
	int Boundary;  //Boundary Condition Constant 0~4
	int CondID;		// Conductor ID , NO Conductor is zero
	float Temp;
	float BackDens;
    int face;
    float area;
}HGA;
typedef struct __Host_Csize_Array{
    int PlasmaRegion;
    float eps_r;
}HCA;
typedef struct __Device_Gsize_Array{

}dev_GA;
typedef struct __Device_Csize_Array{

}dev_CA;
typedef struct _Host_Charged_Particle
{ // Bufer in CPU
	float *den;
	int *MaxPtNumInCell;// SIZE(Gsize)
	int *PtNumInCell; 	// SIZE(Gsize)
	int *CellID;// SIZE(NP_LIM)
    float *x;   // SIZE(NP_LIM)
    float *y;  	// SIZE(NP_LIM)
    float *vx;  // SIZE(NP_LIM)
	float *vy;  // SIZE(NP_LIM)
	float *vz;  // SIZE(NP_LIM)
}HCP;
typedef struct _Device_Charged_Particle
{
	//SIZE = MAXNP
	int CellID;
    float x;  
    float y;
    float vx;
	float vy;
	float vz;
}dev_CP;
typedef struct tag_species
{
	//Constant
	char name[10];
	int  Loadtype;				// density load type
	float x_center,x_fall;		// density load position
	float y_center,y_fall;		// density load position
	int  S_ID;					// species ID
	float InitDens;				//initial density
	float Temp;					//Load temperature
	float np2c;					// np2c
	int Ratio;
	int MAXNP;					// maximum np
	float mass;
	float q;
	float q_density;
	float vti;
	int np;
	float qm;
	float Escale;
	float Ascale;
}Species;
typedef struct fluid{
	char name[10];				// string name
	int   Loadtype;				// density load type
	float x_center,x_fall;		// density load position
	float y_center,y_fall;		// density load position
	float InitDens;				// initial density
	float Temp;					// Temperature
	float mass;
} Fluid;
typedef struct BackGround{
	char name[10];				// string name
	float Ratio;
	float Pres;				    // Pressure
	float Temp;				// Temperature
	float mass;
} BackG;
typedef struct _Collision_Flag
{
    float Flag; // CX = Flag * CX
	float mofM; // target/projectile
	float Th_e; //Threshold energy
	float RR; // Reaction RATE
} CollF;
typedef struct _Argon_Collision_Data
{
    float xe;  // log10(xee)
	float xee; // pow(10,xe)
	float cx_0; //"0.e+Ar>e+Ar"
	float cx_1; //"1.e+Ar>e+Ar*"
	float cx_2; //"2.e+Ar>e+Ar*m"
	float cx_3; //"3.e+Ar>2e+Ar^"
	float cx_4; //"4.e+Ar*m>2e+Ar^"
	float cx_5; //"5.Ar+Ar^>Ar^+Ar"
	float cx_6; //"6.Ar+Ar^>Ar+Ar^"
} ArCollD;
typedef struct _Oxygen_Collision_Data
{
    float xe;  // log10(xee)
	float xee; // pow(10,xe)
	float cx_0; //"0.e+O2>e+O2");
    float cx_1; //"1.e+O2>e+O2*");
    float cx_2; //"2.e+O2>e+O2*");
    float cx_3; //"3.e+O2>e+O2A");
    float cx_4; //"4.e+O2>e+O2B");
    float cx_5; //"5.e+O2>e+O2*");
    float cx_6; //"6.e+O2>OP+O-");
    float cx_7; //"7.e+O2>e+2OP");
    float cx_8; //"8.e+O2>e+OP+OD");
    float cx_9; //"9.e+O2>e+2OD");
    float cx_10; //"10.e+O2>2e+O2^");
    float cx_11; //"11.e+O2>e+OP+O*");
    float cx_12; //"12.e+O2>e+O^+O-");
    float cx_13; //"13.e+O2>2e+O^+OP");
    float cx_14; //"14.e+O2A>2e+O2+");
    float cx_15; //"15.e+O2A>OP+O-");
    float cx_16; //"16.e+O2A>e+O2");
    float cx_17; //"17.e+O2A>e+O2");
    float cx_18; //"18.e+O2A>e+2OP");
    float cx_19; //"19.e+O2A>e+OP+OD");
    float cx_20; //"20.e+O2A>e+2OD");
    float cx_21; //"21.e+O2A>2e+O^+OP");
    float cx_22; //"22.e+O2B>2e+O2^");
    float cx_23; //"23.e+O2B>OP+O-");
    float cx_24; //"24.e+O2B>e+O2");
    float cx_25; //"25.e+O2B>e+O2");
    float cx_26; //"26.e+O2B>e+2O");
    float cx_27; //"27.e+O2B>e+OP+OD");
    float cx_28; //"28.e+O2B>e+2OD");
    float cx_29; //"29.e+O2B>2e+O^+OP");
    float cx_30; //"30.e+O->2e+OP");
    float cx_31; //"31.e+O2^>OP+OD");
    float cx_32; //"32.e+OP>e+OP");
    float cx_33; //"33.e+OP>e+OD");
    float cx_34; //"34.e+OP>e+O*");
    float cx_35; //"35.e+OP>e+O*");
    float cx_36; //"36.e+OP>e+O*");
    float cx_37; //"37.e+OP>2e+O^");
    float cx_38; //"38.e+OP>e+O*");
    float cx_39; //"39.e+OD>2e+O^");
    float cx_40; //"40.e+OD>e+OP");
    float cx_41; //"41.O-+O2>O-+O2");
    float cx_42; //"42.O-+O2>e+OP+O2");
    float cx_43; //"43.O-+OP>e+O2");
    float cx_44; //"44.O-+O2^>OP+O2");
    float cx_45; //"45.O-+O^>2OP");
    float cx_46; //"46.O-+O2A>e+OP+O2");
    float cx_47; //"47.O2^+OP>O2+O^");
    float cx_48; //"48.O2^+O2>O2+O2^");
    float cx_49; //"49.O2^+O2>O2^+O2");
    float cx_50; //"50.O2^+O2>O^+OP+O2");
    float cx_51; //"51.O2^+O2A>O2+O2^");
    float cx_52; //"52.O2^+O2B>O2+O2^");
    float cx_53; //"53.O^+O2>OP+O2^");
    float cx_54; //"54.O^+O2>O^+O2");
    float cx_55; //"55.O^+OP>OP+O^");
    float cx_56; //"56.O^+O2A>O2^+OP");
    float cx_57; //"57.O^+O2B>O2^+OP");
} O2CollD;
typedef struct _ArO2_Collision_Data
{
    float xe;  // log10(xee)
	float xee; // pow(10,xe)
    float cx_0; //0.e+Ar>e+Ar");
    float cx_1; //1.e+Ar>e+Ar*");
    float cx_2; //2.e+Ar>e+Ar*m");
    float cx_3; //3.e+Ar>2e+Ar+");
    float cx_4; //4.e+Ar*m>2e+Ar^");
    float cx_5; //5.e+O2>e+O2");
    float cx_6; //6.e+O2>e+O2*");
    float cx_7; //7.e+O2>e+O2*");
    float cx_8; //8.e+O2>e+O2A");
    float cx_9; //9.e+O2>e+O2B");
    float cx_10; //"10.e+O2>e+O2*");
    float cx_11; //"11.e+O2>OP+O-");
    float cx_12; //"12.e+O2>e+2OP");
    float cx_13; //"13.e+O2>e+OP+OD");
    float cx_14; //"14.e+O2>e+2OD");
    float cx_15; //"15.e+O2>2e+O2+");
    float cx_16; //"16.e+O2>e+OP+O*");
    float cx_17; //"17.e+O2>e+O++O-");
    float cx_18; //"18.e+O2>2e+O^+OP");
    float cx_19; //"19.e+O2A>2e+O2^");
    float cx_20; //"20.e+O2A>OP+O-");
    float cx_21; //"21.e+O2A>e+O2");
    float cx_22; //"22.e+O2A>e+O2");
    float cx_23; //"23.e+O2A>e+2O"); 
    float cx_24; //"24.e+O2A>e+OP+OD");
    float cx_25; //"25.e+O2A>e+2OD");
    float cx_26; //"26.e+O2A>2e+O^+OP");
    float cx_27; //"27.e+O2B>2e+O2^");
    float cx_28; //"28.e+O2B>OP+O-");
    float cx_29; //"29.e+O2B>e+O2");
    float cx_30; //"30.e+O2B>e+O2");
    float cx_31; //"31.e+O2B>e+2O");
    float cx_32; //"32.e+O2B>e+OP+OD");
    float cx_33; //"33.e+O2B>e+2OD");
    float cx_34; //"34.e+O2B>2e+O++OP");
    float cx_35; //"35.e+O->2e+OP");
    float cx_36; //"36.e+O2+>OP+OD ");
    float cx_37; //"37.e+OP>e+OP");
    float cx_38; //"38.e+OP>e+OD");
    float cx_39; //"39.e+OP>e+O*");
    float cx_40; //"40.e+OP>e+O*");
    float cx_41; //"41.e+OP>e+O*");
    float cx_42; //"42.e+OP>2e+O^");
    float cx_43; //"43.e+OP>e+O*");
    float cx_44; //"44.e+OD>2e+O+");
    float cx_45; //"45.e+OD>e+O");
    float cx_46; //"46.O-+O2>O-+O2");
    float cx_47; //"47.O-+O2>e+OP+O2");
    float cx_48; //"48.O-+OP>e+O2");
    float cx_49; //"49.O-+O2^>OP+O2");
    float cx_50; //"50.O-+O^>2OP");
    float cx_51; //"51.O-+O2A>e+OP+O2");  
    float cx_52; //"52.O2^+OP>O2+O^");
    float cx_53; //"53.O2^+O2>O2+O2^");
    float cx_54; //"54.O2^+O2>O2^+O2");   
    float cx_55; //"55.O2^+O2>O^+OP+O2");
    float cx_56; //"56.O2^+O2A>O2+O2^");
    float cx_57; //"57.O2^+O2B>O2+O2^");
    float cx_58; //"58.O2^+Ar>O2+Ar^");
    float cx_59; //"59.O2^+Ar>O2^+Ar^");
    float cx_60; //"60.O^+O2>OP+O2^");
    float cx_61; //"61.O^+O2>O^+O2");
    float cx_62; //"62.O^+OP>OP+O^");
    float cx_63; //"63.O^+O2A>O2^+OP");
    float cx_64; //"64.O^+O2B>O2^+OP");
    float cx_65; //"65.Ar^+Ar>Ar+Ar^");
    float cx_66; //"66.Ar^+Ar>Ar++Ar");
    float cx_67; //"67.Ar^+O2>O2+Ar^");
} ArO2CollD;
#endif
/*
#ifndef TotalDomain
typedef struct tag_species{
	//Constant
	char name[10];
	int  Loadtype;		// density load type
	float x_center,x_fall;		// density load position
	float y_center,y_fall;		// density load position
	int  S_ID;			// number name
	float InitDens;		//initial density
	float Temp;			//Load temperature
	float np2c;			// np2c
	int MAXNP;			// maximum np
	float mass;
	float q;

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
*/
/*
#ifndef TotalDomain
typedef struct fluid{
	char name[10];				// string name
	int   Loadtype;				// density load type
	float x_center,x_fall;		// density load position
	float y_center,y_fall;		// density load position
	float InitDens;				// initial density
	float Temp;					// Temperature
	float mass;
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
	char name[10];				// string name
	float Pres;				    // Pressure
	float Temp;				// Temperature
	float mass;
	//float GasDensity;			// density [m^-3]
	//float GasThermalVelocity;	// thermal velocity
	//float m;					// Mass
	//float **den;				// Temperature gadient density
//	float vt;
//	float sig;
} BackG;
#endif
*/



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




