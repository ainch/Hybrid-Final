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
#define MAX_VIEW_NP 1000000
#define N_MAX 100
#define THREADS_PER_BLOCK 512   
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
} point;
typedef struct __Global_PCG_Const_Set{ // for cuda_Field_SpeedTest
	int A_size;
	int Iter;
	float tol;
	float tol2;
    float rsold;
    float Temp;
    float rnew;
	float alpha;
	float beta;
}DPS_Const;
typedef struct __Global_PCG_DATA_Set{ // for cuda_Field_SpeedTest
	float R;
	float Z;
	float P;
    float AP;
    float M;
}DPS_Data;
typedef struct __Global_Gsize_Array{
	//SIZE = Gsize
	int Boundary;   //Boundary Condition Constant 0~4
	int CondID;		// Conductor ID , NO Conductor is zero
	int Face;
	int DensRegion;
	float Area;
	float Temp;
	float BackDen1;
	float BackVel1;
	float BackDen2;
	float BackVel2;
	float Lap_Pot;
	float Pois_Pot;
	float Ex;
	float Ey;
}GGA;
typedef struct __Global_Csize_Array{
	//SIZE = Csize
    int PlasmaRegion;
    float eps_r;
	float eps;
}GCA;
typedef struct __Host_Charged_Particle
{ // for save and dump load
	int *CellID;// SIZE(NP_LIM)
    float *x;   // SIZE(NP_LIM)
    float *y;  	// SIZE(NP_LIM)
    float *vx;  // SIZE(NP_LIM)
	float *vy;  // SIZE(NP_LIM)
	float *vz;  // SIZE(NP_LIM)
}HCP;
typedef struct __Global_Charged_Particle
{
	//SIZE = nsp * MAXNP
	int CellID;
    float x;  
    float y;
    float vx;
	float vy;
	float vz;
}GCP;
typedef struct __Global_Particle_Gsize_Data
{
	//SIZE = nsp * Gsize
	int PtNumInCell;
	int PtNumMoveInterCell;
    int MaxPtNumInCell;  
	int PtNumMCCInCell;
	int PtNullMCCInCell;
    float den;
	float smt_den;
    float sum_den;
	float ave_den;
	float sigma;
}GPG;
typedef struct __Global_Info_Particle
{
	//Size = [nsp]
	//Constant
	int spnum;
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
	int St_num;					// Start Address Number 
	int End_num;				// End Address Number
	float mass;
	float q;
	float q_density;
	float vti;
	int np;
	float qm;
	float Escale;
	float Ascale;
	float Denscale;
    float MCCscale;
}Species;
typedef struct fluid{
	//Size =[nfsp]  
	char name[10];				// string name
	int   Loadtype;				// density load type
	float x_center,x_fall;		// density load position
	float y_center,y_fall;		// density load position
	float InitDens;				// initial density
	float Temp;					// Temperature
	float mass;
	float Vel;
	float Gamma1;				// Quenching
	float ave_Den;
    int CSS_Flag;				// Continuity equation on or off
	int CSS_Check;				// Cycle to judge convergence
	float CSS_Conver;           // convergence judge value
}Fluid;
typedef struct __Global_Fluid_Csize_Data
{
	float **D;     	
    float **den;		
	float **Source;	
    float **gummel_ax;			// Gummel method coefficient for x direction flux
	float **gummel_bx;			// Gummel method coefficient for x direction flux
	float **gummel_ay;			// Gummel method coefficient for y direction flux
	float **gummel_by;			// Gummel method coefficient for y direction flux
	float **flux_x;				// x direction flux
	float **flux_y;				// y direction flux
}GFC;
typedef struct __Global_Fluid_Gsize_Data
{
	//SIZE = nfsp * Gsize, real size is nfsp * [ngx*ncy]
	float n;
}GFG;
typedef struct Continuity_region_x {
	// SIZE = nfsp * Conti_xnum
	int *x1;
	int *x2;
	int *yy;		
	float *fg1;
	float *fg2;	
}Con_RegionX;
typedef struct Continuity_region_y {
	// SIZE = nfsp * Conti_ynum
	int *xx;
	int *y1;
	int *y2;	
	float *fg1;
	float *fg2;	
}Con_RegionY;
typedef struct BackGround{
	char name[10];				// string name
	float Ratio;
	float Pres;				    // Pressure
	float Temp;				// Temperature
	float mass;
	float InitDens;
}BackG;
typedef struct __Host_History{
	//SIZE = nsp or nfsp
    float *np;
}Hist;
typedef struct __Global_MCC_sigmav
{
	// Size[] = Ar:3, O2:20, ArO2=25
	float val; // value
}MCC_sigmav;
typedef struct __Global_Collision_Flag
{
	// Size[TnRct], TnRct = Total Number of Reaction
    float Flag; // CX = Flag * CX
	float mofM; // target/projectile
	float Th_e; //Threshold energy
	float RR; // Reaction coefficient
} CollF;
typedef struct __Global_Argon_Collision_Data
{
	// Size[N_LOGX], N_LOGX = Number of data
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
typedef struct __Global_Oxygen_Collision_Data
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
typedef struct __Global_ArO2_Collision_Data
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



