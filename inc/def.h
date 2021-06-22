/******************************************/
/* def.h */
#define SNG_MIN 1e-200

#define DECIMAL float
#define DECIMAL_TYPE  "%g"
#define HIST  1000000

#define M_PI 3.14159265358979323846
#define EPS0            8.8542e-12      /* (F/m)  */
#define KB 1.38064852e-23
#define C0 3.0e8

#define SINGLEGPU

#define XYPIC	0
#define RZPIC	1

#define RZ_X_LIM 200

#define DIRICHLET 1
#define NEUMANN 2
#define CONDUCTOR 3
#define DIELECTRIC 4
#define DIRICHLET_LIN 5

#define LEFT    0
#define RIGHT   1 
#define UP      2
#define DOWN    3
#define UL_CORN 4
#define UR_CORN 5
#define LL_CORN 6
#define LR_CORN 7
#define NO_FACE 8

#define QUIET_START	0
#define RANDOM_XY	1
#define RANDOM_RZ	2

#define  NVTS    3.3

#define False 0
#define True 1

#define DBL_MIN    1E-200
#define NperTORR   8.3221e20

#ifndef max
#define max(x,y) (((x)>(y)) ? (x) : (y))
#endif

#ifndef min
#define min(x,y) (((x)<(y)) ? (x) : (y))
#endif

#define ecolsp  0
#define ionsp   1

#define ARGON 0
#define OXYGEN 1
#define ARO2 2