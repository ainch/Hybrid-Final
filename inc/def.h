/******************************************/
/* def.h */
#ifndef __PhysicalConstant__
#define __PhysicalConstant__
#define M_PI 3.14159265358979323846
#define EPS0 8.8542e-12      /* (F/m)  */
#define KB 1.38064852e-23
#define C0 3.0e8
#define SNG_MIN 1e-200
#define DBL_MIN    1E-200
#define NperTORR   8.3221e20
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
enum ParticleLoadInitial {
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

#define  NVTS    3.3




