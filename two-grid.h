/*
 *
 * Two grid MG
 *
 *
 *  */

#include "phg.h"


#define MAX_N_LEVEL 10

typedef struct MG_DATA_ {
  int nlevel;
    GRID *g[MAX_N_LEVEL];
    DOF *dof[MAX_N_LEVEL];
    SOLVER *solver[MAX_N_LEVEL];
    SOLVER *smoother[MAX_N_LEVEL];
    MAT *A[MAX_N_LEVEL];
    MAT *P[MAX_N_LEVEL];
    VEC *x[MAX_N_LEVEL];
    VEC *f[MAX_N_LEVEL];
    VEC *r[MAX_N_LEVEL];
    VEC *fext[MAX_N_LEVEL];
} MG_DATA;

typedef void (*FUNC_BUILD_SYSTEM)(SOLVER *solver, DOF *u_h);


MG_DATA *
mg_create2(GRID *g1,
	   GRID *g2,
	   DOF *u,
	   SOLVER *solver1,
	   char *S1_opts,
	   char *S2_opts,
	   int Cmat_type,
	   FUNC_BUILD_SYSTEM build_system
	   );

MG_DATA *
mg_createN(int nlevel,
	   GRID *g1, 
	   GRID *g2, 
	   GRID *g3, 
	   GRID *g4, 
	   DOF *u,
	   SOLVER *solver1,
	   char *S1_opts,
	   char *S2_opts, 
	   char *S3_opts, 
	   char *S4_opts, 
	   int Cmat_type,
	   FUNC_BUILD_SYSTEM build_system	  
	   );


void mg_solve(MG_DATA *mg, SOLVER *solver, DOF *u);
void mg_pc_proc(void *ctx, VEC *b0, VEC **x0);

void destroy_mg(MG_DATA **mg_ptr);
