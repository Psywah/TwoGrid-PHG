/*
 *
 * Two grid solver
 *
 *  (highlight-regexp "[a-zA-Z]+1" 'hi-yellow)
 *  (highlight-regexp "[a-zA-Z]+2" 'hi-green)
 *
 *  */
#include "two-grid.h"
#include "oct-search.h"
#include "petscmat.h" 


/* for PETSc 3.3 changes */
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=3) || PETSC_VERSION_MAJOR>3
# define MatCreateMPIAIJ        MatCreateAIJ
# define MatCreateMPIBAIJ       MatCreateBAIJ
# define MatCreateMPISBAIJ      MatCreateSBAIJ
# define MatCreateMPIDense      MatCreateDense
#endif	/* PETSC_VERSION_MAJOR==3 && ... */


/* mg params */
static int mg_n_pre_smooth = 1;
static int mg_n_post_smooth = 1;
static int mg_n_coarst_smooth = 1;
static FLOAT mg_correct_damp = 1.;


#define smooth(A, x, f0, nsmooth)	{	\
	smoother->maxit = nsmooth;		\
	phgVecCopy(f0, &smoother->rhs);		\
	phgSolverVecSolve(smoother, FALSE, x);	\
    }

#define ZeroVec(x) phgVecAXPBY(0., NULL, 0., &(x));

#define TIMING_BEGIN if (1) {			\
        time0 = phgGetTime(NULL);               \
    }
#define TIMING_END(desp) if (1) {		\
        phgInfo(0, "   L:%d, T:%0.8lfs, %s\n",  \
                level,				\
                phgGetTime(NULL) - time0,	\
                desp);				\
    }


void
proj_vec(MAT_OP op, FLOAT alpha, MAT *A, VEC *x, FLOAT beta, VEC *y, VEC *fext)
{
    int i;
    
    if (op == MAT_OP_T) {  	/* F2C */
	/* fext = A^T * x  */
	phgMatVec(op, alpha, A, x, 0., &fext);
 
	/* fext -> y */
	if (y != NULL) {
	    for (i = 0; i < y->map->nlocal; i++)
		y->data[i] = beta * y->data[i] + fext->data[i];
	}
    }
    else if (op == MAT_OP_N) {	/* C2F */
	/* x-> fext */
	if (x != NULL)
	    for (i = 0; i < x->map->nlocal; i++)
		fext->data[i] = x->data[i];
	
	/* y = A * fext  */
	phgMatVec(op, alpha, A, fext, beta, &y);	
    }
    else {
	phgError(1, "   Unkonwn MAT op type!!!\n");
    }
}


static FLOAT time0, time1;

void
mg_cycle(MG_DATA *mg, int level, BOOLEAN init_zero)
{
  int max_level = mg->nlevel;
  GRID *g = mg->g[level];
    SOLVER *smoother = mg->smoother[level];
    MAT *A = mg->A[level];
    MAT *P = mg->P[level];
    VEC *x = mg->x[level];
    VEC *f = mg->f[level];
    VEC *r = mg->r[level];
    int i;
 
    if (smoother == NULL) 	/* empty rank */
	return;

    if (level == max_level - 1) { /* timing on finest level */
	MPI_Barrier(g->comm);
        time1 = phgGetTime(NULL);
	phgInfo(0, "\n\n");
    }
    
    if (init_zero)
	ZeroVec(x);

    /*
     * coarsest level solve
     * */
    if (level <= 0) {
	VEC *f0 = phgVecCopy(f, NULL);

	assert(level >= 0);

	/* use smoother on coarsest grid */
	TIMING_BEGIN;
	smooth(A, x, f0, mg_n_coarst_smooth);
	TIMING_END("Coarst smooth");

	phgVecDestroy(&f0);
    }
    /*
     * intermeidate level solve
     * */
    else {
	VEC *f0 = NULL;
	f0 = phgVecCopy(f, NULL);

	/* pre-smooth */
	TIMING_BEGIN;
	smooth(A, x, f0, mg_n_pre_smooth);
	TIMING_END("Pre smooth");

#if 1
	/* Enable coarse grid correction */

	
	/* residual:
	 *   r_{k} = b_{k} - A_{k} * x_{k} */
	phgVecCopy(f0, &r);
	TIMING_BEGIN;
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
	TIMING_END("Residual");

	/* restriction:
	 *   f_{k-1} = P_{k} r_{k} */
	TIMING_BEGIN;
	proj_vec(MAT_OP_T, 1.0, P, r, 0., mg->f[level-1],
		 mg->fext[level-1]);
	TIMING_END("restrict");
	
	/* recursive cycle
	 * n_cycle = 1: V-cycle
	 * n_cycle = 2: W-cycle
	 */
	const int n_cycle = 1;
	for (i = 0; i < n_cycle; i++) {
	    /* Note:
	     *   Not sure in smoothing as done in Alberta is neccessory.
	     *   Without this in smoothing, W(+) cycle may not improve none
	     *   compare to V cycle.
	     * */
	    mg_cycle(mg, level-1, (i==0));
	}


	/* prolongation:
	 *   x_{k} += P_{k}^T x_{k-1} */
	TIMING_BEGIN;
	proj_vec(MAT_OP_N, mg_correct_damp,
		 P, mg->x[level-1], 1.0, x,
		 mg->fext[level-1]);
	TIMING_END("prolong");
#endif

	/* Post smooth */
	TIMING_BEGIN;
	smooth(A, x, f0, mg_n_post_smooth);
	TIMING_END("Post smooth");

	phgVecDestroy(&f0);
    }


    MPI_Barrier(g->comm);
    if (g->rank == 0)
	phgInfo(0, "   Cycle[%2d] cost: %0.8lfs\n",
		level, phgGetTime(NULL) - time1);

    if (level == max_level - 1) {
	MPI_Barrier(g->comm);
        if (g->rank == 0)
	    phgInfo(0, "   Cycle cost: %0.8lfs\n",
		    phgGetTime(NULL) - time1);
    }
    return;
}


void
mg_pc_proc(void *ctx, VEC *b0, VEC **x0)
{
    SOLVER *pc = (SOLVER *)ctx;
    MG_DATA *mg = pc->mat->mv_data[0];
    MAT *A = pc->mat;
    int level = 1;
    VEC *x, *r;

    r = phgMapCreateVec(A->rmap, 1);
    x = phgMapCreateVec(A->rmap, 1);
    phgVecCopy(b0, &r);
    phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);

    /* MG cycle setup */
    phgVecCopy(r, &mg->f[level]);

    /* MG correct on x */
    mg_cycle(mg, 1, TRUE);	/* two level, init zero */

    /* return smoother solution */
    phgVecCopy(mg->x[level], x0);

    phgVecDestroy(&r);
    phgVecDestroy(&x);
    return;
}



void
mg_solve(MG_DATA *mg, SOLVER *solver, DOF *u)
{
    FLOAT r0_norm, r_norm, b_norm;
    int nits = 0, maxit = solver->maxit,
	level, max_level = mg->nlevel;
    MAT *A = solver->mat;
    VEC *x, *r, *b0 = solver->rhs;


    /* work on finest level */
    level = max_level - 1; 
    x = mg->x[level];
    r = phgMapCreateVec(solver->rhs->map, 1);
    phgMapDofToLocalData(solver->rhs->map, 1, &u, x->data);


    /* initial residual */
    r_norm = 0.;
    b_norm = phgVecNorm2(b0, 0, NULL);
    phgVecCopy(b0, &r);
    phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r); 	
    r0_norm = phgVecNorm2(r, 0, NULL);
    if (solver->verb > 1)
	phgPrintf("   initial err %E\n", r0_norm);
    //phgVecDumpMATLAB(r, "r0", "r0_.m");

    /* MG cycle */
    while (nits < maxit) {
	/* MG cycle setup */
	phgVecCopy(b0, &mg->f[level]);

	/* MG correct on x */
	mg_cycle(mg, level, FALSE);
	    
	/* new residual */
	phgVecCopy(b0, &r);
	phgMatVec(MAT_OP_N, -1.0, A, x, 1.0, &r);
	r_norm = phgVecNorm2(r, 0, NULL);
	phgPrintf("   niter: %d, err %E\n",
		  nits, r_norm / r0_norm);

	nits++;
	if ((r_norm / r0_norm) < solver->rtol)
	    break;
    }

    solver->nits = nits;
    solver->residual = r_norm;
    if (nits < maxit) {
	phgPrintf("   MG converge at step: %d\n", nits);
    } else
	phgPrintf("WARNING: maxit attained in OEM solver, "
		  "the result may be inaccurate.\n");

    if (0) {
	VEC *solu = phgMapCreateVec(solver->rhs->map, 1);
	phgMapDofToLocalData(solver->rhs->map, 1, &u, solu->data);
	phgVecDumpMATLAB(solu, "x0", "x0_.m");
	phgVecDumpMATLAB(x, "x", "x_.m");
	phgVecDestroy(&solu);
    }

    phgMapLocalDataToDof(solver->rhs->map, 1, &u, x->data);
    phgVecDestroy(&r);

    return;
}






typedef struct PT_INFO_ {
    int rank;			/* rank */
    INT id;			/* vector index */
    FLOAT x[Dim];		/* xyz */
} PT_INFO;


static void
func_xyz(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
    *(v++) = x;
    *(v++) = y;
    *(v++) = z;
}    


#define TEST_ORDER 2
static void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
#if TEST_ORDER == 1    
    *v = x + 2*y - 3*z;
#else
    *v = x*x + 2*y*y - 3*z*z;
#endif
}    


/* ----------------------------------------
 * 
 * Two grid interp Mat,
 *    for one dof_type
 *
 * ---------------------------------------- */
static MAT *
build_interp_mat(GRID *g, GRID *g2,
		 _OCT_TREE *og, DOF_TYPE *dof_type)
{
    SIMPLEX *e1, *e2;
    DOF *u1, *u2, *coord1;;
    MAP *map1, *map2, *map1_coord, *map2_ext;
    VEC *vec1, *vec2, *vec1_coord, *vec2_ext;
    MAT *matP;
    VEC *vecP;
    MPI_Datatype type;
    FLOAT *xyz, lambda[Dim+1];
    int i, j, k, rank;
    FLOAT (*bboxs)[Dim][2];
    int nprocs2;
    //MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm comm = g->comm;

    /* bcast bboxs of coars grid */
    if (phgRank == 0)
	nprocs2 = g2->nprocs;
    MPI_Bcast(&nprocs2, 1, MPI_INT, 0, comm);
    bboxs = phgCalloc(nprocs2, sizeof(*bboxs));
    if (phgRank == 0)
	memcpy(&bboxs[0][0][0], &og->bboxs[0][0][0], nprocs2 * sizeof(*bboxs));
    MPI_Bcast(&bboxs[0][0][0],
	      nprocs2 * Dim * 2, PHG_MPI_FLOAT, 0, comm);

    for (rank = 0; rank < nprocs2; rank++) {
	phgInfo(0, "rank box %d: [%12.4f, %12.4f] x [%12.4f %12.4f] x [%12.4f %12.4f]\n",
		rank,
		bboxs[rank][0][0],
		bboxs[rank][0][1],
		bboxs[rank][1][0],
		bboxs[rank][1][1],
		bboxs[rank][2][0],
		bboxs[rank][2][1]
		);
    }
	
    
    
    /* DOF1 */
    u1 = phgDofNew(g, dof_type, 1, "u_h", DofNoAction);
    map1 = phgMapCreate(u1, NULL);
    vec1 = phgMapCreateVec(map1, 1);

    /* DOF2 */
    if (g2) {
	u2 = phgDofNew(g2, dof_type, 1, "u2_h", DofNoAction);
	map2 = phgMapCreate(u2, NULL);
	vec2 = phgMapCreateVec(map2, 1);
    }
    /* extend map2 */
    {
	INT nglobal;

	if (phgRank == 0)
	    nglobal = map2->nglobal;
	MPI_Bcast(&nglobal, 1, PHG_MPI_INT, 0, comm);
	map2_ext = phgMapCreateSimpleMap(comm,
					 (g2) ? map2->nlocal : 0,
					 nglobal);
	vec2_ext = phgMapCreateVec(map2_ext, 1);
    }



    /* Mat projection */
    matP = phgMapCreateMat(map1, map2_ext);
    vecP = phgMapCreateVec(map1, 1);
    phgVecDisassemble(vecP);
    phgInfo(0, "map1: %d\n", map1->nlocal);

    /* coord */
    coord1 = phgDofNew(g, dof_type, Dim, "coord1", func_xyz);
    map1_coord = phgMapCreate(coord1, NULL);
    vec1_coord = phgMapCreateVec(map1_coord, 1);
    phgMapDofToLocalData(map1_coord, 1, &coord1, vec1_coord->data);


    /* search */
    BOOLEAN *flag;
    INT nsend, nrecv, //*counts, *displs,
	*scnts, *sdsps, *rcnts, *rdsps;
	
    flag = phgCalloc(map1->nlocal, sizeof(*flag));

    xyz = vec1_coord->data;
    for (i = 0; i < map1->nlocal; i++) { /* vector index */
	BOOLEAN found;
	int iG = map1->partition[map1->rank] + i;    /* v2G */
	phgInfo(3, "dof coord: %12.8e, %12.8e, %12.8e\n", 
		xyz[0], xyz[1], xyz[2]);

	/* local */
	if (og)
	    found = _locate_point(og, xyz, &e2, lambda);
	else
	    found = FALSE;
	if (found) {
	    int N = u2->type->nbas;	
	    INT IG[N];	/* global index */
	    FLOAT *bas, one = 1.;

	    for (j = 0; j < N; j++) 
		IG[j] = phgMapE2G(map2, 0, e2, j);

	    bas = (FLOAT *) u2->type->BasFuncs(u2, e2, 0, -1, lambda);
	    phgMatAddGlobalEntries(matP, 1, &iG, N, IG, bas);
	    phgVecAddGlobalEntry(vecP, 0, iG, one);

	    phgInfo(3, "g1: %d\n", iG);
	    for (j = 0; j < N; j++)
		phgInfo(3, "   %d %12.8e\n", IG[j], bas[j]);

	    flag[i] = TRUE;	/* added */
	} 

	xyz += Dim;
    }
	

    /*
     *
     * send search info
     *
     * */
    scnts = phgCalloc(4 * g->nprocs, sizeof(*scnts));
    //scnts_val = phgCalloc(g->nprocs, sizeof(*scnts));
    sdsps = scnts + g->nprocs;
    rcnts = sdsps + g->nprocs;
    rdsps = rcnts + g->nprocs;
    //counts = phgAlloc(g->nprocs * sizeof(*counts));
    //displs = phgAlloc(g->nprocs * sizeof(*displs));

    /* count send */
    xyz = vec1_coord->data;
    for (i = 0; i < map1->nlocal; i++) { /* vector index */
	if (!flag[i]) {
		    
	    for (rank = 0; rank < nprocs2; rank++) {
		if (rank != g->rank
		    && bboxs[rank][0][0] - EPS < xyz[0]
		    && xyz[0] < bboxs[rank][0][1] + EPS 
		    && bboxs[rank][1][0] - EPS < xyz[1]
		    && xyz[1] < bboxs[rank][1][1] + EPS 
		    && bboxs[rank][2][0] - EPS < xyz[2]
		    && xyz[2] < bboxs[rank][2][1] + EPS) {

		    scnts[rank]++;
		}
	    }
	}
	xyz += Dim;
    }    

    MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, comm);
    nsend = nrecv = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
	sdsps[rank] = nsend;
	rdsps[rank] = nrecv;
	nsend += scnts[rank];
	nrecv += rcnts[rank];
    }

    PT_INFO *sbuf_pt, *rbuf_pt, *pt_info;
    sbuf_pt = phgCalloc(nsend + nrecv, sizeof(*sbuf_pt));
    rbuf_pt = sbuf_pt + nsend;

    /* fill send buf */
    xyz = vec1_coord->data;
    for (i = 0; i < map1->nlocal; i++) { /* vector index */
	if (!flag[i]) {
	    int iG = map1->partition[map1->rank] + i;    /* v2G */
		    
	    for (rank = 0; rank < nprocs2; rank++) {
		if (rank != g->rank
		    && bboxs[rank][0][0] - EPS < xyz[0]
		    && xyz[0] < bboxs[rank][0][1] + EPS 
		    && bboxs[rank][1][0] - EPS < xyz[1]
		    && xyz[1] < bboxs[rank][1][1] + EPS 
		    && bboxs[rank][2][0] - EPS < xyz[2]
		    && xyz[2] < bboxs[rank][2][1] + EPS) {
			
		    pt_info = sbuf_pt + sdsps[rank]++;
		    pt_info->rank = g->rank;
		    pt_info->id = iG;
		    memcpy(pt_info->x, xyz, Dim * sizeof(FLOAT));
		}
	    }
	}
	xyz += Dim;
    }    
	
    /* restore face sdsps[] */
    for (rank = 0; rank < g->nprocs; rank++)
	sdsps[rank] -= scnts[rank];
    

    /* alltoallv */
    MPI_Type_contiguous(sizeof(*sbuf_pt), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuf_pt, scnts, sdsps, type, 
		  rbuf_pt, rcnts, rdsps, type, comm);
    MPI_Type_free(&type);


    /*
     *
     * step 3. remote local search
     *
     * */
    if (og == NULL) 
	assert(nrecv == 0);

    pt_info = rbuf_pt;
    for (i = 0; i < nrecv; i++, pt_info++) {
	BOOLEAN found;
	INT iG = pt_info->id;

	phgInfo(3, "to find %12.8e %d %d\n", 
		pt_info->x[0], pt_info->rank, pt_info->id);
	found = _locate_point(og, pt_info->x, &e2, lambda);

	if (found) {
	    int N = u2->type->nbas;
	    INT IG[N];	/* global index */
	    FLOAT *bas;

	    for (j = 0; j < N; j++)
		IG[j] = phgMapE2G(map2, 0, e2, j);

	    bas = (FLOAT *) u2->type->BasFuncs(u2, e2, 0, -1, lambda);
	    phgMatAddGlobalEntries(matP, 1, &iG, /* add to remote row */
				   N, IG, bas);
	    /* phgVecAddGlobalEntry(vecP, 0, iG, one); */
	    /* phgInfo(0, "vecP: %x %d\n", vecP->O2Gmap); */

	    phgInfo(3, "g1: %d\n", iG);
	    for (j = 0; j < N; j++)
		phgInfo(3, "   %d %12.8e\n", IG[j], bas[j]);

	    pt_info->id = -1; /* marked as found */
	}
    }


    /*
     * Send back found info.
     *  */
    /* alltoallv */
    MPI_Type_contiguous(sizeof(*sbuf_pt), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(rbuf_pt, rcnts, rdsps, type,
		  sbuf_pt, scnts, sdsps, type, comm);
    MPI_Type_free(&type);

	
    xyz = vec1_coord->data;
    for (i = 0; i < map1->nlocal; i++) { /* vector index */
	if (!flag[i]) {
	    int iG = map1->partition[map1->rank] + i;    /* v2G */
		    
	    int found = 0;
	    for (rank = 0; rank < nprocs2; rank++) {
		if (rank != g->rank
		    && bboxs[rank][0][0] - EPS < xyz[0]
		    && xyz[0] < bboxs[rank][0][1] + EPS 
		    && bboxs[rank][1][0] - EPS < xyz[1]
		    && xyz[1] < bboxs[rank][1][1] + EPS 
		    && bboxs[rank][2][0] - EPS < xyz[2]
		    && xyz[2] < bboxs[rank][2][1] + EPS) {
			
		    pt_info = sbuf_pt + sdsps[rank]++;
		    if (pt_info->id == -1)
			found++;
		}
	    }
	    if (found == 0) {
		phgError(1, "can't find point (%12.8e, %12.8e, %12.8e)\n", 
			 pt_info->x[0], 
			 pt_info->x[1], 
			 pt_info->x[2]
			 );
	    }

	    FLOAT value = found;
	    phgVecAddGlobalEntry(vecP, 0, iG, value);
	}
	xyz += Dim;
    }    

    /* restore face sdsps[] */
    for (rank = 0; rank < g->nprocs; rank++)
	sdsps[rank] -= scnts[rank];




    /* Assemble and rescale */
    phgMatAssemble(matP);	/* Costly !!! */
    //phgMatDumpMATLAB(matP, "P0", "P0_.m");
    phgVecAssemble(vecP);
    //phgVecDumpMATLAB(vecP, "vP", "vP_.m");
    if (1) {
	MAT_ROW *row = matP->rows;
	const FLOAT *v = vecP->data;
	FLOAT sum;

	assert(matP->type == PHG_UNPACKED);
	for (i = 0; i < matP->rmap->nlocal; i++) {
	    INT iG = map1->partition[map1->rank] + i;
	    assert(*v > .5);
	    phgInfo(3, "row: %d, scale: %e\n", 
		    iG, *v);

	    sum = 0.;
	    for (j = 0; j < row->ncols; j++) {
		phgInfo(3, "   col: %d, %e\n",
			row->cols[j], row->data[j]);
		row->data[j] /= *v;
		sum += row->data[j];
	    }
	    if (fabs(sum - 1.) > 1e-12)
		phgError(1, "   row[%d] sum %e not one !\n",
			 iG, sum);

	    row++;
	    v++;
	}
    }
	


    phgFree(flag);
    phgFree(scnts);
    phgFree(sbuf_pt);

	
#if 0	
    /* Dump mat */
    phgMatDumpMATLAB(matP, "P", "P_.m");
#endif


#if 1	
    /* Function test */
    phgDofSetDataByFunction(u1, func_u);
    phgMapDofToLocalData(map1, 1, &u1, vec1->data);
    if (g2) {
	phgDofSetDataByFunction(u2, func_u);
	phgMapDofToLocalData(map2, 1, &u2, vec2->data);
	memcpy(vec2_ext->data, vec2->data,
	       map2->nlocal * sizeof(*vec2->data));
    }

#  if TEST_ORDER == 1
    assert(u1->type->order == 1);
#  else
    assert(u1->type->order > 1);
#  endif
	
    //phgVecDumpMATLAB(vec1, "v1", "v1_.m");
    //phgVecDumpMATLAB(vec2, "v2", "v2_.m");
    phgMatVec(MAT_OP_N, -1., matP, vec2_ext, 1., &vec1);
    //phgVecDumpMATLAB(vec1, "v1", "r_.m");
    phgPrintf("\n\n\n---\nMat projection error: %e\n", phgVecNorm2(vec1, 0, NULL));
#endif





    phgVecDestroy(&vecP);

    phgVecDestroy(&vec1_coord);
    phgMapDestroy(&map1_coord);
    phgDofFree(&coord1);

    phgVecDestroy(&vec1);
    phgMapDestroy(&map1);
    phgDofFree(&u1);

    if (g2) {
	phgVecDestroy(&vec2);
	phgMapDestroy(&map2);
	phgDofFree(&u2);
    }
    phgVecDestroy(&vec2_ext);
    phgMapDestroy(&map2_ext);

    phgFree(bboxs);
    
    return matP;
}









/* Petsc mat utils */
Mat 
mat_phg2petsc(MAT *mat)
{
    MAP *rmap = mat->rmap;
    MAP *cmap = mat->cmap;
    PetscInt N = rmap->nglobal, Nlocal = rmap->nlocal;
    PetscInt M = cmap->nglobal, Mlocal = cmap->nlocal;
    PetscInt i, j, k, prealloc, *d_nnz, *o_nnz;
    PetscInt jstart, jend;
    INT d, o;
    Mat A;

    prealloc = PETSC_DECIDE;
    d_nnz = phgAlloc(Nlocal * sizeof(*d_nnz));
    o_nnz = phgAlloc(Nlocal * sizeof(*o_nnz));
#if 0
    /* phgMatGetNnzInRow not right for non-square mat */
    for (i = 0; i < Nlocal; i++) {
	phgMatGetNnzInRow(mat, i, &d, &o);
	d_nnz[i] = d;
	o_nnz[i] = o;
    }
#else
    jstart = cmap->partition[cmap->rank];
    jend   = cmap->partition[cmap->rank + 1];
    j = rmap->partition[rmap->rank];
    for (i = 0; i < rmap->nlocal; i++, j++) {
	const MAT_ROW *row = phgMatGetRow(mat, i);
	int ncols = row->ncols;
	INT *cols  = row->cols;
	FLOAT *values = row->data;
	int d, o;

	d = o = 0;
	for (k = 0; k < ncols; k++) {
	    if (row->cols[k] >= jstart
		&& row->cols[k] < jend)
		d++;
	    else
		o++;
	}
	d_nnz[i] = d;
	o_nnz[i] = o;
    }
#endif
    MatCreateMPIAIJ(rmap->comm, Nlocal, Mlocal, N, M,
		 prealloc, d_nnz, prealloc, o_nnz, &A);


    j = rmap->partition[rmap->rank];
    for (i = 0; i < rmap->nlocal; i++, j++) {
	const MAT_ROW *row = phgMatGetRow(mat, i);
	int ncols = row->ncols;
	INT *cols  = row->cols;
	FLOAT *values = row->data;
    
	PetscInt n = ncols;
	PetscScalar *buffer = phgAlloc(n * sizeof(*buffer));
	for (k = 0; k < n; k++)
	    buffer[k] = (PetscScalar)values[k];
	MatSetValues(A, 1, &j, ncols, cols,
		     buffer, INSERT_VALUES);
	phgFree(buffer);

	phgInfo(3, "   row %5d; %3d %3d %3d\n", i, ncols, d_nnz[i], o_nnz[i]);
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); 
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); 

    
    phgFree(d_nnz);
    phgFree(o_nnz);
    return A;
}



MAT *
mat_petsc2phg(Mat A, MAP *map)
{
    int N, M, Nlocal, Mlocal;
    PetscInt          m,n,rstart,rend;
    PetscInt          row,ncols,j,nrows;
    const PetscInt    *cols;
    const PetscScalar *vals;
    MAT *mat;
    
    MatGetLocalSize(A, &Nlocal, &Mlocal);
    MatGetSize(A, &N, &M);
    assert(map->nlocal == Nlocal);;
    assert(map->nglobal == N);;

    mat = phgMapCreateMat(map, map);
    
    MatGetOwnershipRange(A, &rstart, &rend);
    nrows = 0;
    for (row = rstart; row < rend; row++) {
	MatGetRow(A, row, &ncols, &cols, &vals);

	FLOAT *buffer = phgAlloc(ncols * sizeof(*buffer));
	phgInfo(3, "row: %d\n", row);
	for (j = 0; j < ncols; j++) {
	    phgInfo(3, "   col: %d, %e\n", cols[j], vals[j]);
	    buffer[j] = vals[j];
	}
	phgMatAddGlobalEntries(mat, 1, &row, ncols, cols, buffer);
	phgFree(buffer);
	
	MatRestoreRow(A, row, &ncols, &cols, &vals);
    }    

    phgMatAssemble(mat);
    
    return mat;
}




/* ----------------------------------------
 * 
 * Mg create 2
 *
 * ---------------------------------------- */
MG_DATA *
mg_create2(GRID *g1,
	  GRID *g2,		/* NULL if not in comm */
	  DOF *u,
	  SOLVER *solver1,
	  char *S1_opts,
	  char *S2_opts, 
	  int Cmat_type,
	  FUNC_BUILD_SYSTEM build_system	  
	  )
{
    MG_DATA *mg = phgCalloc(1, sizeof(MG_DATA));
    _OCT_TREE *og = NULL;
    SOLVER *solver2;
    DOF *u2;
    MAT *matP, *A2;
    
    mg->g[1] = g1;
    mg->g[0] = g2;

    /* Build interp mat */
    if (g2)
	og = _build_octree(g2);
    else
	og = NULL;

    /* { */
    /* 	FLOAT X[3] = {0.25, -0.25, 0.5}; */
    /* 	BOOLEAN found = _locate_point(og, X, &e2, lambda); */
    /* 	phgInfo(0, "Test special found: %d\n", found); */
    /* } */

    
    time0 = phgGetTime(NULL);
    matP = build_interp_mat(g1, g2,
			    og, u->type); /* might have more */
    phgPrintf("Build interp mat: %lfs\n", phgGetTime(NULL) - time0);
    

    if (og)
      _destroy_octree(&og);	/* free oct tree */

    

    
    /* fine level */
    mg->A[1] = solver1->mat;
    mg->P[1] = matP;
    mg->x[1] = phgVecCopy(solver1->rhs, NULL);
    mg->f[1] = phgVecCopy(solver1->rhs, NULL);;
    mg->r[1] = phgVecCopy(solver1->rhs, NULL);;
    phgOptionsPush();
    phgOptionsSetOptions(S1_opts);
    mg->smoother[1] = phgMat2Solver(SOLVER_DEFAULT, mg->A[1]);
    mg->smoother[1]->warn_maxit = FALSE;
    mg->smoother[1]->verb = 0;
    phgOptionsPop();



    /* coarse level */
    if (g2) {
	u2 = phgDofNew(g2, u->type, u->dim, "u2", func_u);
	time0 = phgGetTime(NULL);
	solver2 = phgSolverCreate(SOLVER_DEFAULT, u2, NULL);
	phgPrintf("Build mat time: %lfs\n", phgGetTime(NULL) - time0);
    }
    else {
	solver2 = NULL;
	u2 = NULL;
    }
	    
    if (Cmat_type == 0)
	phgPrintf("\n\n\n---\n   A2 = L_h2\n");
    else if (Cmat_type == 1)
	phgPrintf("\n\n\n---\n   A2 = PtAP\n");
    else if (Cmat_type == 2)
	phgPrintf("\n\n\n---\n   A2 = PtAP with fill-in of L_h2 \n");

    /* Discreate at coarse level */
    if (g2) {
	build_system(solver2, u2);
	A2 = solver2->mat;	/* default: discreate A2 */
    }


    if (Cmat_type >= 1) {
	/* A = P^t A P */
	Mat A, P, C;
	MAT *mat;
	A = mat_phg2petsc(solver1->mat);
	P = mat_phg2petsc(matP);
		
	MatPtAP(A, P,
		MAT_INITIAL_MATRIX,
		2.0, /* scall estimate */
		&C);

	if (g2) {
	    mat = mat_petsc2phg(C, solver2->mat->rmap);

	    //phgMatDumpMATLAB(mat, "A2", "A2_.m");
	    //phgMatDumpMATLAB(solver2->mat, "A3", "A3_.m");
	}
	else {
	    mat = NULL;
	}

	if (Cmat_type == 1)	/* Use PtAP */
	    A2 = mat;

	if (Cmat_type == 2) {	/* Use PtAP with fill-in of L_h */
	    /* Not done yet */
	}
    }


    
    
    if (g2) {
	mg->dof[0] = u2;
	mg->solver[0] = solver2;
	mg->A[0] = A2;
	mg->x[0] = phgVecCopy(solver2->rhs, NULL);
	mg->f[0] = phgVecCopy(solver2->rhs, NULL);
	mg->r[0] = phgVecCopy(solver2->rhs, NULL);
	phgOptionsPush();
	phgOptionsSetOptions(S2_opts);
	mg->smoother[0] = phgMat2Solver(SOLVER_DEFAULT, mg->A[0]);
	mg->smoother[0]->warn_maxit = FALSE;
	mg->smoother[0]->verb = 0;
	phgOptionsPop();
    }

    mg->fext[0] = phgMapCreateVec(matP->cmap, 1); /* ext vec */

    
    return mg;
}

/* ----------------------------------------
 * 
 * Mg create N
 *
 * ---------------------------------------- */

MG_DATA *
mg_createN(int nlevel,
	   GRID *g1,
	   GRID *g2,		/* NULL if not in comm */
	   GRID *g3,		/* NULL if not in comm */
	   GRID *g4,		/* NULL if not in comm */
	   DOF *u,
	   SOLVER *solver1,
	   char *S1_opts,
	   char *S2_opts, 
	   char *S3_opts, 
	   char *S4_opts, 
	   int Cmat_type,
	   FUNC_BUILD_SYSTEM build_system	  
	   )
{
  MG_DATA *mg = phgCalloc(1, sizeof(MG_DATA));
  _OCT_TREE *og2 = NULL;
  _OCT_TREE *og3 = NULL;
  _OCT_TREE *og4 = NULL;
  SOLVER *solver2 = NULL;
  SOLVER *solver3 = NULL;
  SOLVER *solver4 = NULL;
  DOF *u2 = NULL, *u3 = NULL, *u4 = NULL;
  MAT *matP1 = NULL, *A2 = NULL;
  MAT *matP2 = NULL, *A3 = NULL;
  MAT *matP3 = NULL, *A4 = NULL;
  int ilevel;

  mg->nlevel = nlevel;

  if (nlevel == 2) {
    mg->g[1] = g1;
    mg->g[0] = g2;
  }
  else if (nlevel == 3) {
    mg->g[2] = g1;
    mg->g[1] = g2;
    mg->g[0] = g3;
  }
  else if (nlevel == 4) {
    mg->g[3] = g1;
    mg->g[2] = g2;
    mg->g[1] = g3;
    mg->g[0] = g4;
  }

  /* Build interp mat */
  if (g2)
    og2 = _build_octree(g2);
  else
    og2 = NULL;

  if (nlevel >= 3) {
    if (g3)
      og3 = _build_octree(g3);
    else
      og3 = NULL;
  }

  if (nlevel >= 4) {
    if (g4)
      og4 = _build_octree(g4);
    else
      og4 = NULL;
  }

  /* { */
  /* 	FLOAT X[3] = {0.25, -0.25, 0.5}; */
  /* 	BOOLEAN found = _locate_point(og, X, &e2, lambda); */
  /* 	phgInfo(0, "Test special found: %d\n", found); */
  /* } */

    
  time0 = phgGetTime(NULL);
  matP1 = build_interp_mat(g1, g2,
			  og2, u->type); /* might have more */
  phgPrintf("Build interp mat: %lfs\n", phgGetTime(NULL) - time0);
  if (og2)
    _destroy_octree(&og2);	/* free oct tree */

  if (nlevel >= 3) {
    time0 = phgGetTime(NULL);
    matP2 = build_interp_mat(g2, g3,
			    og3, u->type); /* might have more */
    phgPrintf("Build interp mat: %lfs\n", phgGetTime(NULL) - time0);
    if (og3)
      _destroy_octree(&og3);	/* free oct tree */
  }
    
  if (nlevel >= 4) {
    time0 = phgGetTime(NULL);
    matP3 = build_interp_mat(g3, g4,
			    og4, u->type); /* might have more */
    phgPrintf("Build interp mat: %lfs\n", phgGetTime(NULL) - time0);
    if (og4)
      _destroy_octree(&og4);	/* free oct tree */
  }

    
  /* fine level */
  ilevel = nlevel -1 ;
  mg->A[ilevel] = solver1->mat;
  mg->P[ilevel] = matP1;
  mg->x[ilevel] = phgVecCopy(solver1->rhs, NULL);
  mg->f[ilevel] = phgVecCopy(solver1->rhs, NULL);;
  mg->r[ilevel] = phgVecCopy(solver1->rhs, NULL);;
  phgOptionsPush();
  phgOptionsSetOptions(S1_opts);
  mg->smoother[ilevel] = phgMat2Solver(SOLVER_DEFAULT, mg->A[ilevel]);
  mg->smoother[ilevel]->warn_maxit = FALSE;
  mg->smoother[ilevel]->verb = 0;
  phgOptionsPop();
  mg->dof[ilevel] = mg->A[ilevel]->cmap->dofs[0];

	    
  if (Cmat_type == 0)
    phgPrintf("\n\n\n---\n   A2 = L_h2\n");
  else if (Cmat_type == 1)
    phgPrintf("\n\n\n---\n   A2 = PtAP\n");
  else if (Cmat_type == 2)
    phgPrintf("\n\n\n---\n   A2 = PtAP with fill-in of L_h2 \n");


  /* coarse level */
  if (g2) {
    u2 = phgDofNew(g2, u->type, u->dim, "u2", func_u);
    time0 = phgGetTime(NULL);
    solver2 = phgSolverCreate(SOLVER_DEFAULT, u2, NULL);
    phgPrintf("Build mat time: %lfs\n", phgGetTime(NULL) - time0);
  }
  else {
    solver2 = NULL;
    u2 = NULL;
  }

  /* Discreate at coarse level */
  if (g2) {
    build_system(solver2, u2);
    A2 = solver2->mat;	/* default: discreate A2 */
  }


  if (nlevel >= 3) {
    /* coarse level */
    if (g3) {
      u3 = phgDofNew(g3, u->type, u->dim, "u3", func_u);
      time0 = phgGetTime(NULL);
      solver3 = phgSolverCreate(SOLVER_DEFAULT, u3, NULL);
      phgPrintf("Build mat time: %lfs\n", phgGetTime(NULL) - time0);
    }
    else {
      solver3 = NULL;
      u3 = NULL;
    }

    /* Discreate at coarse level */
    if (g3) {
      build_system(solver3, u3);
      A3 = solver3->mat;	/* default: discreate A2 */
    }
  }


  if (nlevel >= 4) {
    /* coarse level */
    if (g4) {
      u4 = phgDofNew(g4, u->type, u->dim, "u4", func_u);
      time0 = phgGetTime(NULL);
      solver4 = phgSolverCreate(SOLVER_DEFAULT, u4, NULL);
      phgPrintf("Build mat time: %lfs\n", phgGetTime(NULL) - time0);
    }
    else {
      solver4 = NULL;
      u4 = NULL;
    }

    /* Discreate at coarse level */
    if (g4) {
      build_system(solver4, u4);
      A4 = solver4->mat;	/* default: discreate A2 */
    }
  }



  assert(Cmat_type == 0);
  /* if (Cmat_type >= 1) { */
  /*   /\* A = P^t A P *\/ */
  /*   Mat A, P, C; */
  /*   MAT *mat; */
  /*   A = mat_phg2petsc(solver1->mat); */
  /*   P = mat_phg2petsc(matP); */
		
  /*   MatPtAP(A, P, */
  /* 	    MAT_INITIAL_MATRIX, */
  /* 	    2.0, /\* scall estimate *\/ */
  /* 	    &C); */

  /*   if (g2) { */
  /*     mat = mat_petsc2phg(C, solver2->mat->rmap); */

  /*     //phgMatDumpMATLAB(mat, "A2", "A2_.m"); */
  /*     //phgMatDumpMATLAB(solver2->mat, "A3", "A3_.m"); */
  /*   } */
  /*   else { */
  /*     mat = NULL; */
  /*   } */

  /*   if (Cmat_type == 1)	/\* Use PtAP *\/ */
  /*     A2 = mat; */

  /*   if (Cmat_type == 2) {	/\* Use PtAP with fill-in of L_h *\/ */
  /*     /\* Not done yet *\/ */
  /*   } */
  /* } */


    
    

  ilevel = nlevel - 2;
  mg->fext[ilevel] = phgMapCreateVec(matP1->cmap, 1); /* ext vec */
  if (g2) {
    mg->dof[ilevel] = u2;
    mg->solver[ilevel] = solver2;
    mg->A[ilevel] = A2;
    mg->P[ilevel] = matP2;
    mg->x[ilevel] = phgVecCopy(solver2->rhs, NULL);
    mg->f[ilevel] = phgVecCopy(solver2->rhs, NULL);
    mg->r[ilevel] = phgVecCopy(solver2->rhs, NULL);
    phgOptionsPush();
    phgOptionsSetOptions(S2_opts);
    mg->smoother[ilevel] = phgMat2Solver(SOLVER_DEFAULT, mg->A[ilevel]);
    mg->smoother[ilevel]->warn_maxit = FALSE;
    mg->smoother[ilevel]->verb = 0;
    phgOptionsPop();
  }


  if (nlevel >= 3) {

    ilevel = nlevel - 3;
    mg->fext[ilevel] = phgMapCreateVec(matP2->cmap, 1); /* ext vec */
    if (g3) {
      mg->dof[ilevel] = u3;
      mg->solver[ilevel] = solver3;
      mg->A[ilevel] = A3;
      mg->P[ilevel] = matP3;
      mg->x[ilevel] = phgVecCopy(solver3->rhs, NULL);
      mg->f[ilevel] = phgVecCopy(solver3->rhs, NULL);
      mg->r[ilevel] = phgVecCopy(solver3->rhs, NULL);
      phgOptionsPush();
      phgOptionsSetOptions(S3_opts);
      mg->smoother[ilevel] = phgMat2Solver(SOLVER_DEFAULT, mg->A[ilevel]);
      mg->smoother[ilevel]->warn_maxit = FALSE;
      mg->smoother[ilevel]->verb = 0;
      phgOptionsPop();
    }
  }

  if (nlevel >= 4) {

    ilevel = nlevel - 4;
    mg->fext[ilevel] = phgMapCreateVec(matP3->cmap, 1); /* ext vec */
    if (g4) {
      mg->dof[ilevel] = u4;
      mg->solver[ilevel] = solver4;
      mg->A[ilevel] = A4;
      mg->P[ilevel] = NULL;	/* note */
      mg->x[ilevel] = phgVecCopy(solver4->rhs, NULL);
      mg->f[ilevel] = phgVecCopy(solver4->rhs, NULL);
      mg->r[ilevel] = phgVecCopy(solver4->rhs, NULL);
      phgOptionsPush();
      phgOptionsSetOptions(S4_opts);
      mg->smoother[ilevel] = phgMat2Solver(SOLVER_DEFAULT, mg->A[ilevel]);
      mg->smoother[ilevel]->warn_maxit = FALSE;
      mg->smoother[ilevel]->verb = 0;
      phgOptionsPop();
    }
  }


  /* MG info */
  phgPrintf("\n\n-----\nMG info\n");;
  for (ilevel = nlevel - 1; ilevel >= 0; ilevel--) {
    if (mg->dof[ilevel])
      phgPrintf("MG # Dof[%d]: %10d (%10d)\n", ilevel, 
		DofGetDataCountGlobal(mg->dof[ilevel]),
		DofGetDataCount(mg->dof[ilevel])
		);
  }
  phgPrintf("-----\nEnd MG info\n\n");;

  
    
  return mg;
}




void 
destroy_mg(MG_DATA **mg_ptr)
{
    MG_DATA *mg = *mg_ptr;

    /* Free MG */
    phgVecDestroy(&mg->x[1]);
    phgVecDestroy(&mg->f[1]);
    phgVecDestroy(&mg->r[1]);
    phgSolverDestroy(&mg->smoother[1]);

    if (mg->x[0]) {
	phgVecDestroy(&mg->x[0]);
	phgVecDestroy(&mg->f[0]);
	phgVecDestroy(&mg->r[0]);
	phgSolverDestroy(&mg->smoother[0]);

	phgSolverDestroy(&mg->solver[0]);
	phgDofFree(&mg->dof[0]);
    }	

    /* Free vectors */
    phgMatDestroy(&mg->P[1]);
    phgVecDestroy(&mg->fext[0]);

    phgFree(mg);
    mg_ptr = NULL;
}
	
