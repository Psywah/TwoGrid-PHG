/*
 *
 *
 * Test two grid MG.
 *
 *
 *  (highlight-regexp "[a-zA-Z]+1" 'hi-yellow)
 *  (highlight-regexp "[a-zA-Z]+2" 'hi-green)
 *  */
#include <time.h>
#include "phg.h"
#include "two-grid.h"
#include "pgrid.h"

/* Segv */
#include <signal.h>
#include <execinfo.h> // backtrace(), etc
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>

static FLOAT a = 1.0;

#define TEST_ORDER 2
static void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
#if TEST_ORDER == 1    
    *v = x + 2*y - 3*z;
#else
    *v = x*x + 2*y*y - 3*z*z;
#endif
}    

static void
func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    func_u(x, y, z, value);
    *value *= a;
    
#if TEST_ORDER == 1
    *value += 0;
#else
    *value += - (2 + 4 - 6);
#endif
}





static void
build_linear_system(SOLVER *solver, DOF *u_h)
{
    int i, j;
    GRID *g = u_h->g;
    SIMPLEX *e;

    assert(u_h->dim == 1);
    ForAllElements(g, e) {
	int N = DofGetNBas(u_h, e);	/* number of basises in the element */
	FLOAT A[N][N], rhs[N], buffer[N];
	INT Iu[N];

	/* compute \int \grad\phi_i \cdot \grad\phi_j making use of symmetry */
	for (i = 0; i < N; i++) {
	    Iu[i] = phgSolverMapE2L(solver, 0, e, i);
	    for (j = 0; j <= i; j++) {
		A[j][i] = A[i][j] =
		    /* stiffness */
		    phgQuadGradBasDotGradBas(e, u_h, i, u_h, j, QUAD_DEFAULT) +
		    /* mass */
		    a * phgQuadBasDotBas(e, u_h, i, u_h, j, QUAD_DEFAULT);
	    }
	}

	/* loop on basis functions */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(u_h, e, i, func_u, buffer, rhs+i,
					DOF_PROJ_NONE)) {
		phgSolverAddMatrixEntries(solver, 1, Iu + i, N, Iu, buffer); 
	    }
	    else {	/* interior node */
		/* right hand side = \int f * phi_i */
		rhs[i] = phgQuadFuncDotBas(e, func_f, u_h, i, QUAD_DEFAULT);
		phgSolverAddMatrixEntries(solver, 1, Iu + i, N, Iu, A[i]); 
	    }
	}
	phgSolverAddRHSEntries(solver, N, Iu, rhs);
    }
}




/* BOOLEAN */
/* part_no_action(GRID *g, int nprocs, DOF *weights, FLOAT power); */

/* BOOLEAN */
/* part_no_action(GRID *g, int nprocs, DOF *weights, FLOAT power) */
/* /\* */
/*  * Partioner with no action, mark is set some time else. */
/*  * */
/*  *  *\/ */
/* { */
/*     phgPrintf("part no action\n"); */
/*     return TRUE; */
/* } */



static
void handler(int sig) 
{
  void *array[100];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 100);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);

  MPI_Finalize();             /* Cause other procs waiting,
			       * then manually interupte mpirun.
			       *  */

  fprintf(stderr, "MPI finalized.\n");
  exit(1);
}

void install_error_log()
{
  /* Redirect stderr */
  {
    char name[1000];
    FILE *fp;

    sprintf(name, "err/%04d.log", phgRank);
    fp = fopen(name, "w");
    assert(fp != NULL);
    dup2(fileno(fp), fileno(stderr));
    fclose(fp);
    fprintf(stderr, "redirect stderr.\n");
  }


  // install segv handler
  signal(SIGSEGV, handler);
  signal(SIGABRT, handler);
  signal(SIGINT, handler);
  signal(SIGFPE, handler);
  signal(SIGTERM, handler);
}

//#define partitioner phgPartitionParMetis
#define partitioner phgPartitionSFC

int
main(int argc, char *argv[])
{
    GRID *g;
    SIMPLEX *e;
    static const char *fn = NULL; 
    static const char *fn2 = NULL; 
    static const char *fn3 = NULL; 
    static const char *fn4 = NULL; 
    int pre_refines = 0;
    int nLevel = 2, ilevel;
    DOF *u;
    /* MAP *map;  */
    /* VEC *vec; */
    BOOLEAN debug = FALSE;
    BOOLEAN export = FALSE;
    int Cmat_type = 0;
    int nprocs2 = 1;
    int nprocs3 = 1;
    int nprocs4 = 1;
    char *S_opts = NULL;
    char *S2_opts = NULL;
    char *S3_opts = NULL;
    char *S4_opts = NULL;
    char hostname[256];
    FLOAT time0;

    install_error_log();

    phgOptionsRegisterNoArg("debug", "mpi debug", &debug);
    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&fn);
    phgOptionsRegisterFilename("mesh_file2", "Mesh file2", (char **)&fn2);
    phgOptionsRegisterFilename("mesh_file3", "Mesh file3", (char **)&fn3);
    phgOptionsRegisterFilename("mesh_file4", "Mesh file4", (char **)&fn4);
    phgOptionsRegisterInt("pre_refines", "Pre-refines", &pre_refines);
    phgOptionsRegisterInt("Cmat_type", "Coarse mat type", &Cmat_type);

    phgOptionsRegisterString("S_opts", "smoother options", &S_opts);
    phgOptionsRegisterString("S2_opts", "smoother options", &S2_opts);
    phgOptionsRegisterString("S3_opts", "smoother options", &S3_opts);
    phgOptionsRegisterString("S4_opts", "smoother options", &S4_opts);

    phgOptionsRegisterInt("nprocs2", "Coarse nprocs", &nprocs2);
    phgOptionsRegisterInt("nprocs3", "Coarse nprocs", &nprocs3);
    phgOptionsRegisterInt("nprocs4", "Coarse nprocs", &nprocs4);

    phgOptionsRegisterInt("nlevel", "# MG level", &nLevel);


    phgInit(&argc, &argv);


    if (debug) {
        int _i_ = 0;
        unsigned int t = 5;
        int pid = getpid();

        gethostname(hostname, sizeof(hostname));
        printf("#### Lengweee debug "
               "PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        while (0 == _i_) {
            MPI_Barrier(MPI_COMM_WORLD);
            printf("%d after bar\n", pid);
            fflush(stdout);
            for (t=0; t<50; t++)
                sleep(1);
            printf("%d after sleep\n", pid);
            fflush(stdout);
        }
        printf("### PID %d, ready!\n", getpid());
    }


    g = phgNewGrid(-1);
    if (!phgImport(g, fn, FALSE))
	phgError(1, "can't read file \"%s\".\n", fn);


#if 0    
    /* distribute */
    if (phgBalanceGrid(g, 1.2, 1, NULL, 0.))
	phgPrintf("Repartition mesh, load imbalance: %lg\n",
		  (double)g->lif);
#else
    /* distribute, with refine  */
    {
	//phgPartUserSetFunc(part_no_action);

	time0 = phgGetTime(NULL);
	partitioner(g, MPI_COMM_WORLD, NULL, 0.);
	phgRedistributeGrid(g);
	phgPrintf("* Redist time %d: %lfs\n", 
		  phgNProcs, phgGetTime(NULL) - time0);

	/* pre-refine */
	phgRefineAllElements(g, pre_refines);
	phgPrintf("pre_refines: %d\n", pre_refines);
	phgPrintf("# of elem: %d\n", g->nelem_global);

	time0 = phgGetTime(NULL);
	partitioner(g, MPI_COMM_WORLD, NULL, 0.);
	phgRedistributeGrid(g);
	phgPrintf("* Redist time %d: %lfs\n", 
		  phgNProcs, phgGetTime(NULL) - time0);
    }
#endif
    phgInfo(0, "_nlocal: %d, nglobal: %d\n", g->nvert, g->nvert_global);

    if (export) phgExportVTK(g, "g1.vtk", NULL);

    /* DOF */
    u = phgDofNew(g, DOF_DEFAULT, 1, "u_h", func_u);




#if 0
    /*
     * Serial
     * Locate points test for cube.8x8x8.mesh [-1,1]^3
     *
     *  */
    {
	FLOAT x[3];
	int i, k;
	BOOLEAN found = FALSE;
	SIMPLEX *e;
	FLOAT lambda[Dim+1], v0, v1;

	og = build_octree(g);

	srand (time(NULL));


	for (i = 0; i < 100; i++) { 
	    /* Test many points */

	    for (k = 0; k < Dim; k++) {
		x[k] = 2 * (rand() / (FLOAT) RAND_MAX) - 1.;
	    }

	    found = locate_point(og, x, &e, lambda);

	    phgInfo(0, "locate (%12.8e, %12.8e, %12.8e) found %d\n", 
		    x[0], x[1], x[2], found);
	    if (found) {
		phgInfo(0, "found in elem %d\n", e->index);

		/* value check */
	    
		func_u(x[0], x[1], x[2], &v0);
		phgDofEval(u, e, lambda, &v1);
		phgInfo(0, "eval error: %e\n", fabs(v0 - v1));
	    }

	    if (phgNProcs == 1) {
		assert(found);
	    }
	}
	
    }
#endif



    /*
     * Parallel
     * Locate points test for cube.8x8x8.mesh [-1,1]^3
     *
     *  */
    {
    }




    /*
     * 
     * Two grid MG
     *
     *  */
    {
	GRID *g2 = NULL;
	GRID *g3 = NULL;
	GRID *g4 = NULL;
	GRID *g2_ = NULL;
	GRID *g3_ = NULL;
	GRID *g4_ = NULL;
	SIMPLEX *e2 = NULL;
	SIMPLEX *e3 = NULL;
	SIMPLEX *e4 = NULL;

	g2 = phgNewGrid(-1);
	if (!phgImport(g2, fn2, FALSE))
	    phgError(1, "can't read file \"%s\".\n", fn2);
	if (nLevel >= 3) {
	  g3 = phgNewGrid(-1);
	  if (!phgImport(g3, fn3, FALSE))
	    phgError(1, "can't read file \"%s\".\n", fn3);
	}
	if (nLevel >= 4) {
	  g4 = phgNewGrid(-1);
	  if (!phgImport(g4, fn4, FALSE))
	    phgError(1, "can't read file \"%s\".\n", fn4);
	}

#if 0	
	/* distribute on comm_world */
	if (phgBalanceGrid(g2, 1.2, 1, NULL, 0.))
	    phgPrintf("Repartition mesh, load imbalance: %lg\n",
		      (double)g2->lif);
	if (export) phgExportVTK(g2, "g2.vtk", NULL);
#else
	/* distribute on nprocs */
	MPI_Group MPI_COMM_GROUP, group;
	MPI_Comm comm2;
	MPI_Comm comm3;
	MPI_Comm comm4;
	int i, *ranks;

	{
	  /* Create sub comm */
	  ranks = phgCalloc(nprocs2, sizeof(*ranks));
	  for (i = 0; i < nprocs2; i++)
	    ranks[i] = i;
	
	  MPI_Comm_group(MPI_COMM_WORLD, &MPI_COMM_GROUP);
	  MPI_Group_incl(MPI_COMM_GROUP, nprocs2, ranks, &group);
	  MPI_Comm_create(MPI_COMM_WORLD, group, &comm2);
	  phgFree(ranks);
	  /* TODO!!! */
	  /* free sub group and sub comm, somewhere */


	  {
	    //phgPartUserSetFunc(part_no_action);
	    time0 = phgGetTime(NULL);

	    if (phgRank < nprocs2)
	      partitioner(g2, comm2, NULL, 0.);
	    phgRedistributeGrid(g2);
	    phgPrintf("* Redist time %d: %lfs\n", 
		      phgNProcs, phgGetTime(NULL) - time0);
	  }

	  if (export) phgExportVTK(g2, "g2_first.vtk", NULL);

	    
	  {
	    phgInfo(0, "g2 nelem %d %d\n", phgRank, g2->nelem);
	    if (phgRank >= nprocs2)
	      assert(g2->nelem == 0);

	    if (phgRank < nprocs2) {
	      /* Re import */
	      TET *tet, *t;
	      SIMPLEX *e2;
	      int i;

	      tet = phgCalloc(g2->nelem, sizeof(*tet));
	      t = tet;
	      ForAllElements(g2, e2) {
		for (i = 0; i < NVert; i++)
		  t->verts[i] = e2->verts[i];
		for (i = 0; i < NFace; i++) {
		  t->bound_type[i] = e2->bound_type[i];
		  t->bound_type[i] &= ~OWNER;
		}
		t->e = e2;
		t++;
	      }

	      g2_ = phgImportParallelGrid_(NULL,
					   g2->nvert,
					   g2->nelem,
					   g2->nvert_global,
					   g2->nelem_global,
					   g2->L2Gmap_vert,
					   g2->L2Gmap_elem,
					   g2->verts[0],
					   tet,
					   comm2
					   );
	      phgFree(tet);
	    }
	    else {
	      g2_ = NULL;
	    }

	    phgFreeGrid(&g2);
	    g2 = g2_;
	  }
	  if (g2)
	    phgInfo(0, "_nlocal: %d, nglobal: %d\n", g2->nvert, g2->nvert_global);

	  phgInfo(0, "g2: %x\n", g2);
	  if (g2) 
	    if (export) phgExportVTK(g2, "g2.vtk", NULL);
	}


	if (nLevel >= 3) {
	  /* Create sub comm */
	  ranks = phgCalloc(nprocs3, sizeof(*ranks));
	  for (i = 0; i < nprocs3; i++)
	    ranks[i] = i;
	
	  MPI_Comm_group(MPI_COMM_WORLD, &MPI_COMM_GROUP);
	  MPI_Group_incl(MPI_COMM_GROUP, nprocs3, ranks, &group);
	  MPI_Comm_create(MPI_COMM_WORLD, group, &comm3);
	  phgFree(ranks);
	  /* TODO!!! */
	  /* free sub group and sub comm, somewhere */


	  {
	    //phgPartUserSetFunc(part_no_action);
	    time0 = phgGetTime(NULL);

	    if (phgRank < nprocs3)
	      partitioner(g3, comm3, NULL, 0.);
	    phgRedistributeGrid(g3);
	    phgPrintf("* Redist time %d: %lfs\n", 
		      phgNProcs, phgGetTime(NULL) - time0);
	  }

	  if (export) phgExportVTK(g3, "g3_first.vtk", NULL);

	    
	  {
	    phgInfo(0, "g3 nelem %d %d\n", phgRank, g3->nelem);
	    if (phgRank >= nprocs3)
	      assert(g3->nelem == 0);

	    if (phgRank < nprocs3) {
	      /* Re import */
	      TET *tet, *t;
	      SIMPLEX *e3;
	      int i;

	      tet = phgCalloc(g3->nelem, sizeof(*tet));
	      t = tet;
	      ForAllElements(g3, e3) {
		for (i = 0; i < NVert; i++)
		  t->verts[i] = e3->verts[i];
		for (i = 0; i < NFace; i++) {
		  t->bound_type[i] = e3->bound_type[i];
		  t->bound_type[i] &= ~OWNER;
		}
		t->e = e3;
		t++;
	      }

	      g3_ = phgImportParallelGrid_(NULL,
					   g3->nvert,
					   g3->nelem,
					   g3->nvert_global,
					   g3->nelem_global,
					   g3->L2Gmap_vert,
					   g3->L2Gmap_elem,
					   g3->verts[0],
					   tet,
					   comm3
					   );
	      phgFree(tet);
	    }
	    else {
	      g3_ = NULL;
	    }

	    phgFreeGrid(&g3);
	    g3 = g3_;
	  }
	  if (g3)
	    phgInfo(0, "_nlocal: %d, nglobal: %d\n", g3->nvert, g3->nvert_global);

	  phgInfo(0, "g3: %x\n", g3);
	  if (g3) 
	    if (export) phgExportVTK(g3, "g3.vtk", NULL);
	}



	if (nLevel >= 4) {
	  /* Create sub comm */
	  ranks = phgCalloc(nprocs4, sizeof(*ranks));
	  for (i = 0; i < nprocs4; i++)
	    ranks[i] = i;
	
	  MPI_Comm_group(MPI_COMM_WORLD, &MPI_COMM_GROUP);
	  MPI_Group_incl(MPI_COMM_GROUP, nprocs4, ranks, &group);
	  MPI_Comm_create(MPI_COMM_WORLD, group, &comm4);
	  phgFree(ranks);
	  /* TODO!!! */
	  /* free sub group and sub comm, somewhere */


	  {
	    //phgPartUserSetFunc(part_no_action);
	    time0 = phgGetTime(NULL);

	    if (phgRank < nprocs4)
	      partitioner(g4, comm4, NULL, 0.);
	    phgRedistributeGrid(g4);
	    phgPrintf("* Redist time %d: %lfs\n", 
		      phgNProcs, phgGetTime(NULL) - time0);
	  }

	  if (export) phgExportVTK(g4, "g4_first.vtk", NULL);

	    
	  {
	    phgInfo(0, "g4 nelem %d %d\n", phgRank, g4->nelem);
	    if (phgRank >= nprocs4)
	      assert(g4->nelem == 0);

	    if (phgRank < nprocs4) {
	      /* Re import */
	      TET *tet, *t;
	      SIMPLEX *e4;
	      int i;

	      tet = phgCalloc(g4->nelem, sizeof(*tet));
	      t = tet;
	      ForAllElements(g4, e4) {
		for (i = 0; i < NVert; i++)
		  t->verts[i] = e4->verts[i];
		for (i = 0; i < NFace; i++) {
		  t->bound_type[i] = e4->bound_type[i];
		  t->bound_type[i] &= ~OWNER;
		}
		t->e = e4;
		t++;
	      }

	      g4_ = phgImportParallelGrid_(NULL,
					   g4->nvert,
					   g4->nelem,
					   g4->nvert_global,
					   g4->nelem_global,
					   g4->L2Gmap_vert,
					   g4->L2Gmap_elem,
					   g4->verts[0],
					   tet,
					   comm4
					   );
	      phgFree(tet);
	    }
	    else {
	      g4_ = NULL;
	    }

	    phgFreeGrid(&g4);
	    g4 = g4_;
	  }
	  if (g4)
	    phgInfo(0, "_nlocal: %d, nglobal: %d\n", g4->nvert, g4->nvert_global);

	  phgInfo(0, "g4: %x\n", g4);
	  if (g4) 
	    if (export) phgExportVTK(g4, "g4.vtk", NULL);
	}


#endif




	/* ----------------------------------------
	 * 
	 * Define poisson problem
	 *
	 * ---------------------------------------- */
	SOLVER *solver;
	solver = phgSolverCreate(SOLVER_DEFAULT, u, NULL);

	time0 = phgGetTime(NULL);
	build_linear_system(solver, u);
	phgPrintf("Build mat time: %lfs\n", phgGetTime(NULL) - time0);

	



	
	/* ------------
	 *
	 * Mg create  
	 *
	 * ------------ */
#if 0
	MG_DATA *mg = mg_create2(g, g2, u,
				solver,
				S_opts, S2_opts,
				Cmat_type, build_linear_system);
#else
	MG_DATA *mg = mg_createN(nLevel, g, g2, g3, g4,
				 u, solver,
				 S_opts, S2_opts, S3_opts, S4_opts,
				 Cmat_type, build_linear_system);
#endif


	
	if (0) {
	  /* ----------------------------------------
	   * 
	   * two grid Mg solve 
	   *
	   * ---------------------------------------- */
	  mg_solve(mg, solver, u);
	}
	else {

	  /* ----------------------------------------
	   * 
	   * two grid Mg as PC
	   *
	   * ---------------------------------------- */
	  SOLVER *pc = NULL;
	  pc = phgMat2Solver(SOLVER_GMRES, solver->mat);
	  phgSolverSetPC(solver, pc, mg_pc_proc);
	  pc->mat->mv_data = phgAlloc(sizeof(*pc->mat->mv_data));
	  pc->mat->mv_data[0] = (void *) mg;

	
	  time0 = phgGetTime(NULL);
	  phgSolverSolve(solver, FALSE, u, NULL);
	  phgPrintf("Solving time: %lfs\n", phgGetTime(NULL) - time0);

	  phgSolverDestroy(&solver);
	  phgSolverDestroy(&pc);
	}


	/* free mg */
	destroy_mg(&mg);

	
	phgFreeGrid(&g2);
    }

    



    phgDofFree(&u);
    phgFreeGrid(&g);
    MPI_Barrier(MPI_COMM_WORLD);
    phgPrintf("\n\n\n---\nDone.\n");
    phgFinalize();

    return 0;
}












