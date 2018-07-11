#ifndef PARALELL_GRID_
#define PARALELL_GRID_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct TET_ {
    int verts[4];
    //int neighbours[4];
    BTYPE bound_type[4];
    SIMPLEX *e;
    SHORT region_mark;
} TET;


typedef BOOLEAN (*SUBREGION_FUNC)(GRID *g, SIMPLEX *e);

GRID *phgImportParallelGrid_(GRID *old_grid,
			    int nvert, int nelem, int nvert_global, int nelem_global,
			    INT *L2Gmap_vert, INT *L2Gmap_elem, double *coord, TET *tet,
			    MPI_Comm comm);

GRID *phgGetSubGrid_(GRID *g, SUBREGION_FUNC subregion_func);

void poisson_test(GRID *g);



#ifdef __cplusplus
}
#endif
#endif
