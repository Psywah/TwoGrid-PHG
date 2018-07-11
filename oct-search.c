//#include "vtk-draw.h"
#include "oct-search.h"

/*
 *
 *
 * oct tree implementation
 *
 *
 * */



#if DEBUG_SEARCH
static int verb = 2;
static char LEVEL[8][100] = {"   ",
			     "      ",
			     "         ",
			     "            ",
			     "               ", 
			     "                  ", 
			     "                     "};
#endif


/* Search grid varibles */
static COORD *tet_center;
static FLOAT *tet_inner_diam;
static GRID *g_;		/* current grid */
static _CUBE **cube_list;
static int ncube_list;  	/* # of cube in lists */
static int ncube_size;

/* static functions */
static BOOLEAN test_intersection(SIMPLEX * tet, _CUBE * cube);


#define DRAW_TET_OF__CUBE(e) DRAW_TET_OF__CUBE_(e, 1) 

#define DRAW_TET_OF__CUBE_(e, with_box) {	\
	int j_;					\
	vtkTmpActorsBegin(verb);		\
	if (with_box) {				\
	    vtkSetColor(verb, "red");		\
	    vtkSetTransparent(verb, 1.);	\
	    vtkDrawBox(verb, e->bbox[0]);	\
	}					\
						\
	for (j_ = 0; j_ < e->ntet; j_++) {	\
	    vtkSetColor(verb, "white");		\
	    vtkSetTransparent(verb, .25);	\
	    vtkDrawElement(verb, e->liste[j_]);	\
	}					\
	vtkPause(verb);				\
	vtkTmpActorsClear(verb);		\
    }


#define DRAW_TETS_OF__CUBE_(e, with_box) {	\
	int j_;					\
	vtkTmpActorsBegin(verb);		\
	if (with_box) {				\
	    vtkSetColor(verb, "red");		\
	    vtkSetTransparent(verb, .25);	\
	    vtkDrawBox(verb, e->bbox[0]);	\
	}					\
						\
	for (j_ = 0; j_ < e->ntet; j_++) {	\
	    vtkTmpActorsBegin(verb);		\
	    vtkSetColor(verb, "white");		\
	    vtkSetTransparent(verb, .25);	\
	    vtkDrawElement(verb, e->liste[j_]);	\
						\
	    vtkPause(verb);			\
	    vtkTmpActorsClear(verb);		\
	}					\
	vtkPause(verb);				\
	vtkTmpActorsClear(verb);		\
    }






void _octtree_insert(_CUBE * cube_father)
/* insert children to cube */
{
    int i, j, k, l;
    _CUBE *new_cubes, *cube_child;
    const FLOAT *midx = cube_father->midx;
    int verb = 3;

    new_cubes = phgCalloc(8, sizeof(*new_cubes));

    /* record */
    if (ncube_list >= ncube_size) {
	cube_list =
	    phgRealloc_(cube_list,
			2 * ncube_size * sizeof(_CUBE),
			ncube_size * sizeof(_CUBE));
	ncube_size *= 2;
    }
    cube_list[ncube_list++] = new_cubes;
    
    
    if (verb <= 2)
	DUMP_BOX(cube_father->bbox);
    phgInfo(verb, "midx [%f %f %f]\n", midx[0], midx[1], midx[2]);

    for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
	    for (k = 0; k < 2; k++) {
		//cube_child = calloc(1, sizeof(*cube_child));
		cube_child = new_cubes + 4 * i + 2 * j + k;
		cube_father->children[4 * i + 2 * j + k] = cube_child;
		cube_child->father = cube_father;

		cube_child->level = cube_father->level + 1;

		memcpy(&cube_child->bbox, &cube_father->bbox,
		       sizeof(cube_child->bbox));
		if (i == 0)  cube_child->bbox[0][1] = midx[0];
		if (i == 1)  cube_child->bbox[0][0] = midx[0];
		if (j == 0)  cube_child->bbox[1][1] = midx[1];
		if (j == 1)  cube_child->bbox[1][0] = midx[1];
		if (k == 0)  cube_child->bbox[2][1] = midx[2];
		if (k == 1)  cube_child->bbox[2][0] = midx[2];
		phgInfo(verb, " [%d %d %d]\n", i, j, k);
		if (verb <= 2)
		    DUMP_BOX(cube_child->bbox);

		for (l = 0; l < Dim; l++)
		    cube_child->midx[l] =
			.5 * (cube_child->bbox[l][0] +
			      cube_child->bbox[l][1]);

		FLOAT dx = cube_child->bbox[0][1] - cube_child->bbox[0][0];
		FLOAT dy = cube_child->bbox[1][1] - cube_child->bbox[1][0];
		FLOAT dz = cube_child->bbox[2][1] - cube_child->bbox[2][0];
		cube_child->diam = sqrt(dx * dx + dy * dy + dz * dz);
		cube_child->dh = cube_father->dh / 2.;

		cube_child->size = 10;
		cube_child->liste =
		    phgCalloc(cube_child->size,
			      sizeof(*cube_child->liste));
	    }

}










/* --------------------------------------------------------------------------------
 *            
 *            Traverse subroutines
 *
 * -------------------------------------------------------------------------------- */
#undef CB_ARGS
#undef CB_ARGS0
#if 0
/* carry g in the argument list of callback functions */
#define CB_ARGS(e)	(_OCT_TREE *og, _CUBE *e)
#define CB_ARGS0(g, e)	(og, e)
#else
/* don't carry g in the argument list of callback functions */
#define CB_ARGS(e)	(_CUBE *e)
#define CB_ARGS0(g, e)	(e)
#endif

static _OCT_TREE *og_traverse;	/* used by phgTraverseElements, other functions
				   shouldn't use it */

/* flag controlling whether to call callback function for non-leaf elements */
static BOOLEAN traverse_all_flag = FALSE;

/* Depth of the element in the tree.
   This variable is made global and may be used by the callback functions */
INT _octTreeTraverseDepth = 0, _octTreeTraverseIndex = 0;

static BOOLEAN traverse_element(_CUBE * e, BOOLEAN(*cb) CB_ARGS(e))
/*
 *  Stop traversing if returns FALSE. 
 *  If TRUE all the time, traverse the whole tree.
 *  */
{
    BOOLEAN is_leaf = TRUE;
    int ichild = 0;

    for (ichild = 0; ichild < 8; ichild++) {
	if (e->children[ichild] != NULL) {
	    is_leaf = FALSE;
	    _octTreeTraverseDepth++;
	    if (!traverse_element(e->children[ichild], cb))
		return FALSE;
	    _octTreeTraverseDepth--;
	}
    }

    if (!is_leaf && !traverse_all_flag) {
	return TRUE;
    } else {
	BOOLEAN ret = cb CB_ARGS0(og_traverse, e);
	Unused(og_traverse);
	_octTreeTraverseIndex++;
	return ret;
    }
}

void _octTreeTraverseElements(_OCT_TREE * og, BOOLEAN(*cb) CB_ARGS(e))
/* function traversing (traverse_all_flag ? leaf : all) elements.
 * Stop traversing if callback function returns FALSE */
{
    INT i;

    if (og->nroot == 0)
	return;

    og_traverse = og;
    _octTreeTraverseIndex = 0;
    for (i = 0; i < og->nroot; i++) {
	_octTreeTraverseDepth = 0;
	if (!traverse_element(og->roots + i, cb))
	    break;
    }
}

void _octTreeTraverseAllElements(_OCT_TREE * og, BOOLEAN(*cb) CB_ARGS(e))
/* function traversing all (including non leaf) elements */
{
    traverse_all_flag = TRUE;
    _octTreeTraverseElements(og, cb);
    traverse_all_flag = FALSE;
}



static void divide_herds(_CUBE * cube_father)
/* divide the element list to childrens */
{
    SIMPLEX *tet;
    _CUBE *cube;
    int ichild, i;


    for (ichild = 0; ichild < 8; ichild++) {
	cube = cube_father->children[ichild];

	for (i = 0; i < cube_father->ntet; i++) {
	    tet = cube_father->liste[i];

	    if (test_intersection(tet, cube)) {

		/* found */
		if (cube->ntet >= cube->size) {
		    cube->liste =
			phgRealloc_(cube->liste,
				    2 * cube->size * sizeof(SIMPLEX *),
				    cube->size * sizeof(SIMPLEX *));
		    cube->size *= 2;
		}

		cube->liste[cube->ntet++] = tet;
	    }
	    /* end of intersect */
	}

	phgInfo(3, "  cube level: %d, ntet: %d\n", cube->level,
		cube->ntet);
    }

}





static BOOLEAN 
traverse_insert(_CUBE * e, BOOLEAN(*cb) CB_ARGS(e))
/*
 * 1) if cube has small # of tets, no need to divide
 * 2) if cube is geomatrically small enough, no need to divide
 * 3) if refine level is too deep
 *
 * 
 * Remark: Bad divide condition:
 *    if cube has many tets, keep divide.
 *    A cube could have 100+ tets, and one of its children still have
 *    the same # of tets, which makes divide goes to endless loop.
 * 
 *  */
{
    int ichild = 0;
    int i;

    if (e->ntet < MAX_NELEM)
	return TRUE;

    /* cube size small enough */
    for (i = 0; i < e->ntet; i++) {
#if 0
	/* Use inner sphr diam */
	if (e->dh > 2 * tet_inner_diam[e->liste[i]->index]) {
	    break;
	}
#else
	/* Use outer sphr diam */
	if (e->dh > .25 * phgGeomGetDiameter(g_, e->liste[i])) {
	    break;
	}
#endif
    }
    if (i >= e->ntet) {
	int  verb = 3;
	//phgWarning("Cube too small.\n");
	DRAW_TET_OF__CUBE(e);

	return TRUE;
    }


    /* reach max # level */
    if (_octTreeTraverseDepth == MAX_NLEVEL) {
	int  verb = 3;
	phgWarning("Too many levels\n");

	DRAW_TET_OF__CUBE(e);
	return TRUE;
    }

#if 0
    /* Refine cube does not help.
     *
     * Remark: obselete.
     *   Use # of tets to test divide is bad.
     * */
    if (e->father != NULL && e->ntet < MAX_NELEM0	/* above MAX_NELEM0(60),
							 * keep dividing */
	&& e->ntet == e->father->ntet) {  /* divide not improving */
	int  verb = 3;
	phgWarning("Refine not helping\n");

	DRAW_TET_OF__CUBE(e);
	return TRUE;
    }
#endif


    _octtree_insert(e);
    divide_herds(e);

    phgInfo(3, "Insert cube\n");

    for (ichild = 0; ichild < 8; ichild++) {
	_octTreeTraverseDepth++;
	if (!traverse_insert(e->children[ichild], cb))
	    return FALSE;
	_octTreeTraverseDepth--;
    }

    return TRUE;
}





/* const static FLOAT * */
/* geomXYZ2Lambda(const FLOAT *J, FLOAT x, FLOAT y, FLOAT z) */
/* { */
/*     static FLOAT lambda[Dim + 1]; */

/*     lambda[0] = J[0] * x + J[1] * y + J[2] * z + J[3]; */
/*     J += 4; */
/*     lambda[1] = J[0] * x + J[1] * y + J[2] * z + J[3]; */
/*     J += 4; */
/*     lambda[2] = J[0] * x + J[1] * y + J[2] * z + J[3]; */
/*     J += 4; */
/*     lambda[3] = J[0] * x + J[1] * y + J[2] * z + J[3]; */

/*     return lambda; */
/* } */





static BOOLEAN 
test_intersection(SIMPLEX * tet, _CUBE * cube)
/* Test tetra cube intersection */
{
    int i, j, k, v, l, f;;
    int verb = 3;

    /*
     * Step 1. distance test
     * */
    FLOAT dir[Dim], dist;
    for (k = 0; k < Dim; k++)
	dir[k] = //tet->center[k]
	    tet_center[tet->index][k]
	    - cube->midx[k];
    dist = sqrt(INNER_PRODUCT(dir, dir));
    if (dist > 1.5 * phgGeomGetDiameter(g_, tet)
	     + 1.5 * cube->diam + EPS) {
	phgInfo(3, "   Too far to intersect.\n");
	return FALSE;
    }

    /*
     * Step 2. tet vert in cube test
     * */
    for (v = 0; v < NVert; v++) {
	FLOAT *xe = g_->verts[tet->verts[v]];
	if (cube->bbox[0][0] - EPS < xe[0]
	    && xe[0] < cube->bbox[0][1] + EPS
	    && cube->bbox[1][0] - EPS < xe[1]
	    && xe[1] < cube->bbox[1][1] + EPS
	    && cube->bbox[2][0] - EPS < xe[2]
	    && xe[2] < cube->bbox[2][1] + EPS) {
	    phgInfo(3, "   Tet vert in cube.\n");
	    return TRUE;
	}
    }

    /*
     * Step 3. cube vert in tet test
     * */
    for (i = 0; i < 2; i++)
	for (j = 0; j < 2; j++)
	    for (k = 0; k < 2; k++) {
		FLOAT xc[3] = { cube->bbox[0][i],
				cube->bbox[1][j],
				cube->bbox[2][k]
		};
		const FLOAT *lambda =
		    //geomXYZ2Lambda(tet->J[0], xc[0], xc[1], xc[2]);
		    phgGeomXYZ2Lambda(g_, tet, xc[0], xc[1], xc[2]);
		if (LambdaInElement(lambda)) {
		    phgInfo(3, "   Cube vert in Tet.\n");
		    return TRUE;
		}
	    }


    /*
     * Step 4. tet seg & cube face test
     *
     * */
    for (l = 0; l < NEdge; l++) {
	FLOAT *x0 = g_->verts[tet->verts[GetEdgeVertex(l, 0)]];
	FLOAT *x1 = g_->verts[tet->verts[GetEdgeVertex(l, 1)]];

	int xy[3][2] = { {1, 2}, {0, 2}, {0, 1} };
	FLOAT t, xt[3];		/* t is the parameters to discribe line:
				 * x = x0 + t * (x1 - x0)  */

	//xy[k][0];
	for (k = 0; k < 3; k++) {	/* face_{xyz} */
	    for (f = 0; f < 2; f++) {	/* face +1,-1 */

		if (fabs(x1[k] - x0[k]) < 1e-12) {
		    /* seg is parallel to face_{xyz} */
		    continue;
		}

		t = (cube->bbox[k][f] - x0[k]) / (x1[k] - x0[k]);

		if (t > 1 + EPS0 || t < -EPS0) {
		    /* out of segment */
		    continue;
		}

		xt[0] = x0[0] + (x1[0] - x0[0]) * t;
		xt[1] = x0[1] + (x1[1] - x0[1]) * t;
		xt[2] = x0[2] + (x1[2] - x0[2]) * t;

		/* test xt in [xy] */
		int d0 = xy[k][0], d1 = xy[k][1];
		if (cube->bbox[d0][0] - EPS < xt[d0]
		    && xt[d0] < cube->bbox[d0][1] + EPS
		    && cube->bbox[d1][0] - EPS < xt[d1]
		    && xt[d1] < cube->bbox[d1][1] + EPS) {

		    phgInfo(3, "inter pts: %f %f %f, %e\n",
			    xt[0], xt[1], xt[2], t);
		    //DUMP_BOX(cube->bbox);

		    phgInfo(3, "   Tet seg intersect cube face.\n");
		    return TRUE;
		}

	    }			/* end face +1,-1 */
	}			/* end face_{xyz} */
    }				/* end elem edge */

    /*
     * Step 5. cube seg & tet face test
     *
     * */
    static int cube_edge[12][2] = {
	{0, 1}, {1, 2}, {2, 3}, {3, 0},
	{4, 5}, {5, 6}, {6, 7}, {7, 4},
	{0, 4}, {1, 5}, {2, 6}, {3, 7}
    };
#define BX(i, j) cube->bbox[i][j]
    FLOAT cube_vert[8][3] = {
	{BX(0, 0), BX(1, 0), BX(2, 0)},
	{BX(0, 1), BX(1, 0), BX(2, 0)},
	{BX(0, 1), BX(1, 1), BX(2, 0)},
	{BX(0, 0), BX(1, 1), BX(2, 0)},
	{BX(0, 0), BX(1, 0), BX(2, 1)},
	{BX(0, 1), BX(1, 0), BX(2, 1)},
	{BX(0, 1), BX(1, 1), BX(2, 1)},
	{BX(0, 0), BX(1, 1), BX(2, 1)}
    };
#undef BX

    for (l = 0; l < 12; l++) {
	FLOAT *x0 = cube_vert[cube_edge[l][0]];
	FLOAT *x1 = cube_vert[cube_edge[l][1]];

	FLOAT L0[Dim+1], L1[Dim+1], L[Dim+1];
	const FLOAT *lambda_tmp;

	lambda_tmp = 
	    //geomXYZ2Lambda(tet->J[0], x0[0], x0[1], x0[2]);
	    phgGeomXYZ2Lambda(g_, tet, x0[0], x0[1], x0[2]);
	memcpy(L0, lambda_tmp, (Dim+1) * sizeof(FLOAT));
	lambda_tmp = 
	    //geomXYZ2Lambda(tet->J[0], x1[0], x1[1], x1[2]);
	    phgGeomXYZ2Lambda(g_, tet, x1[0], x1[1], x1[2]);
	memcpy(L1, lambda_tmp, (Dim+1) * sizeof(FLOAT));

	FLOAT t;		/* Line: L0 * t + L1 * (1-t)
				 * intersect (0, x, x, x), (x, 0, x, x), etc */

	for (i = 0; i < NVert; i++) {
	    if (fabs(L1[i] - L0[i]) < 1e-12) {
		/* seg is parallel to face_{i} */
		continue;
	    }

	    t = (0 - L0[i]) / (L1[i] - L0[i]);
	    if (t > 1 + EPS0 || t < -EPS0) {
		/* out of segment */
		continue;
	    }

	    for (k = 0; k < Dim+1; k++)
		L[k] = L0[k] + (L1[k] - L0[k]) * t;
	    
	    if (LambdaInElement(L)) {
		return TRUE;
	    }
	} /* end tet face */
    }	  /* end cube edge */



    return FALSE;
}




/* ------------------------------------------------------ */
/* Geom functions */
static void 
r8vec_cross_3d ( double v1[3], double v2[3] , double v3[3])
{
  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

static double 
r8vec_length ( int dim_num, double x[] )
{
    int i;
    double value;

    value = 0.0;
    for ( i = 0; i < dim_num; i++ )  {
	value = value + pow ( x[i], 2 );
    }
    value = sqrt ( value );

    return value;
}

static double 
r8mat_det_4d ( double a[4*4] )
{
  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;
}

static void 
tetrahedron_insphere_3d(double tetra[3 * 4], double *r, double pc[3])
{
#define DIM_NUM 3

    double b[4 * 4];
    double gamma;
    int i;
    int j;
    double l123;
    double l124;
    double l134;
    double l234;
    double n123[DIM_NUM];
    double n124[DIM_NUM];
    double n134[DIM_NUM];
    double n234[DIM_NUM];
    double v21[DIM_NUM];
    double v31[DIM_NUM];
    double v41[DIM_NUM];
    double v32[DIM_NUM];
    double v42[DIM_NUM];
    double v43[DIM_NUM];

    for (i = 0; i < DIM_NUM; i++) {
	v21[i] = tetra[i + 1 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
    }
    for (i = 0; i < DIM_NUM; i++) {
	v31[i] = tetra[i + 2 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
    }
    for (i = 0; i < DIM_NUM; i++) {
	v41[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 0 * DIM_NUM];
    }
    for (i = 0; i < DIM_NUM; i++) {
	v32[i] = tetra[i + 2 * DIM_NUM] - tetra[i + 1 * DIM_NUM];
    }
    for (i = 0; i < DIM_NUM; i++) {
	v42[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 1 * DIM_NUM];
    }
    for (i = 0; i < DIM_NUM; i++) {
	v43[i] = tetra[i + 3 * DIM_NUM] - tetra[i + 2 * DIM_NUM];
    }

    r8vec_cross_3d(v21, v31, n123);
    r8vec_cross_3d(v41, v21, n124);
    r8vec_cross_3d(v31, v41, n134);
    r8vec_cross_3d(v42, v32, n234);

    l123 = r8vec_length(DIM_NUM, n123);
    l124 = r8vec_length(DIM_NUM, n124);
    l134 = r8vec_length(DIM_NUM, n134);
    l234 = r8vec_length(DIM_NUM, n234);

    for (i = 0; i < DIM_NUM; i++) {
	pc[i] = (l234 * tetra[i + 0 * DIM_NUM]
		 + l134 * tetra[i + 1 * DIM_NUM]
		 + l124 * tetra[i + 2 * DIM_NUM]
		 + l123 * tetra[i + 3 * DIM_NUM])
	    / (l234 + l134 + l124 + l123);
    }

    for (j = 0; j < 4; j++) {
	for (i = 0; i < DIM_NUM; i++) {
	    b[i + j * 4] = tetra[i + j * DIM_NUM];
	}
	b[3 + j * 4] = 1.0;
    }

    gamma = fabs(r8mat_det_4d(b));

    *r = gamma / (l234 + l134 + l124 + l123);

    return;
#undef DIM_NUM
}


/* ------------------------------------------------------ */
/* plot */
BOOLEAN 
tranverse_plot_all CB_ARGS(e) 
{
    int verb = 3;
    DRAW_TET_OF__CUBE(e);    
    return TRUE;
}


static int max_ntet;
static int min_ntet;
static int all_ntet;

BOOLEAN 
tranverse_plot_cube CB_ARGS(e) 
{
    int verb = 3, j;

    vtkSetColor(verb, "white");	
    vtkSetTransparent(verb, .5);	
    vtkDrawBox(verb, e->bbox[0]);

    verb = 3;
    phgInfo(verb, "cube[%5d] [%5d]: ",  
	    _octTreeTraverseIndex, e->ntet);
    for (j = 0; j < e->ntet; j++) {
	phgInfo(verb, "%5d ", e->liste[j]->index);
    }
    phgInfo(verb, "\n");

    /* ntet range */
    if (e->ntet > max_ntet)
	max_ntet = e->ntet;

    if (e->ntet < min_ntet)
	min_ntet = e->ntet;

    all_ntet += e->ntet;

    return TRUE;
}


/* ------------------------------------------------------ */
/* main functions */
_OCT_TREE *
_build_octree(GRID * g)
{
    _OCT_TREE *og;
    FLOAT dx, dy, dz;
    SIMPLEX *tet;
    _CUBE *cube;
    int i, j, k;


    /*************/
    /* Draw mesh */
    /*************/
    if (1) {
	int verb = 0;
	vtkInit(verb, g);
	vtkSetTransparent(verb, .15);
	vtkDrawMesh(verb);
	vtkSetTransparent(verb, 1.);
	vtkSetColor(verb, "blue");
	vtkDrawRemoteFace(verb);
	vtkSetColor(verb, "white");
	vtkSetTransparent(verb, 1.);
	//vtkPause(0);
    }


    g_ = g;

    /* 
     * Step 1.
     *    build all oct roots.
     * */
    og = phgCalloc(1, sizeof(*og));
    og->nx = NX;
    og->ny = NY;
    og->nz = NZ;
    og->nroot = NX * NY * NZ;
    og->roots = phgCalloc(og->nroot, sizeof(*og->roots));

    /* compute local bound box */
    if (g->nvert == 0) {
	for (j = 0; j < Dim; j++) {
	    og->bbox[j][0] = og->bbox[j][1] = 0.;
	}
    } else {
	FLOAT f;

	for (j = 0; j < Dim; j++) {
	    og->bbox[j][0] = og->bbox[j][1] = g->verts[0][j];
	}
	for (i = 1; i < g->nvert; i++) {
	    COORD *c = g->verts + i;
	    for (j = 0; j < Dim; j++) {
		if ((f = (*c)[j]) < og->bbox[j][0]) og->bbox[j][0] = f;
		if (f > og->bbox[j][1]) og->bbox[j][1] = f;
	    }
	}
    }


    phgInfo(0, "Linear box resize\n");
    for (k = 0; k < Dim; k++) {
	og->bbox[k][0] -= 10 * EPS;
	og->bbox[k][1] += 10 * EPS;
    }

    dx = (og->bbox[0][1] - og->bbox[0][0]) / NX;
    dy = (og->bbox[1][1] - og->bbox[1][0]) / NY;
    dz = (og->bbox[2][1] - og->bbox[2][0]) / NZ;
    og->dx = dx;
    og->dy = dy;
    og->dz = dz;


    for (i = 0; i < NX; i++)
	for (j = 0; j < NY; j++)
	    for (k = 0; k < NZ; k++) {
		cube = og->roots + (i * NY + j) * NZ + k;

		cube->bbox[0][0] = og->bbox[0][0] + dx * i;
		cube->bbox[0][1] = og->bbox[0][0] + dx * (i + 1);
		cube->midx[0] = og->bbox[0][0] + dx * (i + .5);

		cube->bbox[1][0] = og->bbox[1][0] + dy * j;
		cube->bbox[1][1] = og->bbox[1][0] + dy * (j + 1);
		cube->midx[1] = og->bbox[1][0] + dy * (j + .5);

		cube->bbox[2][0] = og->bbox[2][0] + dz * k;
		cube->bbox[2][1] = og->bbox[2][0] + dz * (k + 1);
		cube->midx[2] = og->bbox[2][0] + dz * (k + .5);

		{
		    FLOAT dx_ = cube->bbox[0][1] - cube->bbox[0][0];
		    FLOAT dy_ = cube->bbox[1][1] - cube->bbox[1][0];
		    FLOAT dz_ = cube->bbox[2][1] - cube->bbox[2][0];
		    cube->diam = sqrt(dx_ * dx_ + dy_ * dy_ + dz_ * dz_);
		    cube->dh = dx_;
		}

		cube->size = 10;
		cube->liste = phgCalloc(cube->size, sizeof(*cube->liste));
	    }

    FLOAT (*bboxs)[Dim][2];
    bboxs = phgCalloc(g->nprocs, sizeof(*bboxs));
    MPI_Allgather(&og->bbox, sizeof(*bboxs), MPI_BYTE,
		  bboxs, sizeof(*bboxs), MPI_BYTE, g->comm);
    og->bboxs = bboxs;


    /* 
     * Step 2.
     *    allocate tets on the oct tree.
     * */

    /* init element center */
    tet_center = phgCalloc(g->nelem, sizeof(*tet_center));
    tet_inner_diam = phgCalloc(g->nelem, sizeof(*tet_inner_diam));
    ForAllElements(g, tet) {
	FLOAT xmid[Dim] = { 0, 0, 0 };
	for (i = 0; i < NVert; i++) {
	    for (k = 0; k < Dim; k++)
		xmid[k] += .25 * g->verts[tet->verts[i]][k];
	}
	memcpy(tet_center[tet->index], xmid, sizeof(*tet_center));


	FLOAT X[3*4], r, pc[Dim];
	for (i = 0; i < NVert; i++)
	    for (k = 0; k < Dim; k++)
		X[k+Dim*i] = g->verts[tet->verts[i]][k];
	tetrahedron_insphere_3d(X, &r, pc);
	phgInfo(3, "inner r: %e\n", r);
	tet_inner_diam[tet->index] = r;
    }




    phgInfo(0, "Elem in Root cubes\n");
    ForAllElements(g, tet) {
	int verb = 3;		/* tet intersect cube */

	for (i = 0; i < og->nroot; i++) {
	    cube = &og->roots[i];

	    if (test_intersection(tet, cube)) {

		vtkTmpActorsBegin(verb);
		vtkSetColor(verb, "white");
		vtkSetTransparent(verb, .25);
		vtkDrawElement(verb, tet);
		vtkSetColor(verb, "green");
		vtkSetTransparent(verb, .25);
		vtkDrawBox(verb, cube->bbox[0]);
		vtkPause(verb);

		vtkTmpActorsClear(verb);

		/* found */
		phgInfo(verb, "Found %d %d.\n", tet->index,
			cube - og->roots);
		if (cube->ntet >= cube->size) {
		    cube->liste =
			phgRealloc_(cube->liste,
				    2 * cube->size * sizeof(SIMPLEX *),
				    cube->size * sizeof(SIMPLEX *));
		    cube->size *= 2;
		}

		cube->liste[cube->ntet++] = tet;
	    }			/* end of intersect */
	}			/* end of oct tree roots */
    }				/* end of tets */


    /* Check tets in cube  */
    for (i = 0; i < og->nroot; i++) {
	int verb = 3;
	cube = &og->roots[i];

	vtkTmpActorsBegin(verb);
	vtkSetColor(verb, "red");
	vtkSetTransparent(verb, .75);
	vtkDrawBox(verb, cube->bbox[0]);

	phgInfo(verb, "cube_%5d [%5d]: ", i, cube->ntet);
	for (j = 0; j < cube->ntet; j++) {
	    phgInfo(verb, "%5d ", cube->liste[j]->index);
	    vtkSetColor(verb, "white");
	    vtkSetTransparent(verb, .25);
	    vtkDrawElement(verb, cube->liste[j]);
	}
	phgInfo(verb, "\n");

	vtkPause(verb);
	vtkTmpActorsClear(verb);
    }


    /* cube list */
    ncube_list = 0;
    ncube_size = 8;
    cube_list = phgCalloc(ncube_size, sizeof(*cube_list));
    
    /* tranvers insert */
    for (i = 0; i < og->nroot; i++) {
	_octTreeTraverseDepth = 0;
	if (!traverse_insert(og->roots + i, NULL))
	    break;
    }



    /*
     * 
     * tranvers plot
     *
     * */
    if (0) {
	/* Plot all tets on cube*/
	traverse_all_flag = FALSE;
	_octTreeTraverseElements(og, tranverse_plot_all);
    }
    
    if (1) {
	/* All leaf cubes */
	int verb = 3;
	vtkTmpActorsBegin(verb);

	max_ntet = -100000;
	min_ntet = 100000;
	all_ntet = 0;

	traverse_all_flag = FALSE;
	_octTreeTraverseElements(og, 
				tranverse_plot_cube);
	og->ncube = _octTreeTraverseIndex;
	//vtkPause(verb);

	vtkTmpActorsClear(verb);

	phgInfo(0, "ntet per cube: [%d %d], ave %f\n", 
		min_ntet, max_ntet, 
		all_ntet / (FLOAT) og->ncube);
    }

    og->tet_center = tet_center;
    og->tet_inner_diam = tet_inner_diam;
    tet_center = NULL;
    tet_inner_diam = NULL;
    g_ = NULL;

    og->cube_list = cube_list;
    og->ncube_list = ncube_list;
    cube_list = NULL;
    ncube_list = 0;
    ncube_size = 0;

    
    og->g = g;

    return og;
}



BOOLEAN
_locate_point(_OCT_TREE *og, FLOAT *pts, SIMPLEX **e_ptr, FLOAT *lambda)
/* Locate points using oct tree.
 * If found, return true and the element and lambda.
 * */
{
    int i, j, k, l, iroot, ichild;
    FLOAT dx = og->dx, dy = og->dy, dz = og->dz;
    static BOOLEAN visit_tet[1000];

    //assert(g_ == og->g);	/* search grid one by one. */
    assert(g_ == NULL);
    g_ = og->g;

    tet_center = og->tet_center;
    tet_inner_diam = og->tet_inner_diam;

    for (k = 0; k < Dim; k++)
	if (pts[k] < og->bbox[k][0] - EPS 
	    || pts[k] > og->bbox[k][1] + EPS) {
	    phgInfo(3, "Out of box\n");

	    if (g_->nprocs == 1) {
		int verb = 0;
		vtkSetColor(verb, "red");
		vtkSetTransparent(verb, 1.);
		vtkDrawPoint(verb, pts, "  feet");
		vtkPause(verb);

		phgInfo(0, "Out of box, no found\n");
		show_V_(0, pts, Dim, "point", "%30.15f");
		phgError(-1, "Not found???\n");
	    }
	    tet_center = NULL;
	    tet_inner_diam = NULL;
	    g_ = NULL;
	    return FALSE;
	}

    i = (int) ((pts[0] - og->bbox[0][0]) / dx);
    j = (int) ((pts[1] - og->bbox[1][0]) / dy);
    k = (int) ((pts[2] - og->bbox[2][0]) / dz);
    if (i == NX) i = NX-1;
    if (j == NY) j = NY-1;
    if (k == NZ) k = NZ-1;

    iroot = (i * NY + j) * NZ + k;

    /* search oct tree */
    _CUBE *cube = og->roots + iroot;
    _octTreeTraverseIndex = 0;
    while (TRUE) {
	int verb = 3;

	vtkSetColor(verb, "white");
	vtkSetTransparent(verb, .75);
	vtkDrawBox(verb, cube->bbox[0]);

	if (cube->children[0] != NULL) {

	    /* go deeper,  locate child */
	    i = 0, j = 0, k = 0;
	    if (pts[0] > cube->midx[0]) i = 1;
	    if (pts[1] > cube->midx[1]) j = 1;
	    if (pts[2] > cube->midx[2]) k = 1;

	    ichild = 4 * i + 2 * j + k;
	    cube = cube->children[ichild];
	    _octTreeTraverseIndex++;

	} else {
	    /* leaf cube */

#if 0
	    /* ----------------------
	     *
	     *
	     *   Search tets neighs 
	     *
	     *
	     * ---------------------- */
	    static BOOLEAN tested[1000];
	    int to_test = -1;	/* tet to test */

	    bzero(tested, sizeof(tested));
	    while (TRUE) {
		if (to_test < 0) {
		    for (j = 0; j < cube->ntet; j++) {
			if (tested[j] == FALSE) {/* only test untouched */
			    to_test = j;
			    break;
			}
		    }
		    if (j == cube->ntet) /* all element tested */
			break;
		}
		
		assert(tested[to_test] == FALSE);
		SIMPLEX *tet = cube->liste[to_test];
		const FLOAT *lambda_tmp;

		tested[to_test] = TRUE;
		lambda_tmp = phgGeomXYZ2Lambda(g_, tet, 
					       pts[0], pts[1], pts[2]
					       );

		phgInfo(3, "cube tet[%2d]\n", to_test);

		if (LambdaInElement(lambda_tmp)) { /* lie in test replaced in geom. */
		    /* found */
		    *e_ptr = tet;
		    memcpy(lambda, lambda_tmp, 
			   (Dim+1) * sizeof(FLOAT));
		    tet_center = NULL;
		    tet_inner_diam = NULL;
		    g_ = NULL;
		    return TRUE;
		}
		else {
		    /* not found, goto neigh */
		    FLOAT minL = 2.;
		    int   imin;
		    for (k = 0; k < Dim+1; k++) {
			if (lambda_tmp[k] < minL) {
			    minL = lambda_tmp[k];
			    imin = k;
			}
		    } /* end find min */

		    /* goto next neigh of face[imin] */
		    if (HasLocalNeighbour(tet, imin)) {
			SIMPLEX *tet_ngh = 
			    (SIMPLEX *)tet->neighbours[imin];

			for (l = 0; l < cube->ntet; l++) {
			    if (tet_ngh->index == cube->liste[l]->index
				&& tested[l] == FALSE) {
				to_test = l;  /* neigh in list and not tested*/
				break;
			    }
			}
			if (l == cube->ntet) 
			    to_test = -1; /* neigh not in list */
		    } else {
			/* no local neigh, keep on search in list */
			to_test = -1;
		    }

		} /* end not found */


	    } /* end main loop */
#else
	    /* ----------------------------
	     *
	     *
	     *   Search tets one by one
	     *
	     *
	     * ---------------------------- */
	    verb = 3;
	    phgInfo(verb, "----------\n");
	    phgInfo(verb, "%30.15e %30.15e %30.15e, ntet: %d\n", 
		    pts[0], pts[1], pts[2], cube->ntet); 

	    BOOLEAN watch = FALSE;
#  if 0
	    if (0) {
		/* Debug:
		 *   Watch particular point */
		FLOAT tmp[3] = {7.7903907113750800175679556841545903012047e-07,
				-4.5733511807839383189744353330752346664667e-01,
				1.0000102611938554986181770800612866878510e+00};
		for (k = 0; k < Dim; k++) 
		    tmp[k] -= pts[k];
		if (sqrt(INNER_PRODUCT(tmp, tmp) < 1e-12))
		    watch = TRUE;
		if (watch) 
		    phgInfo(0, "watching...\n");
	    }
#  endif


	    for (j = 0; j < cube->ntet; j++) {
		SIMPLEX *tet = cube->liste[j];
		const FLOAT *lambda_tmp;
		
		if (watch) 
		    phgInfo(0, "search elem: %d\n", tet->index);

		phgInfo(verb, "\nsearch %d\n", j);
		lambda_tmp = phgGeomXYZ2Lambda(g_, tet, 
					       pts[0], pts[1], pts[2]
					       );

		if (watch) {
		    phgInfo(0, "   lambda: %30.15f %30.15f %30.15f %30.15f\n", 
			    lambda_tmp[0], lambda_tmp[1], lambda_tmp[2], lambda_tmp[3]);
		    phgInfo(0, "   lie in test: %d\n", LambdaInElement(lambda_tmp));
		}

		phgInfo(verb, "  cube tet[%2d]\n", j);
		show_V_(verb, lambda_tmp, Dim+1, "lambda", "%f");

		if (LambdaInElement(lambda_tmp)) { 
		    /* found */
		    *e_ptr = tet;
		    memcpy(lambda, lambda_tmp, 
			   (Dim+1) * sizeof(FLOAT));
		    tet_center = NULL;
		    tet_inner_diam = NULL;
		    g_ = NULL;
		    return TRUE;
		}
	    } /* end main loop */
	    phgInfo(verb, "  Not found in cube\n");
#endif


	    /*
	     *
	     *
	     * Locate pts for one proc.
	     * Must found
	     *
	     *
	     * */
	    if (g_->nprocs == 1) {
		int verb = 0;
		vtkSetColor(verb, "red");
		vtkSetTransparent(verb, 1.);
		vtkDrawPoint(verb, pts, "  feet");
		vtkPause(verb);

		show_V_(0, pts, Dim, "not found", "%30.15e");
		

		SIMPLEX *tet;
		ForAllElements(g_, tet) {
		    const FLOAT *lambda_tmp;
		    lambda_tmp = 
			phgGeomXYZ2Lambda(g_, tet, 
					  pts[0], pts[1], pts[2] 
					  );

		    if (LambdaInElement(lambda_tmp)) {
			vtkSetColor(verb, "yellow");
			vtkSetTransparent(verb, 1.);
			vtkDrawElement(verb, tet);
			
			/* check */
			printf("inter test %d (%d)\n", 
			       test_intersection(tet, cube), 
			       tet->index);

			show_V_(0, lambda_tmp, Dim+1, "crv lambda", "%f");
			break;
		    }
		}
		vtkSetColor(verb, "red");
		vtkSetTransparent(verb, .15);
		vtkDrawBox(verb, cube->bbox[0]);	
		vtkPause(verb);

		
		DRAW_TETS_OF__CUBE_(cube, 1);
		/* error */
		show_V_(0, pts, Dim, "point", "%30.15f");
		phgError(-1, "Not found???\n");
	    }

	    /* not found */
	    *e_ptr = NULL;
	    tet_center = NULL;
	    tet_inner_diam = NULL;
	    g_ = NULL;
	    return FALSE;
	} 

    }

    tet_center = NULL;
    tet_inner_diam = NULL;
    g_ = NULL;
    return FALSE;
}



static BOOLEAN
free_cube(_CUBE *e)
{
    phgFree(e->liste);
    
    return TRUE;
}

void 
_destroy_octree(_OCT_TREE **og_ptr)
{
    _OCT_TREE *og = *og_ptr;
    int i;

    //assert(g_ == og->g);
    assert(g_ == NULL);

    phgFree(tet_center);
    tet_center = NULL;

    phgFree(tet_inner_diam);
    tet_inner_diam = NULL;

    /* Recursively free cubes */
    _octTreeTraverseAllElements(og, free_cube);
    for (i = 0; i < og->ncube_list; i++)
	phgFree(og->cube_list[i]);
    phgFree(og->cube_list);

    phgFree(og->bboxs);
    phgFree(og->roots);

    phgFree(og);
    *og_ptr = NULL;
}




