/*
 * Method of characteristic.
 *
 * */

#ifndef OCT_SEARCH_H


#include "phg.h"
//#include "vtk-draw.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>


/***************/
/* Parameters  */
/***************/

#define EPS_LIE_IN   2.5E-12
#define EPS          1e-13
#define EPS0         1e-14		/* normalized */

#define NX 3			/* root cubes in x-dir */
#define NY 3
#define NZ 3
#define MAX_NELEM  8		/* Cube do NOT need to divide if it has only 8 tets */
#define MAX_NELEM0 60		/* Cube do NOT need to divide if it has so many tets,
				 *   Not used*/
#define MAX_NLEVEL 9


/*******************/
/* Debug utilities */
/*******************/

#define DEBUG_SEARCH 1		/* print control */


/**********************************/
/* Structure & function prototype */
/**********************************/


typedef FLOAT BBOX[Dim][2];

typedef struct _CUBE_ {
    struct _CUBE_ *children[2 * 2 * 2];			/**< Pointers to children */
    struct _CUBE_ *father;	/* father */
    int level;			/* oct tree level */
    FLOAT dh;			/* element size */
    FLOAT bbox[Dim][2];		/* bounding box */
    FLOAT midx[Dim];		/* mid point */
    FLOAT diam;

    int ntet, size;
    SIMPLEX **liste;

} _CUBE;

typedef struct _OCT_TREE_ {
    int nx, ny, nz;
    int ncube;
    int nroot;
    FLOAT dx, dy, dz;
    _CUBE *roots;
    FLOAT bbox[Dim][2];		/* bounding box */
    FLOAT (*bboxs)[Dim][2];

    GRID *g;			/* working tree */
    COORD *tet_center;		/* extra tets info */
    FLOAT *tet_inner_diam;	/* extra tets info */
    _CUBE **cube_list;		/* record new cubes alloc */
    INT ncube_list;

} _OCT_TREE;

/* static COORD *tet_center; */
/* static FLOAT *tet_inner_diam; */
/* static GRID *g_;		/\* current grid *\/ */
/* static _CUBE **cube_list; */
/* static int ncube_list;  	/\* # of cube in lists *\/ */
/* static int ncube_size; */


typedef void (*FUNC_BDRY)(FLOAT x, FLOAT y, FLOAT z, FLOAT t, FLOAT *value);

/* oct tree search */
_OCT_TREE *_build_octree(GRID * g);
BOOLEAN _locate_point(_OCT_TREE *og, FLOAT *pts, SIMPLEX **e_ptr, FLOAT *lambda);
void _destroy_octree(_OCT_TREE **og_ptr);




/**************/
/* Misc utils */
/**************/
#define COORD_SIZE (Dim*sizeof(FLOAT))
#define LAMBDA_SIZE ((Dim+1)*sizeof(FLOAT))





#define LIE_IN(x) ((x) >= (- EPS_LIE_IN) && (x) <= (1. + EPS_LIE_IN)) 
#define LambdaInElement(lambda) (LIE_IN(lambda[0]) && LIE_IN(lambda[1]) && \
				 LIE_IN(lambda[2]) && LIE_IN(lambda[3]))		

#define SQUARE(x) ((x)*(x))
#define INNER_PRODUCT(p, q)			\
    (*(p    ) * *(q    ) +			\
     *(p + 1) * *(q + 1) +			\
     *(p + 2) * *(q + 2))
#define VEC_LENGTH(v) sqrt(INNER_PRODUCT(v, v))

#define BUF_ADD_NUM 100
#define INITIAL_BUF(buf, ndata, pos, buf_size) {			\
	pos = 0;							\
	buf = phgAlloc((buf_size = BUF_ADD_NUM) * ndata * sizeof(*buf)); \
    }

#define REALLOC_BUF(buf, ndata, pos, buf_size, add_size)		\
    if ((pos + add_size - 1)>= buf_size) {				\
	size_t buf_size0_ = buf_size;					\
	buf_size *= 2;							\
	DebugPrint(verb+1, "Resize " #buf " list: %d --> %d\n",		\
		   buf_size0_ * ndata, buf_size * ndata);		\
	buf = phgRealloc_(buf, buf_size * ndata * sizeof(*buf),		\
			  buf_size0_ * ndata * sizeof(*buf));		\
    }
/* Utils */
#define DUMP_BOX(bbox)					\
    phgInfo(0, "bbox: [%f %f] x [%f %f] x [%f %f]\n",	\
	    bbox[0][0], bbox[0][1],		\
	    bbox[1][0], bbox[1][1],		\
	    bbox[2][0], bbox[2][1]		\
	    );



/*******************/
/* Debug utilities */
/*******************/

extern BOOLEAN debug_moc_actived;
#if DEBUG_SEARCH
# define DebugPrint(level, fmt, ...)					\
    if (debug_moc_actived)						\
	phgInfo((level), "%s" fmt, LEVEL[(level+1)], ##__VA_ARGS__);

# define show_fV(level, vec, vec_n, desp)		\
    if (debug_moc_actived) {				\
	show_V_(level, vec, vec_n, desp, "%24.12E");	\
    } 
# define show_iV(level, vec, vec_n, desp)		\
    if (debug_moc_actived)				\
	show_V_(level, vec, vec_n, desp, "%10d")

#if 0
# define show_V_(level, vec, vec_n, desp, fmt) { int i_;		\
	FILE *fp = phgGetInfoFile(level);				\
	if (debug_moc_actived && fp != NULL) {				\
	    phgInfo(level, "%s%-10s:(%3d) L%5d: [",			\
		    LEVEL[level+1], desp, (vec_n), __LINE__);		\
	    for (i_ = 0; i_ < (vec_n); i_++) {				\
		fprintf(fp, fmt", ", *(vec + i_));			\
		if (vec_n > 6 && (i_+1) % 6 == 0)			\
		    fprintf(fp, "\n%s%s", LEVEL[5], LEVEL[level]);	\
	    }								\
	    fprintf(fp, "]\n");						\
	}								\
    }
#else 
# define show_V_(level, vec, vec_n, desp, fmt) { int i_;		\
	phgInfo(level, "%s%-10s:(%3d) L%5d: [",				\
		LEVEL[level+1], desp, (vec_n), __LINE__);		\
	for (i_ = 0; i_ < (vec_n); i_++) {				\
	    phgInfo(level, fmt", ", *(vec + i_));			\
	    if (vec_n > 6 && (i_+1) % 6 == 0)				\
		phgInfo(level, "\n%s%s", LEVEL[5], LEVEL[level]);	\
	}								\
	phgInfo(level, "]\n");						\
    }								
#endif			

# define show_eV_(level, vec, vec_n, desp) { int i_;		\
	FILE *fp = phgGetInfoFile(level);			\
	if (debug_moc_actived && fp != NULL) {			\
	    phgInfo(level, "%s%-15s:(%3d) L%5d: [",		\
		    LEVEL[level], desp, (vec_n), __LINE__);	\
	    for (i_ = 0; i_ < (vec_n); i_++) {			\
		fprintf(fp, "%4d, ", vec[i_]->index);		\
		if (vec_n > 6 && (i_+1) % 6 == 0)		\
		    fprintf(fp, "\n%s", LEVEL[level]);		\
	    }							\
	    fprintf(fp, "]\n");					\
	}							\
    }			

# define CheckLambdaN(lambda, len) {					\
	FLOAT sum = 0.;							\
	int k_;								\
	for (k_ = 0; k_ < len; k_++)					\
	    sum += lambda[k_];						\
	if(fabs(sum-1.) > EPS_LAM_SUM) {				\
	    show_fV(1, lambda, Dim+1, "!!! Err in lambda");		\
	    phgInfo(0, "sum lambda(%d) err = %E !!!", len, 1.-sum);	\
	}								\
    }

# define CheckLambda(lambda) CheckLambdaN(lambda, 4)
//# define CheckLambda(lambda) 
#else 
# define DebugPrint(level, fmt, ...) ;
# define show_fV(level, vec, vec_n, desp) ;
# define show_iV(level, vec, vec_n, desp) ;
# define show_V_(level, vec, vec_n, desp, fmt) ;
# define show_eV_(level, vec, vec_n, desp) ;
# define CheckLambdaN(lambda, len) ;
# define CheckLambda(lambda) ;
#endif /* DEBUG_SEARCH */

#ifndef DOF_SCALE
# define DOF_SCALE(u, description) {					\
	char trimed[100];						\
	strncpy(trimed, __FUNCTION__, 8);				\
	trimed[8]='\0';							\
	phgPrintf("   ------------------------------\n"			\
		  "   %-10s  /* %s */\n"				\
		  "     func: %-10s, line: %03d\n"			\
		  ,#u":", description,					\
		  trimed, __LINE__);					\
	phgPrintf("     [ %16.8e, %16.8e, %16.8e ] ( L2, L1, Linf); \n",	\
		  phgDofNormL2(u), phgDofNormL1(u), 			\
		  phgDofNormInftyVec(u)					\
		  );							\
    }
#endif



/*************/
/* Vtk Debug */
/*************/
/* Get this using script
 * 1. define
 * 2. void
 * */
#if 0 				/* Vtk Debug */
#  warning --------- Vtk debug enabled -----------

# if 0
#   define VTK_ACT(func) fprintf(stderr, "   %s works!\n", func);
# else
#   define VTK_ACT(func)
# endif

#  define vtkDrawEdge(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawEdge(__VA_ARGS__);}}
#  define vtkDrawElement(verb, ...)     {if ((verb) <= phgVerbosity) { vtk.DrawElement(__VA_ARGS__);}}
#  define vtkDrawFace(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawFace(__VA_ARGS__);}}
#  define vtkDrawLine(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawLine(__VA_ARGS__);}}
#  define vtkDrawRemoteFace(verb, ...)  {if ((verb) <= phgVerbosity) { vtk.DrawRemoteFace(__VA_ARGS__);}}
#  define vtkDrawMesh(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawMesh(__VA_ARGS__);}}
#  define vtkDrawPoint(verb, ...)       {if ((verb) <= phgVerbosity) { vtk.DrawPoint(__VA_ARGS__);}}
#  define vtkhi(verb, ...)              {if ((verb) <= phgVerbosity) { vtk.hi(__VA_ARGS__);}}
#  define vtkInit(verb, ...)            {if ((verb) <= phgVerbosity) { vtk.Init(__VA_ARGS__);}}
#  define vtkPause(verb, ...)           {					\
	if ((verb) <= phgVerbosity) {					\
	    fprintf(stderr, "*%2d* vtk pause \n\tfile:%-25s line:%5d\n\tfunc:%-25s\n\n", \
		    g->rank, __FILE__, __LINE__, __FUNCTION__);		\
	    vtk.Pause(__VA_ARGS__);					\
	}								\
    }
#  define vtkSetColor(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.SetColor(__VA_ARGS__);}}
#  define vtkSetTransparent(verb, ...)  {if ((verb) <= phgVerbosity) { vtk.SetTransparent( __VA_ARGS__ );}}
#  define vtkTmpActorsBegin(verb, ...)  {if ((verb) <= phgVerbosity) { vtk.TmpActorsBegin( __VA_ARGS__ );}}
#  define vtkTmpActorsClear(verb, ...)  {if ((verb) <= phgVerbosity) { vtk.TmpActorsClear( __VA_ARGS__ );}}

#elif 1
#  warning --------- Vtk debug disabled -----------
#  define vtkDrawEdge(verb, ...)
#  define vtkDrawElement(verb, ...)
#  define vtkDrawFace(verb, ...) 
#  define vtkDrawLine(verb, ...) 
#  define vtkDrawRemoteFace(verb, ...)
#  define vtkDrawMesh(verb, ...) 
#  define vtkDrawPoint(verb, ...) 
#  define vtkDrawBox(verb, ...) 
#  define vtkhi(verb, ...) 
#  define vtkInit(verb, ...) 
#  define vtkPause(verb, ...) 
#  define vtkSetColor(verb, ...) 
#  define vtkSetTransparent(verb, ...) 
#  define vtkTmpActorsBegin(verb, ...) 
#  define vtkTmpActorsClear(verb, ...) 
#endif /* Vtk Debug */

#define DUMP_ELEMENT(e, end) {				\
	phgInfo(0, "Element search failed: \n");	\
	phgInfo(0, "Crd=[\n");				\
	phgInfo(0, "%24.12E, %24.12E, %24.12E;\n",	\
		g->verts[e->verts[0]][0],		\
		g->verts[e->verts[0]][1],		\
		g->verts[e->verts[0]][2]);		\
	phgInfo(0, "%24.12E, %24.12E, %24.12E;\n",	\
		g->verts[e->verts[1]][0],		\
		g->verts[e->verts[1]][1],		\
		g->verts[e->verts[1]][2]);		\
	phgInfo(0, "%24.12E, %24.12E, %24.12E;\n",	\
		g->verts[e->verts[2]][0],		\
		g->verts[e->verts[2]][1],		\
		g->verts[e->verts[2]][2]);		\
	phgInfo(0, "%24.12E, %24.12E, %24.12E;\n",	\
		g->verts[e->verts[3]][0],		\
		g->verts[e->verts[3]][1],		\
		g->verts[e->verts[3]][2]);		\
	phgInfo(0, "%24.12E, %24.12E, %24.12E;\n",	\
		end[0],					\
		end[1],					\
		end[2]);				\
	phgInfo(0, "];\n");				\
    }



#define OCT_SEARCH_H
#endif
