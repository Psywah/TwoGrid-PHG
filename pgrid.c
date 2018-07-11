/* Parallel grid import
 *
 * 1. Input:
 *    num of verts
 *    num of elems
 *    elems->verts[4] (local)
 *    elem local index
 *    vert map of local to global 
 *    elem map of local to global
 *    elem face mark:
 *      [interior|boundary|remote]
 *    
 *
 *
 * 2. Output:
 *
 *
 * */

#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "phg.h"
#include "pgrid.h"
//#include "ins.h"


#define SHOW_VEC(vec, n) {				\
	phgInfo(0, #vec "   %d\n", __LINE__);		\
	INT i_;						\
	for (i_ = 0; i_ < n; i_ ++) {			\
	    phgInfo(0, "   [%3d]: %5d\n", i_, vec[i_]);	\
	}						\
    }


#define MEM_DUP(dest, src, count)				\
    phgInfo(2, "# mem copy: (%s) <=== (%s)\n", #src, #dest);	\
    if(src == NULL) {						\
	dest = NULL;						\
	phgInfo(2, "#   NULL copyed\n");			\
    } else {							\
	dest = phgCalloc((count), sizeof(src[0]));		\
	memcpy(dest, src, (count)*sizeof(src[0]));		\
	phgInfo(2, "#   %d copyed\n", (count));			\
    }





/* --------------------------------------------------------------------------------
 *
 *
 *
 *    Set edge index
 *    
 *
 *
 * --------------------------------------------------------------------------------
 */
typedef INT EdgeVerts[2];
typedef struct EdgeVertsSort_ {
    INT v0;			/* vert */
    INT v1;			/* vert */
    INT ridx;			/* recv index */
    INT sidx;			/* uniq index */
} EdgeVertsSort;


static int
phgCompEdge(const void *p0, const void *p1)
/* compares two 'INT's (used in qsort and bsearch) */
{
    int i;
    EdgeVerts *vts0 = (EdgeVerts *) p0;
    EdgeVerts *vts1 = (EdgeVerts *) p1;

    i = (*vts0)[0] - (*vts1)[0];
    
    if (i < 0) {
	return -1;
    } else if (i > 0) {
	return 1;
    } else {
	i = (*vts0)[1] - (*vts1)[1];
	return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
    }

}

static int
phgCompEdge2(const void *p0, const void *p1)
/* compares two 'INT's (used in qsort and bsearch) */
{
    int i;
    EdgeVertsSort *vts0 = (EdgeVertsSort *) p0;
    EdgeVertsSort *vts1 = (EdgeVertsSort *) p1;

    i = (*vts0).v0 - (*vts1).v0;
    
    if (i < 0) {
	return -1;
    } else if (i > 0) {
	return 1;
    } else {
	i = (*vts0).v1 - (*vts1).v1;
	return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
    }

}

#define SORT_VERT(v0, v1) {			\
	int vv;					\
	if (v0 > v1) {				\
	    vv = v0; v0 = v1; v1 = vv;		\
	}					\
    }



/*
 * Set edge index
 * 1. set local edge index
 * 2. mark remote vert then edge
 * 3. alltoall determine global index remote edge
 *
 * */
static void
set_edge_index(GRID *g) 
{
    INT i, j, k;
    SIMPLEX *e;

    /* ------------------------------------------------------------
     * 
     * Step 1. set local edge index
     *
     * ------------------------------------------------------------ */
    /* Note:
     *  phgUpdateEdges(g); is not safe.
     *  Since the local grid is not big, just use slow O(N^2) algorithm
     * */
    EdgeVerts *edge2vert;
    INT count, size;

    count = 0;
    edge2vert = phgCalloc((size = 100), sizeof(EdgeVerts));
    ForAllElements(g, e) {
	for (i = 0; i < NEdge; i++) {
	    INT v0 = e->verts[GetEdgeVertex(i, 0)];
	    INT v1 = e->verts[GetEdgeVertex(i, 1)];
	    INT vv;

	    SORT_VERT(v0, v1);

	    if (count >= size) {
		edge2vert = phgRealloc_(edge2vert,
					2*size * sizeof(EdgeVerts),
					size * sizeof(EdgeVerts));
		size *= 2;
	    }
	    edge2vert[count][0] = v0;
	    edge2vert[count][1] = v1;
	    count++;
	}
    }
    
    /* remove duplicate */
    qsort(edge2vert, count, sizeof(EdgeVerts), phgCompEdge);

    {
	INT i0 = 0;
	for (i = i0+1; i < count; i++) {
	    int cmp = phgCompEdge(edge2vert + i0,
				    edge2vert + i);
	    if (cmp < 0) {
		i0++;
		memcpy(edge2vert + i0,
		       edge2vert + i,
		       sizeof(EdgeVerts));
	    }
	}
	count = i0 + 1;
    }
    g->nedge = count;
    phgInfo(0, "nedge: %d\n", g->nedge);

    ForAllElements(g, e) {
	for (i = 0; i < NEdge; i++) {
	    INT v0 = e->verts[GetEdgeVertex(i, 0)];
	    INT v1 = e->verts[GetEdgeVertex(i, 1)];

	    SORT_VERT(v0, v1);
	    
	    EdgeVerts vts, *p;
	    vts[0] = v0;
	    vts[1] = v1;
	    p = bsearch(&vts, edge2vert,
			count, sizeof(EdgeVerts),
			phgCompEdge);
	    assert(p != NULL);

	    e->edges[i] = p - edge2vert;
	}
    }
    


    /* ------------------------------------------------------------
     * 
     * Step 2. mark remote edge
     *
     * ------------------------------------------------------------ */
    BYTE *vert_mark, *edge_mark;
    INT *edge_local, *edge_remote;
    INT n_edge_local = 0, n_edge_remote = 0;

    vert_mark = phgCalloc(g->nvert, sizeof(*vert_mark));
    edge_mark   = phgCalloc(count, sizeof(*edge_mark));
    edge_local  = phgCalloc(count, sizeof(*edge_local));
    edge_remote = phgCalloc(count, sizeof(*edge_remote));
    g->L2Gmap_edge = phgCalloc(count, sizeof(*g->L2Gmap_edge));
    
    for (i = 0; i < count; i++)
	g->L2Gmap_edge[i] = -1;

    ForAllElements(g, e) {
	for (i = 0; i < NFace; i++) {
	    if (e->bound_type[i] & REMOTE) {
		int v0 = e->verts[GetFaceVertex(i, 0)];
		int v1 = e->verts[GetFaceVertex(i, 1)];
		int v2 = e->verts[GetFaceVertex(i, 2)];
		vert_mark[v0] = 1;
		vert_mark[v1] = 1;
		vert_mark[v2] = 1;
	    }
	}
    }
    
    for (i = 0; i < count; i++) {
	INT v0 = edge2vert[i][0];
	INT v1 = edge2vert[i][1];
	if (vert_mark[v0] == 1
	    && vert_mark[v1] == 1) { /* both vert remote */
	    edge_remote[n_edge_remote++] = i;
	    edge_mark[i] = 1;	     /* edge remote */
	} else {
	    edge_local[n_edge_local++] = i;
	    edge_mark[i] = 0;	     /* edge local */
	}
    }

    INT rank, nsend, nrecv, 
	*scnts, *sdsps, *rcnts, *rdsps; 
    EdgeVerts *bsend, *brecv;
    scnts = phgCalloc(4 * g->nprocs, sizeof(*scnts));
    sdsps = scnts + g->nprocs;
    rcnts = sdsps + g->nprocs;
    rdsps = rcnts + g->nprocs;

    INT nd, nm, nn;
    for (i = 0; i < count; i++) {
	if (edge_mark[i] == 1) {
	    nn = g->L2Gmap_vert[edge2vert[i][0]]
		+ g->L2Gmap_vert[edge2vert[i][1]];
	    nd = nn / g->nprocs;
	    nm = nn % g->nprocs;

	    rank = nm;
	    scnts[rank]++;
	}
    }

    MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, g->comm);

    nsend = nrecv = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
	sdsps[rank] = nsend;
	rdsps[rank] = nrecv;
	nsend += scnts[rank];
	nrecv += rcnts[rank];
    }
    
    bsend = phgAlloc((nsend + nrecv) * sizeof(*bsend));
    brecv = bsend + nsend;

    for (i = 0; i < count; i++) {
	if (edge_mark[i] != 0) {
	    nn = g->L2Gmap_vert[edge2vert[i][0]]
		+ g->L2Gmap_vert[edge2vert[i][1]];
	    nd = nn / g->nprocs;
	    nm = nn % g->nprocs;

	    rank = nm;
	    INT V0 = g->L2Gmap_vert[edge2vert[i][0]];
	    INT V1 = g->L2Gmap_vert[edge2vert[i][1]];

	    SORT_VERT(V0, V1);

	    bsend[sdsps[rank]][0] = V0;
	    bsend[sdsps[rank]][1] = V1;
	    sdsps[rank]++;
	}
    }

    /* restore sdsps[] */
    for (rank = 0; rank < g->nprocs; rank++)
	sdsps[rank] -= scnts[rank];
    
    /* send lists */
    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(*bsend), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(bsend, scnts, sdsps, type,
		 brecv, rcnts, rdsps, type, g->comm);

    /* process received lists */
    EdgeVertsSort *recv_sort;
    recv_sort = phgCalloc(nrecv, sizeof(*recv_sort));
    for (i = 0; i < nrecv; i++) {
	recv_sort[i].v0 = brecv[i][0];
	recv_sort[i].v1 = brecv[i][1];
	recv_sort[i].ridx = i;
    }

    qsort(recv_sort, nrecv, sizeof(*recv_sort), phgCompEdge2);


    /* count edge number */
    INT remote_edge_count = 0;
    if (nrecv > 0) {
	INT i0 = 0;
	recv_sort[0].sidx = 0;
	for (i = 1; i < nrecv; i++) {
	    int cmp = phgCompEdge2(recv_sort + i - 1,
				       recv_sort + i);
	    if (cmp < 0) {
		i0++;
	    }
	    recv_sort[i].sidx = i0;
	}
	remote_edge_count = i0 + 1;
    }


    /* Reassign local edge number */
    INT remote_start = 0, n_remote_all = 0;
    INT local_start = 0, n_local_all = 0;
    {
	INT n_edge_locals[g->nprocs];
	MPI_Allgather(&n_edge_local, 1, PHG_MPI_INT, 
		      &n_edge_locals, 1, PHG_MPI_INT, g->comm);
	
	for (i = 0; i < g->rank; i++)
	    local_start += n_edge_locals[i];
	for (i = 0; i < g->nprocs; i++)
	    n_local_all += n_edge_locals[i];
    }

    /* Reassign remote edge number */
    {
	INT remote_edge_counts[g->nprocs];
	MPI_Allgather(&remote_edge_count, 1, PHG_MPI_INT, 
		      &remote_edge_counts, 1, PHG_MPI_INT, g->comm);
	
	for (i = 0; i < g->rank; i++)
	    remote_start += remote_edge_counts[i];
	remote_start += n_local_all;
	for (i = 0; i < g->nprocs; i++)
	    n_remote_all += remote_edge_counts[i];
    }

    for (i = 0; i < nrecv; i++) {
	INT i0 = recv_sort[i].ridx;
	brecv[i0][0] = brecv[i0][0] + brecv[i0][1]; /* sum check */
	brecv[i0][1] = recv_sort[i].sidx + remote_start;
    }

    /* send back */
    MPI_Alltoallv(brecv, rcnts, rdsps, type,
    		 bsend, scnts, sdsps, type, g->comm);

    
    for (i = 0; i < n_edge_local; i++) {
	g->L2Gmap_edge[edge_local[i]] = local_start + i;
    }

    for (i = 0; i < count; i++) {
	if (edge_mark[i] == 0) { /* local */
	}
	/* remote */
	else {
	    nn = g->L2Gmap_vert[edge2vert[i][0]]
		+ g->L2Gmap_vert[edge2vert[i][1]];
	    nd = nn / g->nprocs;
	    nm = nn % g->nprocs;

	    rank = nm;
	    assert(bsend[sdsps[rank]][0] == nn); /* check sum */
	    g->L2Gmap_edge[i] = bsend[sdsps[rank]][1]; /* global index */
	    sdsps[rank]++;
	}
    }
    /* restore sdsps[] */
    for (rank = 0; rank < g->nprocs; rank++)
	sdsps[rank] -= scnts[rank];


    g->nedge_global = n_local_all + n_remote_all;
    phgInfo(0, "nedge global: %d\n", g->nedge_global);
}








/* --------------------------------------------------------------------------------
 *
 *
 *   Set face index
 *
 *
 * -------------------------------------------------------------------------------- 
 * */
typedef INT FaceVerts[3];
typedef struct FaceVertsSort_ {
    INT v0;			/* vert */
    INT v1;			/* vert */
    INT v2;			/* vert */
    INT ridx;			/* recv index */
    INT sidx;			/* uniq index */
} FaceVertsSort;


static int
phgCompFace(const void *p0, const void *p1)
/* compares two 'INT's (used in qsort and bsearch) */
{
    int i;
    FaceVerts *vts0 = (FaceVerts *) p0;
    FaceVerts *vts1 = (FaceVerts *) p1;

    i = (*vts0)[0] - (*vts1)[0];
    
    if (i < 0) {
	return -1;
    } else if (i > 0) {
	return 1;
    } else {
	i = (*vts0)[1] - (*vts1)[1];
	if (i < 0) {
	    return -1;
	} else if (i > 0) {
	    return 1;
	} else {
	    i = (*vts0)[2] - (*vts1)[2];
	    return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
	}
    }

}

static int
phgCompFace2(const void *p0, const void *p1)
/* compares two 'INT's (used in qsort and bsearch) */
{
    int i;
    FaceVertsSort *vts0 = (FaceVertsSort *) p0;
    FaceVertsSort *vts1 = (FaceVertsSort *) p1;

    i = (*vts0).v0 - (*vts1).v0;
    
    if (i < 0) {
	return -1;
    } else if (i > 0) {
	return 1;
    } else {
	i = (*vts0).v1 - (*vts1).v1;

	if (i < 0) {
	    return -1;
	} else if (i > 0) {
	    return 1;
	} else {
	    i = (*vts0).v2 - (*vts1).v2;
	    return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
	}
    }

}

#define SORT_VERT3(v0, v1, v2) {		\
	INT vv;					\
	if (v0 > v1) {				\
	    vv = v0; v0 = v1; v1 = vv;		\
	}					\
	if (v0 > v2) {				\
	    vv = v0; v0 = v2; v2 = vv;		\
	}					\
	if (v1 > v2) {				\
	    vv = v1; v1 = v2; v2 = vv;		\
	}					\
	assert(v0 != v1);			\
	assert(v1 != v2);			\
	assert(v2 != v0);			\
    }


/*
 * Set face index
 * 1. set local face index
 * 2. mark remote vert then face
 * 3. alltoall determine global index remote face
 *
 * */
static void
set_face_index(GRID *g) 
{
    INT i, j, k;
    SIMPLEX *e;

    /* ------------------------------------------------------------
     * 
     * Step 1. set local face index
     *
     * ------------------------------------------------------------ */
    /* Note:
     *  phgUpdateFaces(g); is not safe.
     *  Since the local grid is not big, just use slow O(N^2) algorithm
     * */
    FaceVerts *face2vert;
    INT count, size;

    count = 0;
    face2vert = phgCalloc((size = 100), sizeof(FaceVerts));
    ForAllElements(g, e) {
	for (i = 0; i < NFace; i++) {
	    INT v0 = e->verts[GetFaceVertex(i, 0)];
	    INT v1 = e->verts[GetFaceVertex(i, 1)];
	    INT v2 = e->verts[GetFaceVertex(i, 2)];

	    SORT_VERT3(v0, v1, v2);
	    
	    if (count >= size) {
		face2vert = phgRealloc_(face2vert,
					2*size * sizeof(FaceVerts),
					size * sizeof(FaceVerts));
		size *= 2;
	    }
	    face2vert[count][0] = v0;
	    face2vert[count][1] = v1;
	    face2vert[count][2] = v2;
	    count++;
	}
    }
    
    /* remove duplicate */
    qsort(face2vert, count, sizeof(FaceVerts), phgCompFace);

    {
	INT i0 = 0;
	for (i = i0+1; i < count; i++) {
	    int cmp = phgCompFace(face2vert + i0,
				    face2vert + i);
	    if (cmp < 0) {
		i0++;
		if (i != i0)
		    memcpy(face2vert + i0,
			   face2vert + i,
			   sizeof(FaceVerts));
	    }
	}
	count = i0 + 1;
    }
    g->nface = count;
    phgInfo(0, "pgrid nface: %d\n", g->nface);

    ForAllElements(g, e) {
	for (i = 0; i < NFace; i++) {
	    INT v0 = e->verts[GetFaceVertex(i, 0)];
	    INT v1 = e->verts[GetFaceVertex(i, 1)];
	    INT v2 = e->verts[GetFaceVertex(i, 2)];

	    SORT_VERT3(v0, v1, v2);
	    
	    FaceVerts vts, *p;
	    vts[0] = v0;
	    vts[1] = v1;
	    vts[2] = v2;
	    p = bsearch(&vts, face2vert,
			count, sizeof(FaceVerts),
			phgCompFace);
	    assert(p != NULL);

	    e->faces[i] = p - face2vert;
	}
    }
    


    /* ------------------------------------------------------------
     *
     * Step 2. mark remote face
     *
     * ------------------------------------------------------------ */
    BYTE *vert_mark, *face_mark;
    INT *face_local, *face_remote;
    INT n_face_local = 0, n_face_remote = 0;

    face_mark   = phgCalloc(count, sizeof(*face_mark));
    face_local  = phgCalloc(count, sizeof(*face_local));
    face_remote = phgCalloc(count, sizeof(*face_remote));
    g->L2Gmap_face = phgCalloc(count, sizeof(*g->L2Gmap_face));
    
    for (i = 0; i < count; i++)
    	g->L2Gmap_face[i] = -1;

    ForAllElements(g, e) {
    	for (i = 0; i < NFace; i++) {
    	    if (e->bound_type[i] & REMOTE) {
		face_mark[e->faces[i]] = 1;	     /* face remote */
    	    }
    	}
    }
    
    for (i = 0; i < count; i++) {
    	if (face_mark[i] == 1) {
    	    face_remote[n_face_remote++] = i;
    	} else {
    	    face_local[n_face_local++] = i;
    	}
    }

    INT rank, nsend, nrecv,
    	*scnts, *sdsps, *rcnts, *rdsps;
    FaceVerts *bsend, *brecv;
    scnts = phgCalloc(4 * g->nprocs, sizeof(*scnts));
    sdsps = scnts + g->nprocs;
    rcnts = sdsps + g->nprocs;
    rdsps = rcnts + g->nprocs;

    INT nd, nm, nn;
    for (i = 0; i < count; i++) {
    	if (face_mark[i] == 1) {
    	    nn = g->L2Gmap_vert[face2vert[i][0]]
    		+ g->L2Gmap_vert[face2vert[i][1]] 
		+ g->L2Gmap_vert[face2vert[i][2]];
    	    nd = nn / g->nprocs;
    	    nm = nn % g->nprocs;

	    /* phgInfo(0, "remote face mark: %6d %6d %6d \n", */
	    /* 	    g->L2Gmap_vert[face2vert[i][0]], */
	    /* 	    g->L2Gmap_vert[face2vert[i][1]], */
	    /* 	    g->L2Gmap_vert[face2vert[i][2]] */
	    /* 	    ); */
    	    rank = nm;
    	    scnts[rank]++;
    	}
    }

    MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, g->comm);

    nsend = nrecv = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
    	sdsps[rank] = nsend;
    	rdsps[rank] = nrecv;
    	nsend += scnts[rank];
    	nrecv += rcnts[rank];
    }
    
    bsend = phgAlloc((nsend + nrecv) * sizeof(*bsend));
    brecv = bsend + nsend;

    for (i = 0; i < count; i++) {
    	if (face_mark[i] == 1) {
    	    nn = g->L2Gmap_vert[face2vert[i][0]]
    		+ g->L2Gmap_vert[face2vert[i][1]]
		+ g->L2Gmap_vert[face2vert[i][2]];
    	    nd = nn / g->nprocs;
    	    nm = nn % g->nprocs;

    	    rank = nm;
    	    INT V0 = g->L2Gmap_vert[face2vert[i][0]];
    	    INT V1 = g->L2Gmap_vert[face2vert[i][1]];
    	    INT V2 = g->L2Gmap_vert[face2vert[i][2]];

	    SORT_VERT3(V0, V1, V2);

    	    bsend[sdsps[rank]][0] = V0;
    	    bsend[sdsps[rank]][1] = V1;
    	    bsend[sdsps[rank]][2] = V2;
    	    sdsps[rank]++;
    	}
    }

    /* restore sdsps[] */
    for (rank = 0; rank < g->nprocs; rank++)
    	sdsps[rank] -= scnts[rank];
    
    /* send lists */
    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(*bsend), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(bsend, scnts, sdsps, type,
    		 brecv, rcnts, rdsps, type, g->comm);

    /* process received lists */
    FaceVertsSort *recv_sort;
    recv_sort = phgCalloc(nrecv, sizeof(*recv_sort));
    for (i = 0; i < nrecv; i++) {
    	recv_sort[i].v0 = brecv[i][0];
    	recv_sort[i].v1 = brecv[i][1];
    	recv_sort[i].v2 = brecv[i][2];
    	recv_sort[i].ridx = i;
    }

    qsort(recv_sort, nrecv, sizeof(*recv_sort), phgCompFace2);

    /* for (i = 0; i < nrecv; i++) { */
    /* 	phgInfo(0, "   remote face sort: %6d %6d %6d \n", */
    /* 		recv_sort[i].v0,  */
    /* 		recv_sort[i].v1, */
    /* 		recv_sort[i].v2); */
    /* } */

    /* count face number */
    INT remote_face_count = 0;
    if (nrecv > 0) {
    	INT i0 = 0;
    	recv_sort[0].sidx = 0;
    	for (i = 1; i < nrecv; i++) {
    	    int cmp = phgCompFace2(recv_sort + i - 1,
    				       recv_sort + i);
    	    if (cmp < 0) {
    		i0++;
    	    }
    	    recv_sort[i].sidx = i0;
    	}
    	remote_face_count = i0 + 1;
    }
    assert(2 * remote_face_count == nrecv); /* each remote face has two dups */
    /* phgInfo(0, "remote face count: %d %d\n", */
    /* 	    remote_face_count, nrecv); */

    /* Reassign local face number */
    INT remote_start = 0, n_remote_all = 0;
    INT local_start = 0, n_local_all = 0;
    {
    	INT n_face_locals[g->nprocs];
    	MPI_Allgather(&n_face_local, 1, PHG_MPI_INT,
    		      &n_face_locals, 1, PHG_MPI_INT, g->comm);
	
    	for (i = 0; i < g->rank; i++)
    	    local_start += n_face_locals[i];
    	for (i = 0; i < g->nprocs; i++)
    	    n_local_all += n_face_locals[i];
    }

    /* Reassign remote face number */
    {
    	INT remote_face_counts[g->nprocs];
    	MPI_Allgather(&remote_face_count, 1, PHG_MPI_INT,
    		      &remote_face_counts, 1, PHG_MPI_INT, g->comm);
	
    	for (i = 0; i < g->rank; i++)
    	    remote_start += remote_face_counts[i];
    	remote_start += n_local_all;
    	for (i = 0; i < g->nprocs; i++)
    	    n_remote_all += remote_face_counts[i];
    }

    for (i = 0; i < nrecv; i++) {
    	INT i0 = recv_sort[i].ridx;
    	brecv[i0][0] = brecv[i0][0] + brecv[i0][1] + brecv[i0][2]; /* sum check */
    	brecv[i0][1] = recv_sort[i].sidx + remote_start;
    }

    /* send back */
    MPI_Alltoallv(brecv, rcnts, rdsps, type,
    		 bsend, scnts, sdsps, type, g->comm);

    
    for (i = 0; i < n_face_local; i++) {
    	g->L2Gmap_face[face_local[i]] = local_start + i;
    }

    for (i = 0; i < count; i++) {
    	if (face_mark[i] == 0) { /* local */
    	}
    	/* remote */
    	else {
    	    nn = g->L2Gmap_vert[face2vert[i][0]]
    		+ g->L2Gmap_vert[face2vert[i][1]]
		+ g->L2Gmap_vert[face2vert[i][2]];
    	    nd = nn / g->nprocs;
    	    nm = nn % g->nprocs;

    	    rank = nm;
    	    assert(bsend[sdsps[rank]][0] == nn); /* check sum */
    	    g->L2Gmap_face[i] = bsend[sdsps[rank]][1]; /* global index */
    	    sdsps[rank]++;
    	}
    }
    /* restore sdsps[] */
    for (rank = 0; rank < g->nprocs; rank++)
    	sdsps[rank] -= scnts[rank];


    g->nface_global = n_local_all + n_remote_all;
    phgInfo(0, "nface global: %d\n", g->nface_global);

}





#if 0
#define GET_HERE(desp)						\
    printf("PHG: get %-15s L%5d: %s\n", __FILE__, __LINE__, desp);	
#else
#define GET_HERE(desp)	
#endif

void phg_init_(int *fComm)
{
    int argc = 1;
    char *argv_[4] = {
	"phg",
	NULL
    };
    char **argv = argv_;
    BOOLEAN debug = FALSE;
    static BOOLEAN phg_initialized = FALSE;

    MPI_Comm comm = MPI_Comm_f2c(*fComm);
    assert(comm == MPI_COMM_WORLD);

    if (!phg_initialized) {
	printf("PHG init\n");
	phgOptionsRegisterNoArg("debug", "mpi debug", &debug);
    
	GET_HERE("func");
	phgInit(&argc, &argv);
	phgOptionsShowUsed();
	GET_HERE("init");

	/* Pauss to attach gdb */
	if (FALSE && debug) {
	    char hostname[256];
	    int _i_ = 0;
	    unsigned int t = 5;
	    int pid = getpid();

	    gethostname(hostname, sizeof(hostname));
	    printf("#### Lengweee debug "
		   "PID %d on %s ready for attach\n", getpid(), hostname);
	    fflush(stdout);
	    while (0 == _i_) {
		MPI_Barrier(MPI_COMM_WORLD);
		printf("%d after bar %d\n", pid, _i_);
		fflush(stdout);
		sleep(100);
		printf("%d after sleep\n", pid);
		fflush(stdout);
	    }
	    printf("### PID %d, ready!\n", getpid());
	}

	phg_initialized = TRUE;
    }	/* end initialized */

}


GRID *
phgImportParallelGrid_(GRID *old_grid,
		      int nvert,
		      int nelem, 
		      int nvert_global,
		      int nelem_global,
		      INT *L2Gmap_vert,
		      INT *L2Gmap_elem,
		      double *coord,
		      TET *tet,
		      MPI_Comm comm
		      )
{
    int i, j, k;
    SIMPLEX *e;
    static int phg_counter = 0;


    /* printf("\n=================\n"); */
    /* printf("PHG Import pgrid times: %d\n", phg_counter++); */
    /* printf("=================\n"); */

    phgInfo(0, "nvert: %5d, nelem: %5d\n",
	    nvert, nelem);
    phgInfo(0, "nvert_global: %5d, nelem_global: %5d\n",
	    nvert_global, nelem_global);
    /* phgInfo(0, "nedge_global: %5d, nface_global: %5d\n", */
    /* 	    g0->nedge_global, g0->nface_global); */


    GRID *g = phgNewGrid(-1); 
    GET_HERE("new grid");

    /* Set MPI */
    MPI_Comm_size (comm, &g->nprocs );
    MPI_Comm_rank (comm, &g->rank );
    MPI_Comm_dup(comm, &g->comm);

    /* --------------------------------------------------------------------------------
     * 
     * Only one proc
     *
     * -------------------------------------------------------------------------------- */
    if (g->nprocs == 1) {
	g->nleaf = nelem;
	g->nvert = nvert;
	g->nedge = 0;
	g->nface = 0;
	g->nelem = nelem;

	g->nvert_global = nvert;
	g->nedge_global = 0;
	g->nface_global = 0;
	g->nelem_global = nelem;

	g->nroot = nelem;
	g->ntree = nelem;

    
	//MEM_DUP(g->verts, coord, nvert); /* Other copy */
	g->verts = phgCalloc(g->nvert, sizeof(*g->verts));
	for (i = 0; i < g->nvert; i++) {
	    g->verts[i][0] = coord[i*3  ];
	    g->verts[i][1] = coord[i*3+1];
	    g->verts[i][2] = coord[i*3+2];
	    /* phgInfo(0, "%d %f %f %f\n", i, */
	    /*        coord[i*3  ], coord[i*3+1], coord[i*3+2]); */
	}
	g->L2Gmap_vert = NULL;
	g->L2Gmap_elem = NULL;
    

	/*
	 * Build elements
	 *
	 *  */
	g->roots = phgNewElements(nelem);
	g->elems = phgCalloc(g->nelem, sizeof(*g->elems));

	e = g->roots;
	for (i = 0; i < nelem; i++, e++) {
	    g->elems[i] = e;
	    e->index = i;
	    e->region_mark = tet[i].region_mark;

	    for (j = 0; j < NVert; j++)
		e->verts[j] = tet[i].verts[j];

	    /* Ordering */
	    {
		INT V0, V1, V2, V3,
		    v0, v1, v2, v3;

		V0 = e->verts[0]; //GlobalVertexP(g, e->verts[0]);
		V1 = e->verts[1]; //GlobalVertexP(g, e->verts[1]);
		V2 = e->verts[2]; //GlobalVertexP(g, e->verts[2]);
		V3 = e->verts[3]; //GlobalVertexP(g, e->verts[3]);
		v0 = v1 = v2 = v3 = 0;
		(V0 > V1) ? v0++ : v1++;
		(V0 > V2) ? v0++ : v2++;
		(V0 > V3) ? v0++ : v3++;
		(V1 > V2) ? v1++ : v2++;
		(V1 > V3) ? v1++ : v3++;
		(V2 > V3) ? v2++ : v3++;
		e->ordering = v0 | (v1 << 2) | (v2 << 4) | (v3 << 6);
	    }

	    /* boundary type, neigh */
	    for (j = 0; j < NFace; j++) {
		if (tet[i].bound_type[j] & INTERIOR) {
		    /* interior */
		    e->bound_type[j] = INTERIOR; 
		    //assert(tet[i].neighbours[j] >= 0);
		    //e->neighbours[j] = g->roots
		    //   + tet[i].neighbours[j];
		} else {
		    /* bdry */
		    e->bound_type[j] =
			tet[i].bound_type[j] & BDRY_MASK;
		    
		    /* assert(tet[i].neighbours[j] == -1); */
		    /* e->neighbours[j] = NULL; */
		}
	    }
	}
	GET_HERE("init elem");

	phgUpdateNeighbours(g);

	/* ================================================================================
	 *
	 *  Build edge & face indices 
	 * 
	 * ================================================================================ */
	//phgUpdateNeighbours(g);

	phgUpdateEdges(g);
	//set_edge_index(g);
	GET_HERE("edge index");

	phgUpdateFaces(g);
	//set_face_index(g);
	GET_HERE("face index");


	phgCheckConformity(g);
	GET_HERE("rneigh");
    }
    /* --------------------------------------------------------------------------------
     * 
     * Multiply procs
     *
     * -------------------------------------------------------------------------------- */
    else {
	g->nleaf = nelem;
	g->nvert = nvert;
	g->nedge = -1;		/* set afterwards */
	g->nface = -1;		/* set afterwards */
	g->nelem = nelem;

	g->nvert_global = nvert_global;
	g->nedge_global = -1;	/* set afterwards */
	g->nface_global = -1;	/* set afterwards */
	g->nelem_global = nelem_global;

	g->nroot = nelem;
	g->ntree = nelem;

    
	//MEM_DUP(g->verts, coord, nvert); /* Other copy */
	g->verts = phgCalloc(g->nvert, sizeof(*g->verts));
	for (i = 0; i < g->nvert; i++) {
	    g->verts[i][0] = coord[i*3  ];
	    g->verts[i][1] = coord[i*3+1];
	    g->verts[i][2] = coord[i*3+2];
	    /* phgInfo(0, "%d %f %f %f\n", i, */
	    /*        coord[i*3  ], coord[i*3+1], coord[i*3+2]); */
	}
	MEM_DUP(g->L2Gmap_vert,  L2Gmap_vert, nvert); /* Other copy */
	MEM_DUP(g->L2Gmap_elem,  L2Gmap_elem, nelem); /* Other copy */

    
    

	/*
	 * Build elements
	 *
	 *  */
	g->roots = phgNewElements(nelem);
	g->elems = phgCalloc(g->nelem, sizeof(*g->elems));
	int num_face_type[3] = {0, 0, 0};

	e = g->roots;
	for (i = 0; i < nelem; i++, e++) {
	    g->elems[i] = e;
	    e->index = i;
	    e->region_mark = tet[i].region_mark;

	    for (j = 0; j < NVert; j++)
		e->verts[j] = tet[i].verts[j];

	    /* Ordering */
	    {
		INT V0, V1, V2, V3,
		    v0, v1, v2, v3;

		V0 = g->L2Gmap_vert[e->verts[0]]; //GlobalVertexP(g, e->verts[0]);
		V1 = g->L2Gmap_vert[e->verts[1]]; //GlobalVertexP(g, e->verts[1]);
		V2 = g->L2Gmap_vert[e->verts[2]]; //GlobalVertexP(g, e->verts[2]);
		V3 = g->L2Gmap_vert[e->verts[3]]; //GlobalVertexP(g, e->verts[3]);
		v0 = v1 = v2 = v3 = 0;
		(V0 > V1) ? v0++ : v1++;
		(V0 > V2) ? v0++ : v2++;
		(V0 > V3) ? v0++ : v3++;
		(V1 > V2) ? v1++ : v2++;
		(V1 > V3) ? v1++ : v3++;
		(V2 > V3) ? v2++ : v3++;
		e->ordering = v0 | (v1 << 2) | (v2 << 4) | (v3 << 6);
	    }

	    /* boundary type, neigh */
	    for (j = 0; j < NFace; j++) {
		if (tet[i].bound_type[j] & REMOTE) {
		    /* remote */
		    e->bound_type[j] = (INTERIOR | REMOTE);
		    // assert(tet[i].neighbours[j] == -1);
		    //e->neighbours[j] = NULL;

		    if (tet[i].bound_type[j] & BDRY_MASK) { /* interior face mark */
			e->bound_type[j] |= 
			    (tet[i].bound_type[j] & BDRY_MASK);
		    }
		} else if (tet[i].bound_type[j] & INTERIOR) {
		    /* interior */
		    e->bound_type[j] = INTERIOR; 
		    // assert(tet[i].neighbours[j] >= 0);
		    /* e->neighbours[j] = g->roots */
		    /* 	+ tet[i].neighbours[j]; */

		    if (tet[i].bound_type[j] & BDRY_MASK) { /* interior face mark */
			e->bound_type[j] |= 
			    (tet[i].bound_type[j] & BDRY_MASK);
		    }
		} else {
		    /* bdry */
		    e->bound_type[j] = tet[i].bound_type[j];
		    
		    // assert(tet[i].neighbours[j] == -1);
		    // e->neighbours[j] = NULL;
		}
	    }

	}
	GET_HERE("init elem");

	phgUpdateNeighbours(g);

	{
	    int b[3];
	    memcpy(b, num_face_type, sizeof(num_face_type));
	    MPI_Allreduce(b, num_face_type, 3, PHG_MPI_INT, MPI_SUM, g->comm);
	    phgInfo(0, "Face types: %d %d %d\n", num_face_type[0],
		    num_face_type[1], num_face_type[2]);
	}


	/* ================================================================================
	 *
	 *  Build edge & face indices 
	 * 
	 * ================================================================================ */

	set_edge_index(g);
	if ((g->flags & EDGE_FLAG)) {
	    phgInfo(0, "Rebuild edge flag\n");
	    /* Warning: set index to global */
	    ForAllElements(g, e)
		for (i = 0; i < NEdge; i++)
		    e->edges[i] = g->L2Gmap_edge[e->edges[i]];
	    phgRebuildL2Gmap_(g, EDGE);
	}
	GET_HERE("edge index");

	set_face_index(g);
	if ((g->flags & FACE_FLAG)) {
	    phgInfo(0, "Rebuild face flag\n");
	    /* Warning: set index to global */
	    ForAllElements(g, e)
		for (i = 0; i < NFace; i++)
		    e->faces[i] = g->L2Gmap_face[e->faces[i]];
	    phgRebuildL2Gmap_(g, FACE);
	}
	GET_HERE("face index");

	phgInfo(0, "Update RNeigh\n");
	phgUpdateRNeighbours_(g);
	GET_HERE("rneigh");
    }


    /* ================================================================================
     *
     *  Update Boundary type
     * 
     * ================================================================================ */
    phgInfo(0, "Update bdry types\n");
    phgUpdateBoundaryTypes(g);
    GET_HERE("bdry type");


    /* --------------------------------------------------------------------------------
     *
     *    Rest 
     *
     * --------------------------------------------------------------------------------
     * */
    // phgExportMedit(g, "imported.mesh");


    /* Check */
    ForAllElements(g, e) {
	for (i = 0; i < NEdge; i++)
	    assert(e->edges[i] >= 0
		   && e->edges[i] < g->nedge);
	for (i = 0; i < NFace; i++)
	    assert(e->faces[i] >= 0
		   && e->faces[i] < g->nface);
    }


    phgGeomInit_(g, TRUE);	
    GET_HERE("geom init");


    //test_poisson(g);
    GET_HERE("test poisson");


    /* check bdry type */
    /* { */
    /* 	int num_face_type[3] = {0, 0, 0}; */
    /* 	ForAllElements(g, e) { */
    /* 	    for (j = 0; j < NFace; j++) { */
    /* 		if (e->bound_type[j] & BC_BOTTOM) { */
    /* 		    num_face_type[0]++; */
    /* 		} else if (e->bound_type[j] & BC_TOP) { */
    /* 		    num_face_type[1]++; */
    /* 		} else if (e->bound_type[j] & BDRY_MASK) { */
    /* 		    num_face_type[2]++; */
    /* 		} */
    /* 	    } */
    /* 	} */

    /* 	{ */
    /* 	    int b[3]; */
    /* 	    memcpy(b, num_face_type, sizeof(num_face_type)); */
    /* 	    MPI_Allreduce(b, num_face_type, 3, PHG_MPI_INT, MPI_SUM, g->comm); */
    /* 	    phgInfo(0, "Face types: %d %d %d\n", num_face_type[0], */
    /* 		    num_face_type[1], num_face_type[2]); */
    /* 	} */
    /* } */



    MPI_Barrier(comm);

    /* /\* */
    /*  * Free old grid. */
    /*  * */
    /*  *  *\/ */
    /* GRID *g_old = (GRID *) old_grid; */
    /* if (old_grid != NULL) { */
    /* 	printf("PHG Free old grid\n"); */
    /* 	//phgDumpGrid(g_old); */
    /* 	phgFreeGrid(&g_old); */
    /* } */

    return (GRID *)g;
}

typedef struct {
    INT		index;	/* global index */ 
    INT		lindex;	/* local index */ 
    SHORT	rank;	/* process owning the entry */
    FLOAT       coord[3];
} B_LIST;


static B_LIST *blist;

typedef struct {
    INT		index;	/* global index */ 
    INT		lindex;	/* local index */ 
    INT         Gindex;	/* Global index */
    FLOAT       coord[3];
} B_LIST2;


static int
blist_comp_indices(const void *p0, const void *p1)
{
    int i;
    B_LIST *b0 = blist + *((const INT *)p0), *b1 = blist + *((const INT *)p1);

    return b0->index - b1->index; /* compare global index */
}


/*
 *
 * Get sub grid. 
 *
 *
 *  */
GRID *
phgGetSubGrid_(GRID *g, SUBREGION_FUNC subregion_func)
{
    SIMPLEX *e;
    int i, j;
    INT nvert_owned, nvert_local, nvert_global, *L2Gmap_vert;
    INT nelem_owned, nelem_local, nelem_global, *L2Gmap_elem;
    int *vert_owned_map, *vert_local_map, *vert_father_map, 
	*elem_owned_map;
    int *neighs;
    int nelem, nelems[g->nprocs],
	myrank, nprocs, ranks[g->nprocs];
    NEIGHBOUR_DATA *neigh_data;


    DOF *dof_mark = phgDofNew(g, DOF_DG0, 1, "sub region mark", DofNoAction);
    phgDofSetDataByValue(dof_mark, 0.);
    ForAllElements(g, e) {
	if (subregion_func(g, e)) {
	    dof_mark->data[e->index] = 1.;
	}
    }
    neigh_data = phgDofInitNeighbourData(dof_mark, NULL);      
    neighs = phgCalloc(g->nelem * NFace, sizeof(*neighs));
    ForAllElements(g, e) {
	if (subregion_func(g, e)) {
	    for (i = 0; i < NFace; i++) {
		if (e->bound_type[i] & INTERIOR) {
		    FLOAT mark = *phgDofNeighbourData(neigh_data, e, i, 0, NULL);
		    if (mark > .5)
			neighs[e->index * NFace + i] = 1;
		}
	    }
	}
    }
    phgDofReleaseNeighbourData(&neigh_data);
    phgDofFree(&dof_mark);

    /* 0. Create sub communicator */
    nelem = 0;
    ForAllElements(g, e) {
	if (subregion_func(g, e)) {
	    nelem++;
	}
    }
    MPI_Allgather(&nelem, 1, PHG_MPI_INT, 
		  nelems, 1, PHG_MPI_INT, g->comm);

    nprocs = 0;
    for (i = 0; i < g->nprocs; i++)
	if (nelems[i] > 0) {
	    ranks[nprocs++] = i;
	}


    MPI_Group orig_group, new_group; 
    MPI_Comm new_comm; 
    MPI_Comm_group(g->comm, &orig_group); 
    MPI_Group_incl(orig_group, nprocs, ranks, &new_group);
    MPI_Comm_create(g->comm, new_group, &new_comm); 

    if (nelem == 0)
	return NULL;

    MPI_Comm_rank(new_comm, &myrank); 
    phgInfo(0, "subregion rank: %d -> %d\n", g->rank, myrank);


    /* 1. build nvert, nvert_global, L2Gmap_vert.
     *    vert_local_map[g->nvert]: index of new local vert
     *    vert_owned_map[g->nvert]: index of new owned vert
     * */
    vert_owned_map = phgCalloc(g->nvert, sizeof(*vert_owned_map));
    vert_local_map = phgCalloc(g->nvert, sizeof(*vert_local_map));
    for (i = 0; i < g->nvert; i++) {
	vert_owned_map[i] = -1;
	vert_local_map[i] = -1;
    }

    nvert_local = 0;
    ForAllElements(g, e) {
	if (subregion_func(g, e)) {
	    for (i = 0; i < NVert; i++) {
		assert(e->verts[i] < g->nvert);
		if (vert_local_map[e->verts[i]] == -1) {
		    vert_local_map[e->verts[i]] = nvert_local++;
		}
	    }
	}
    }

    double *coord;
    coord = phgCalloc(nvert_local * Dim, sizeof(*coord));

    ForAllElements(g, e) {
	if (subregion_func(g, e)) {
	    for (i = 0; i < NVert; i++) {
		if (vert_local_map[e->verts[i]] != -1) {
		    INT ivert = vert_local_map[e->verts[i]];
		    coord[3*ivert  ] = g->verts[e->verts[i]][0];
		    coord[3*ivert+1] = g->verts[e->verts[i]][1];
		    coord[3*ivert+2] = g->verts[e->verts[i]][2];
		}
	    }
	}
    }


    /*
     * Build vert ownership,
     * 
     */
    SHORT *owner_rank_vert;
    INT *owner_index_vert;
    BOOLEAN *owner_flag_vert;
    owner_rank_vert = phgAlloc(g->nvert * sizeof(*owner_rank_vert));
    owner_index_vert = phgAlloc(g->nvert * sizeof(*owner_index_vert));
    owner_flag_vert = phgCalloc(g->nvert, sizeof(*owner_flag_vert));

    for (i = 0; i < g->nvert; i++) {
	owner_rank_vert[i] = -1;
	owner_index_vert[i] = -1;
    }

    {
	/* count entries which have REMOTE bit set */
	int nsend, nrecv, *scnts, *sdsps, *rcnts, *rdsps;
	int rank, ii, jj, n, n0;
	scnts = phgCalloc(4 * nprocs, sizeof(*scnts));
	sdsps = scnts + nprocs;
	rcnts = sdsps + nprocs;
	rdsps = rcnts + nprocs;
	for (i = 0; i < g->nvert; i++) {
	    if (vert_local_map[i] != -1 &&
		g->types_vert[i] & REMOTE) {
		j = GlobalVertex(g, i);
		rank = j % nprocs;
		scnts[rank]++;
	    }
	}

	MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, new_comm);

	nsend = nrecv = 0;
	for (rank = 0; rank < nprocs; rank++) {
	    sdsps[rank] = nsend;
	    rdsps[rank] = nrecv;
	    nsend += scnts[rank];
	    nrecv += rcnts[rank];
	}
	phgInfo(3, "nsend = %d, nrecv = %d, ntotal = %d\n", nsend, nrecv,
		g->nvert);
	
	/* collect entries which have REMOTE bit */
	B_LIST *bsend, *brecv, *bl = NULL;
	bsend = phgAlloc((nsend + nrecv) * sizeof(*bsend));
	brecv = bsend + nsend;

	for (i = 0; i < g->nvert; i++) {
	    if (vert_local_map[i] != -1 &&
		g->types_vert[i] & REMOTE) {
		j = GlobalVertex(g, i);
		rank = j % nprocs;
		bl = bsend + (sdsps[rank]++);

		bl->index  = j;	/* father global idx */
		bl->lindex = i;	/* father local index */
		bl->rank   = myrank; /* rank */
		bl->coord[0] = g->verts[i][0];
		bl->coord[1] = g->verts[i][1];
		bl->coord[2] = g->verts[i][2];
	    }
	}
	
	for (ii = 0; ii < nsend; ii++) {
	    bl = bsend + ii;
	    phgInfo(3, "-vert :%5d,%5d rank: %d\n", 
		    bl->lindex, bl->index, bl->rank);
	}

	/* restore sdsps[] */
	for (rank = 0; rank < nprocs; rank++)
	    sdsps[rank] -= scnts[rank];

	/* send lists */
	MPI_Datatype type;
	MPI_Type_contiguous(sizeof(*bsend), MPI_BYTE, &type);
	MPI_Type_commit(&type);
	MPI_Alltoallv(bsend, scnts, sdsps, type, brecv, rcnts, rdsps, type, new_comm);

	/* process received lists */
	INT *rind;		/* sorted indices of brecv[] */
	rind = phgAlloc(nrecv * sizeof(*rind));
	for (ii = 0; ii < nrecv; ii++)
	    rind[ii] = ii;
	blist = brecv;
	qsort(rind, nrecv, sizeof(*rind), blist_comp_indices);
	for (ii = 0; ii < nrecv; ) {
	    for (n = ii + 1;
		 n < nrecv && blist_comp_indices(rind + ii, rind + n) == 0;
		 n++);

#warning Min rank get ownership: Vert and edge
	    /* the entries in the range [ii, n) correspond to the same entity */
	    rank = 9999999;		/* the owner process rank */
	    n0 = -1;		/* the local index in the owner process */
	    for (jj = ii; jj < n; jj++) {
		bl = brecv + rind[jj];
		if (bl->rank < rank) { /* min rank */
		    rank = bl->rank;
		    n0 = bl->lindex;
		}
		if (FALSE) {
		    assert(fabs(bl->coord[0] - brecv[rind[ii]].coord[0]) < 1e-14);
		    assert(fabs(bl->coord[1] - brecv[rind[ii]].coord[1]) < 1e-14);
		    assert(fabs(bl->coord[2] - brecv[rind[ii]].coord[2]) < 1e-14);
		}
	    }
	    if (rank == 9999999)
		phgError(1, "rank=%d, index=%d\n",
			 bl->rank, bl->index);

	    /* Note: at this stage data in the range [ii, jj) are no longer
	     * needed, they are replaced with the owner index, owner rank and
	     * updated BTYPE value which will be sent back to incoming process */
	    for (jj = ii; jj < n; jj++) {
		bl = brecv + rind[jj];
		bl->rank = rank;	/* replace with owner's rank */
		bl->index = n0;	/* replace with owner's father local index */
	    }
	    ii = n;
	}
	phgFree(rind);

	/* send back lists */
	MPI_Alltoallv(brecv, rcnts, rdsps, type, bsend, scnts, sdsps, type, new_comm);
	MPI_Type_free(&type);

	/* process received list */
	for (ii = 0; ii < nsend; ii++) {
	    bl = bsend + ii;
	    if (bl->rank == myrank && bl->index == bl->lindex) 
		owner_flag_vert[bl->lindex] = TRUE;
	    
	    owner_rank_vert[bl->lindex] = bl->rank;
	    owner_index_vert[bl->lindex] = bl->index;
	    phgInfo(3, "+vert :%5d,%5d rank: %2d [%2d]\n", 
		    bl->lindex, bl->index, bl->rank, myrank);
	}

	phgFree(bsend);
	phgFree(scnts);
    }

    nvert_owned = 0;
    ForAllElements(g, e) {
	if (subregion_func(g, e)) {
	    for (i = 0; i < NVert; i++) {
		if (vert_owned_map[e->verts[i]] == -1
		    && ((!(g->types_vert[e->verts[i]] & REMOTE)) /* local */
			|| owner_flag_vert[e->verts[i]]) /* or remote owned */
		    ) {
		    vert_owned_map[e->verts[i]] = nvert_owned++;
		}
	    }
	}
    }

    INT nvert_owned_all[nprocs], vert_start = 0, elem_start;
    for (i = 0; i < g->nvert; i++) {
	if (vert_local_map[i] != -1)
	    phgInfo(3, "   Vert %5d, %3d %3d %3d %3d (%d|%d)\n", i, 
		    vert_owned_map[i], vert_local_map[i], 
		    owner_rank_vert[i], owner_index_vert[i], 
		    owner_flag_vert[i], g->types_vert[i] & REMOTE);
    }


    MPI_Allgather(&nvert_owned, 1, PHG_MPI_INT, 
		  nvert_owned_all, 1, PHG_MPI_INT, new_comm);
    phgInfo(3, "owned vert: %d\n", nvert_owned_all[myrank]);

    nvert_global = 0;
    for (i = 0; i < nprocs; i++) {
	if (myrank == i)
	    vert_start = nvert_global;
	nvert_global += nvert_owned_all[i];
    }
    
    /* build L2G map */
    L2Gmap_vert = phgCalloc(nvert_local, sizeof(*L2Gmap_vert));
    /* owned vert */
    for (i = 0; i < g->nvert; i++) {
	if (vert_owned_map[i] != -1) {	/* owned */
	    assert(vert_local_map[i] != -1);
	    L2Gmap_vert[vert_local_map[i]] =
		vert_start + vert_owned_map[i];
	}
    }

    /* local not owned */
    B_LIST2 *sbuf, *rbuf, *bl;
    INT *scnts, *sdsps, *rcnts, *rdsps,
	nsend, nrecv;
    int rank, index;
    scnts = phgCalloc(4 * nprocs, sizeof(*scnts));
    sdsps = scnts + nprocs;
    rcnts = sdsps + nprocs;
    rdsps = rcnts + nprocs;

    for (i = 0; i < g->nvert; i++) {
	if (vert_local_map[i] != -1 && 
	    vert_owned_map[i] == -1) {
	    rank = owner_rank_vert[i];
	    index = owner_index_vert[i];
	    phgInfo(3, "non owned %d  rank %d index %d\n", 
		    i, rank, index);
	    assert(rank >= 0);
	    scnts[rank]++;
	}
    }

    MPI_Alltoall(scnts, 1, MPI_INT, rcnts, 1, MPI_INT, new_comm);
    nsend = nrecv = 0;
    for (rank = 0; rank < nprocs; rank++) {
	sdsps[rank] = nsend;
	rdsps[rank] = nrecv;
	nsend += scnts[rank];
	nrecv += rcnts[rank];
    }

    sbuf = phgCalloc(nsend + nrecv, sizeof(*sbuf));
    rbuf = sbuf + nsend;
    for (i = 0; i < g->nvert; i++) {
	if (vert_local_map[i] != -1 && 
	    vert_owned_map[i] == -1) {
	    rank = owner_rank_vert[i];
	    index = owner_index_vert[i];
	    bl = sbuf + sdsps[rank]++;

	    bl->lindex = i; 	/* local index */
	    bl->index = index; 	/* remote index */
	    bl->Gindex = -1;
	    phgInfo(3, "Vert1: %5d(%5d) %5d %5d\n",
	    	    bl->lindex, i, bl->index, bl->Gindex);
	}
    }
    for (rank = 0; rank < nprocs; rank++)
	sdsps[rank] -= scnts[rank]; /* restore */

    /* send lists */
    MPI_Datatype type;
    MPI_Type_contiguous(sizeof(*sbuf), MPI_BYTE, &type);
    MPI_Type_commit(&type);
    MPI_Alltoallv(sbuf, scnts, sdsps, type,
		 rbuf, rcnts, rdsps, type, new_comm);

    /* local index -> global index */
    for (i = 0; i < nrecv; i++) {
    	index = rbuf[i].index;
    	//phgInfo(0, "  own %5d %d\n", index, vert_owned_map[index]);
    	assert (vert_owned_map[index] != -1); 	/* owned */
    	rbuf[i].Gindex = vert_start + vert_owned_map[index];
    	phgInfo(3, "  remote %5d -> %5d\n", index, rbuf[i].Gindex);
    }

    MPI_Alltoallv(rbuf, rcnts, rdsps, type,
		 sbuf, scnts, sdsps, type, new_comm);

    for (i = 0; i < g->nvert; i++) {
	if (vert_local_map[i] != -1 && 
	    vert_owned_map[i] == -1) {

	    rank = owner_rank_vert[i];
	    bl = sbuf + sdsps[rank]++;
	    L2Gmap_vert[vert_local_map[i]] 
		= bl->Gindex;

	    phgInfo(3, "Vert3: %5d(%5d) %5d %5d\n", 
		    bl->lindex, i, bl->index, bl->Gindex);
	}
    }


    /* Map subgrid vert to grid vert */
    vert_father_map = phgCalloc(nvert_local, sizeof(*vert_father_map));
    for (i = 0; i < g->nvert; i++)
	if (vert_local_map[i] != -1)
	    vert_father_map[vert_local_map[i]] = i;

    for (i = 0; i < nvert_local; i++) {
	assert(L2Gmap_vert[i] != -1);
	phgInfo(3, "GVERT %4d %4d (%f %f %f)\n", 
		vert_father_map[i], L2Gmap_vert[i], 
		g->verts[vert_father_map[i]][0],
		g->verts[vert_father_map[i]][1],
		g->verts[vert_father_map[i]][2]
		);
    }






    /* 2. build nelem, nelem_global, L2Gmap_elem.
     * */
    MPI_Allgather(&nelem, 1, PHG_MPI_INT, 
		  nelems, 1, PHG_MPI_INT, new_comm);
    nelem_global = 0;
    for (i = 0; i < nprocs; i++) {
	if (myrank == i)
	    elem_start = nelem_global;
	nelem_global += nelems[i];
    }
    L2Gmap_elem = phgCalloc(nelem, sizeof(*L2Gmap_elem));
    for (i = 0; i < nelem; i++)
	L2Gmap_elem[i] = elem_start + i; 

    /* Map of elem index to subgrid elem index */
    elem_owned_map = phgCalloc(g->nelem, sizeof(*elem_owned_map));
    for (i = 0; i < g->nelem; i++)
	elem_owned_map[i] = -1;
    nelem = 0;
    ForAllElements(g, e) {
	if (subregion_func(g, e)) {
	    elem_owned_map[e->index] = nelem++;
	}
    }
    


    /*
     *
     * Connections
     *
     * */

    TET *tet, *t;
    tet = phgCalloc(nelem, sizeof(*tet));
    t = tet;
    ForAllElements(g, e) {
	SIMPLEX *e_neigh;
	if (subregion_func(g, e)) {
	    for (i = 0; i < NVert; i++)
		t->verts[i] = vert_local_map[e->verts[i]];
	    for (i = 0; i < NFace; i++) {

		if (neighs[e->index * NFace + i]) {/* have neigh, interior */
		    if (e->bound_type[i] & REMOTE) {/* remote */
			t->bound_type[i] = (REMOTE | INTERIOR);
			// t->neighbours[i] = -1;
		    }
		    else {
			t->bound_type[i] = INTERIOR;
			// assert((e_neigh = e->neighbours[i]) != NULL);
			// t->neighbours[i] = elem_owned_map[e_neigh->index];
			// assert(elem_owned_map[e_neigh->index] != -1);
		    }
		} else {	/* boundary */
		    t->bound_type[i] = e->bound_type[i] & BDRY_MASK;
		    // t->neighbours[i] = -1;
		}

	    }
	    t->e = e;
	    t++;
	}
    }


    GRID *g_new = phgImportParallelGrid_(NULL, 
			  nvert_local,
			  nelem, 
			  nvert_global,
			  nelem_global,
			  L2Gmap_vert,
			  L2Gmap_elem,
			  coord,
			  tet,
			  new_comm
			  );

    return g_new;
}








/* --------------------------------------------------------------------------------
 *
 *   
 *   check routines
 *   
 *
 * -------------------------------------------------------------------------------- */
static FLOAT a = 1.0;

static void
func_u_test(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    *value = Cos(2. * M_PI * x / 1000) * Cos(2. * M_PI * y / 1000) * Cos(2. * M_PI * z / 1000);
}

static void
func_f_test(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    func_u_test(x, y, z, value);
    *value = 12. * M_PI * M_PI * *value / 1000 / 1000 + a * *value;
}



static void
build_linear_system(SOLVER *solver, DOF *u_h, DOF *f_h)
{
    int i, j;
    GRID *g = u_h->g;
    SIMPLEX *e;

    assert(u_h->dim == 1);
    ForAllElements(g, e) {
	int N = DofGetNBas(u_h, e);	/* number of basises in the element */
	FLOAT A[N][N], rhs[N], buffer[N];
	INT I[N];

	/* compute \int \grad\phi_i \cdot \grad\phi_j making use of symmetry */
	for (i = 0; i < N; i++) {
	    I[i] = phgSolverMapE2L(solver, 0, e, i);
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
	    if (phgDofDirichletBC(u_h, e, i, func_u_test, buffer, rhs+i,
					DOF_PROJ_NONE)) {
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, buffer); 
	    }
	    else {	/* interior node */
		/* right hand side = \int f * phi_i */
		phgQuadDofTimesBas(e, f_h, u_h, i, QUAD_DEFAULT, rhs + i);
		phgSolverAddMatrixEntries(solver, 1, I + i, N, I, A[i]); 
	    }
	}
	phgSolverAddRHSEntries(solver, N, I, rhs);
    }
}



void poisson_test(GRID *g)
{
    DOF *u_h = phgDofNew(g, DOF_DEFAULT, 1, "u_h", DofInterpolation);
    phgDofSetDirichletBoundaryMask(u_h, BDRY_MASK);
    DOF *f_h = phgDofNew(g, DOF_DEFAULT, 1, "f_h",  func_f_test);
    SOLVER *solver = phgSolverCreate(SOLVER_DEFAULT, u_h, NULL);
    FLOAT L2error;

    phgPrintf("  DOF: %d, unknowns: %d, Dirichlet bdry: %d\n",
	      DofGetDataCountGlobal(u_h), solver->rhs->map->nglobal,
	      solver->rhs->map->bdry_nglobal);
    build_linear_system(solver, u_h, f_h);
    phgSolverSolve(solver, TRUE, u_h, NULL);
    phgSolverDestroy(&solver);

    // phgExportVTK(g, "ice.vtk", u_h, NULL);

    DOF *u = phgDofNew(g, DOF_ANALYTIC, 1, "u", func_u_test);
    phgDofAXPY(-1.0, u, &u_h);
    L2error = phgDofNormL2(u_h);
    phgPrintf("L2 error = %0.6le\n",
	      (double)L2error);

}



