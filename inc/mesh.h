/* Mesh declarations
 * uses gmsh structures and linked list.
 * 
 */

#include "list.h"
#include <math.h>
#include <gsl/gsl_linalg.h>

#ifndef MESH_H
#define MESH_H

typedef struct _elem_t{                             
                                            
    int  npe;                         
    int  ngp;                         
    int  *nodel;
    int  *nodeg;
    void *prop;

}elem_t;

typedef struct _node_t{                             
                                            
    double coor[3];
    list_t elemvL;
    list_t elemsL;
    
}node_t;

typedef struct _mesh_t{                 

    int nelemv;                            
    int nelems;
    int nnodes;
    int nghost;
    node_t *node;
    elem_t *elemv;
    elem_t *elems;

}mesh_t;

typedef int (*cpyelem_t) (node_list_t *elem_nl, elem_t *elem);
typedef int (*cpynode_t) (node_list_t *node_nl, node_t *node);

int mesh_alloc(list_t *list_nodes, list_t *list_ghost, cpynode_t cpynode, list_t *list_elemv, cpyelem_t cpyelemv, list_t *list_elems, cpyelem_t cpyelems, mesh_t *mesh);
int mesh_renum(mesh_t *mesh, int *loc2gold, int *loc2gnew);
int mesh_neigh(mesh_t *mesh, int *loc2gnew);
int mesh_vnorm(double *vec, int n, double *mod);
int mesh_calcarea(mesh_t *mesh,elem_t *elem,int dim,double *area);
int mesh_calcarea_ev(mesh_t *mesh,int e,int dim,double *area);
int mesh_calcvolu(mesh_t *mesh,int e,int dim,double *volu);
int mesh_calcvolu_T4(double v1[3],double v2[3],double v3[3],double v4[3],double *vol);
int mesh_calcarea_T3(double v1[3],double v2[3],double v3[3],double *area);
int mesh_calcarea_Q4(double v1[3],double v2[3],double v3[3],double v4[3],double *area);
int mesh_vcross(double *v1, double *v2, double *vr);
int elem_cmp(void *a, void *b);

#endif
