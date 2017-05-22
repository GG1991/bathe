/* Gmsh declarations
 * 
 * 
 */

#ifndef GMSH_H
#define GMSH_H

typedef struct _gmshN_t{
    
    int    n;
    double coor[3];
    
}gmshN_t;

typedef struct _gmshE_t{
    
    int    node[8];
    int    npe;
    int    gmshid;
    int    elemv;
    
}gmshE_t;

typedef struct _gmshP_t{
    
    char   name[32];
    int    dim;
    int    gmshid;
    
}gmshP_t;


int gmsh_read(char *file,char *efile,char *nfile,int rank,int dim, list_t *list_nodes, list_t *list_ghost, list_t
*list_elemv, list_t *list_elems, list_t *list_phyce,int **loc2gold,int **loc2gnew,int **npp,int nproc);
int gmsh_readnodes(char *file,int nproc,char *nfile,int rank,list_t *list_nodes);
int gmsh_detghosts(list_t *list_nodes,list_t *list_elemv,list_t *list_ghost);
int gmsh_readghosts(char *file,list_t *list_ghost);
int gmsh_readelemv(char *file,int nproc,char *nfile,int rank,int dim,list_t *list_elemv);
int gmsh_readelems(char *file,int dim,list_t *list_elemv,list_t *list_elems);
int gmsh_readphys(char *file,list_t *list_elemv,list_t *list_elems,list_t *list_phyce);
int gmsh_reenumerate(char *nfile,int rank,list_t *list_nodes,list_t *list_ghost,int **loc2gold,int **loc2gnew,int **npp,int nproc);
int gmsh_elems_belongs(list_t *list_elemv,gmshE_t *elems,int *nelemv);
int gmsh_phys_belongs(list_t *list_elemv,list_t *list_elems,gmshP_t *phys);
int gmsh_isvol(int code,int dim);
int gmsh_npe(int code);
int gmsh_nodcmp(void *a, void *b);

#endif
