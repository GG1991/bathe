/*  Copy properties functions for elemv and elems
 * 
 *  These functions tells how the list_nodes, list_ghost, list_elemv and list_elems should be copied inside the mesh_t structure.
 *  
 *  They use the global lists : list_t list_physe 
 *                              list_t list_mater
 *                              list_t list_bound 
 * 
 */

#include "funct.h"

int cpynode (node_list_t *node_nl, node_t *node){

  int d;

  if(!node || !node_nl){
    return 1;
  }

  for(d=0;d<3;d++){
    node->coor[d]=((gmshN_t *)(node_nl->data))->coor[d];
  }

  return 0;    

}


int cpyelemv (node_list_t *elem_nl, elem_t *elemv){

  int           gmshid;
  int           d;

  char          name[32];

  node_list_t * onode;

  pv_t        * pv;

  pvl_t       * pvl;

  if(!elem_nl || !elemv){
    return 1;
  }

  elemv->npe   = ((gmshE_t *)(elem_nl->data))->npe;
  elemv->ngp   = ((gmshE_t *)(elem_nl->data))->npe;
  elemv->nodel = (int *)calloc(elemv->npe,sizeof(int));
  elemv->nodeg = (int *)calloc(elemv->npe,sizeof(int));

  for(d=0;d<elemv->npe;d++)
    elemv->nodeg[d]=((gmshE_t *)(elem_nl->data))->node[d]-1;

  onode = list_physe.head;
  gmshid = ((gmshE_t *)(elem_nl->data))->gmshid;
  while(gmshid != ((gmshP_t*)(onode->data))->gmshid){
    onode=onode->next;
    if(!onode){
      return 1;
    }
  }

  strcpy(name,((gmshP_t*)(onode->data))->name);
  onode = list_mater.head;
  while(strcmp(name,((pvl_t*)(onode->data))->name) != 0){
    onode=onode->next;
    if(!onode)
    {
      PetscPrintf(PETSC_COMM_WORLD,"lst2msh.c:phys %s not found in materials.\n",name); 
      return 1;
    }
  }

  pvl = (pvl_t*)onode->data;

  elemv->prop=(pv_t*)calloc(1,sizeof(pv_t));
  pv=(pv_t*)elemv->prop;
  pv->mattyp = pvl->mattyp;
  pv->params = (double*)calloc(pvl->nparvar,sizeof(double));

  return 0;

}

int cpyelems (node_list_t *elem_nl, elem_t *elems){

  int d;
  ps_t *ps;

  if(!elem_nl || !elems)
    return 1;
  elems->npe = ((gmshE_t *)(elem_nl->data))->npe;
  elems->ngp = ((gmshE_t *)(elem_nl->data))->npe;
  elems->nodel=(int *)calloc(elems->npe,sizeof(int));
  elems->nodeg=(int *)calloc(elems->npe,sizeof(int));
  for(d=0;d<elems->npe;d++)
    elems->nodeg[d]=((gmshE_t *)(elem_nl->data))->node[d]-1;
  elems->prop=(ps_t*)calloc(1,sizeof(ps_t));
  ps = (ps_t*)elems->prop;
  ps->gmshid = ((gmshE_t*)elem_nl->data)->gmshid;
  ps->elemv = ((gmshE_t*)elem_nl->data)->elemv;

  return 0; 
}
