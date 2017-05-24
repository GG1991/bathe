/*  Routines for applying boundary conditions
 * 
 */

#include "bathe.h"

int bth_set_boundary(void){

    /* For every boundary on the list list_bound searchs for 
     * the physical entity associated with the same name and assigns 
     * the dim value for each bound.
     * For each boundary search for the nodes that shears the elements 
     * that have the same physical entity associated (gmshid). The nodes 
     * are saved on the list "nl" and are delete in the last part those
     * who are repeated
     *
     * Generates a vector dir_index with the dirichlet index positions 
     */
    
    int e,n,i,j,k,h,d,p,gmshid,index;
    node_list_t * nb,*np,*nba,*nbb,*nna,*nnb;
    list_t dir_index_list;
    list_t neu_index_list;

    nb=list_bound.head;
    while(nb)
    {
      np=list_physe.head;
      while(np){
        if(!strcmp(((bou_t*)nb->data)->nam ,((gmshP_t*)np->data)->name))
          break;
        np=np->next;
      }
      if(!np){
        PetscPrintf(PETSC_COMM_WORLD,"bound.c:boundary %s has no phys entity.\n",((bou_t*)nb->data)->nam);
        return 1;
      }
      gmshid=((gmshP_t*)np->data)->gmshid;
      ((bou_t*)nb->data)->dim=((gmshP_t*)np->data)->dim;
      for(e=0;e<mesh.nelems;e++){
        if(gmshid==((ps_t*)mesh.elems[e].prop)->gmshid){
          list_insert_se(&((bou_t*)nb->data)->esl,(void*)&e);
          for(n=0;n<mesh.elems[e].npe;n++)
            list_insert_se(&((bou_t*)nb->data)->nl,(void*)&mesh.elems[e].nodeg[n]);
        }
      }
      nb=nb->next;
    }
    
    for(i=list_bound.sizelist;i>0;i--){
      nbb=list_bound.head;
      for(j=0;j<i-1;j++){
        nbb=nbb->next;
      }
      nba=list_bound.head;
      for(j=0;j<i-1;j++){
        nna=((bou_t*)nba->data)->nl.head;
        for(h=0;h<((bou_t*)nba->data)->nl.sizelist;h++){
          nnb=((bou_t*)nbb->data)->nl.head;
          for(k=0;k<((bou_t*)nbb->data)->nl.sizelist;k++){
            if(*(int*)nnb->data==*(int*)nna->data){
              if(list_del(&((bou_t*)nbb->data)->nl,nnb)){
                return 1;
              }
              break;
            }
            nnb=nnb->next;
          }
          nna=nna->next;
        }
        nba=nba->next;
      }
    }

    //==============================    
    //   DIRICHLET INDECES
    //==============================    

    list_init(&dir_index_list,sizeof(int),cmp_int);   
    nba=list_bound.head;
    while(nba){
      nna=((bou_t*)nba->data)->nl.head;
      while(nna){
        for(d=0;d<3;d++){
          p=(1<<(2-d));
          index = *(int*)nna->data * DIM + d;
          if( ((((bou_t*)nba->data)->kin)&p)==p )
            list_insertlast(&dir_index_list,(void*)&index);
        }     
        nna=nna->next;
      }
      nba=nba->next;
    }
    n_dir=dir_index_list.sizelist;
    dir_index = (int*)calloc(n_dir,sizeof(int));
    dir_value = (double*)calloc(n_dir,sizeof(double));
    dir_zeros = (double*)calloc(n_dir,sizeof(double));
    nna=dir_index_list.head;

    i=0;
    while(nna){
        dir_index[i] = *(int*)nna->data;
        dir_zeros[i] = 0.0;
        i++;
        nna=nna->next;
    }

    //==============================    
    //   NEWMANN INDECES
    //==============================    

    list_init(&neu_index_list,sizeof(int),cmp_int);   
    nba=list_bound.head;
    while(nba){
      nna=((bou_t*)nba->data)->nl.head;
      while(nna){
        for(d=0;d<3;d++){
          p=(1<<(2-d));
          index = *(int*)nna->data * DIM + d;
          if( ((((bou_t*)nba->data)->kin)&p)==0 )
            list_insertlast(&neu_index_list,(void*)&index);
        }     
        nna=nna->next;
      }
      nba=nba->next;
    }
    n_neu=neu_index_list.sizelist;
    neu_index = (int*)calloc(n_neu,sizeof(int));
    neu_value = (double*)calloc(n_neu,sizeof(double));
    nna=neu_index_list.head;

    i=0;
    while(nna){
        neu_index[i] = *(int*)nna->data;
        i++;
        nna=nna->next;
    }

    return 0;

}

int bth_set_dirichlet(Vec * u){

  int i,d,p;
  node_list_t *nb, *nn;

  i=0;
  nb=list_bound.head;
  while(nb)
  {
    nn=((bou_t*)nb->data)->nl.head;
    while(nn)
    {
      for(d=0;d<3;d++){

        p=1<<(2-d);

        if( ((((bou_t*)nb->data)->kin)&p)==p ){

          switch(d){

            case 0:
              f1d_eval(calcu.t,((bou_t*)nb->data)->fx,&dir_value[i]);
              break;

            case 1:
              f1d_eval(calcu.t,((bou_t*)nb->data)->fy,&dir_value[i]);
              break;

            case 2:
              f1d_eval(calcu.t,((bou_t*)nb->data)->fz,&dir_value[i]);
              break;

            default:
              return 1;

          }
          i++;
        }
      }     
      nn=nn->next;
    }
    nb=nb->next;
  }
  VecSetValues(*u,n_dir,dir_index,dir_value,INSERT_VALUES);
  VecAssemblyBegin(*u);
  VecAssemblyEnd(*u);

  return 0;

}

int bth_set_neumann(Vec * f, double t){

  int i,d,p;
  node_list_t *nb, *nn;

  i=0;
  nb=list_bound.head;
  while(nb)
  {
    nn=((bou_t*)nb->data)->nl.head;
    while(nn)
    {
      for(d=0;d<3;d++){

        p=1<<(2-d);

        if( ((((bou_t*)nb->data)->kin)&p)==0 ){

          switch(d){

            case 0:
              f1d_eval(t,((bou_t*)nb->data)->fx,&neu_value[i]);
              break;

            case 1:
              f1d_eval(t,((bou_t*)nb->data)->fy,&neu_value[i]);
              break;

            case 2:
              f1d_eval(t,((bou_t*)nb->data)->fz,&neu_value[i]);
              break;

            default:
              return 1;

          }
          i++;
        }
      }     
      nn=nn->next;
    }
    nb=nb->next;
  }
  VecSetValues(*f,n_neu,neu_index,neu_value,INSERT_VALUES);
  VecAssemblyBegin(*f);
  VecAssemblyEnd(*f);

  return 0;

}

int cmp_nod(void *a, void *b){
    if ( *(int*)a > *(int*)b){
        return 1;
    }else if(*(int*)a == *(int*)b){
        return 0;
    }else{
        return -1;
    }
}
