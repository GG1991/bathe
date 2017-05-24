#include "bathe.h"

int bth_extforce(Vec *f_ext){

  return 0;

}


int bth_intforce(Vec *f_int, Vec *u){

  int     e;
  int     error;
  int     npe;
  int     ix[NPE*DIM];

  double  vec_e[NPE*DIM];

  Vec     u_local;

  VecSet(*f_int,0.0);
  VecGhostUpdateBegin(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*u,&u_local);

  for(e=0;e<mesh.nelemv;e++){

    error = bth_elem_intforce( e, &npe, ix, vec_e, &u_local );
    if(error){
      PetscPrintf(PETSC_COMM_WORLD,"bth_force.c.\n"); 
      return 1;
    }

    npe = mesh.elemv[e].npe;
    VecSetValues(*f_int,npe*DIM,ix,vec_e,ADD_VALUES);

  }

  VecAssemblyBegin(*f_int);
  VecAssemblyEnd(*f_int);

  return 0;

}
