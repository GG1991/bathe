
/* Rutines for assembly K matrix */

#include "bathe.h"

int bth_calc_k(Mat * K, Vec *u){

  int     e;
  int     error;
  int     npe;
  int     ix[NPE*DIM];
  double  K_e[NPE*DIM*NPE*DIM];
  Vec     u_local;

  MatZeroEntries(*K);
  VecGhostUpdateBegin(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*u,&u_local);

  for(e=0;e<mesh.nelemv;e++){

    error = bth_elem_k( e, &npe, ix, K_e, &u_local);
    if(error){
      PetscPrintf(PETSC_COMM_WORLD,"bth_calc_k.c.\n"); 
      return 1;
    }

    npe = mesh.elemv[e].npe;
    MatSetValues(*K,npe*DIM,ix,npe*DIM,ix,K_e,ADD_VALUES);

  }

  MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY);

  /* Dirichlet BC */

  MatAssemblyBegin(*K,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(*K,MAT_FINAL_ASSEMBLY);

  return 0;

}

