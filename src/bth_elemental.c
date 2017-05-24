/* Assembly functions for Tangent matrix and residue*/

#include "bathe.h"

int bth_elem_intforce(int e, int *npe, int ix[NPE*DIM], double vec_e[NPE*DIM], Vec * u_local){

  int           i,j,k,d;
  int           gp;
  int           ii;
  int           ngp;
  int           input_i[10];
  int           error;

  double        det;
  double        disp[NPE][DIM];
  double        coor[NPE][DIM];
  double        strain[VOI];
  double        stress[VOI];
  double        B[VOI][NPE*DIM];
  double        der[NPE][DIM];
  double        jac[DIM][DIM];
  double        ijac[DIM][DIM];
  double    *** ode;
  double      * wp;

  *npe = mesh.elemv[e].npe;
  ngp  = mesh.elemv[e].ngp;

  memset(vec_e,0.0,NPE*DIM*sizeof(double));

  fem_calwei(*npe,DIM,&wp);
  fem_calode(*npe,DIM,&ode);

  for(i=0;i<*npe;i++){
    for(d=0;d<DIM;d++){
      ii=mesh.elemv[e].nodel[i]*DIM+d;
      ix[i*DIM+d]=mesh.elemv[e].nodeg[i]*DIM+d;
      coor[i][d]=mesh.node[mesh.elemv[e].nodel[i]].coor[d];
      VecGetValues(*u_local,1,&ii,&disp[i][d]);
    }
  }

  for(gp=0;gp<ngp;gp++){
    
    fem_caljac3(coor,ode,*npe,gp,jac);

    fem_invjac3(jac,ijac,&det);
    if(det<0.0)
      return 1;

    fem_calder3(ijac,*npe,gp,ode,der);

    /* Strain at Gauss point */
    error = bth_strain(*npe,disp,der,B,strain);
    if(error){
      return 1;
    }

    /* Stress at Gauss point */
    input_i[0] = gp;
    error = bth_property( e, input_i, NULL, strain, NULL, stress, CALC_STRESS, NULL, NULL );
    if(error){
      return 1;
    }

    /* Calculating vec_e */
    for(i=0;i<*npe;i++){
      for(j=0;j<DIM;j++){
        for(k=0;k<VOI;k++){
          vec_e[i*DIM+j] += B[k][i*DIM+j]*stress[k] * det * wp[gp];
        }
      }
    }

  }/*gp loop*/ 
  
  return 0;

}


int bth_elem_k(int e, int * npe, int ix[NPE*DIM], double K_e[NPE*DIM], Vec * u_local){

  int           i,j,k,d;
  int           gp;
  int           ii;
  int           ngp;
  int           error;
  int           input_i[10];

  double        det;
  double        disp[NPE][DIM];
  double        coor[NPE][DIM];
  double        strain[VOI];
  double        stress[VOI];
  double        c_t[VOI][VOI];
  double        B[VOI][NPE*DIM];
  double        Baux[VOI][NPE*DIM];
  double        der[NPE][DIM];
  double        jac[DIM][DIM];
  double        ijac[DIM][DIM];
  double    *** ode;
  double      * wp;

  *npe = mesh.elemv[e].npe;
  ngp  = mesh.elemv[e].ngp;

  memset(K_e,0.0,NPE*DIM*NPE*DIM*sizeof(double));
  fem_calwei(*npe,DIM,&wp);
  fem_calode(*npe,DIM,&ode);

  for(i=0;i<*npe;i++){
    for(d=0;d<DIM;d++){
      ii=mesh.elemv[e].nodel[i]*DIM+d;
      ix[i*DIM+d]=mesh.elemv[e].nodeg[i]*DIM+d;
      coor[i][d]=mesh.node[mesh.elemv[e].nodel[i]].coor[d];
      VecGetValues(*u_local,1,&ii,&disp[i][d]);
    }
  }

  for(gp=0;gp<ngp;gp++){
    
    fem_caljac3(coor,ode,*npe,gp,jac);
    fem_invjac3(jac,ijac,&det);

    if(det<0.0)
      return 1;

    fem_calder3(ijac,*npe,gp,ode,der);

    /* Strain at Gauss point */
    error = bth_strain(*npe,disp,der,B,strain);
    if(error){
      return 1;
    }
    
    /* Stress at Gauss point */
    input_i[0] = gp;
    error = bth_property( e, input_i, NULL, strain, c_t, stress, CALC_TANGEN, NULL, NULL );
    if(error){
      return 1;
    }

    /* Build K_e */

    for(i=0;i<VOI;i++){
      for(j=0;j<(*npe)*DIM;j++){
        Baux[i][j]=0.0;
        for(k=0;k<VOI;k++){
          Baux[i][j] += c_t[i][k]*B[k][j];
        }
      }
    }

    for(i=0;i<(*npe)*DIM;i++){
      for(j=0;j<(*npe)*DIM;j++){
        for(k=0;k<VOI;k++){
          K_e[i*((*npe)*DIM)+j] += B[k][i]*Baux[k][j] * det * wp[gp];
        }
      }
    }

  }/*gp loop*/ 
  
  return 0;

}
