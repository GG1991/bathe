/* evolutes material properties on elements */

#include "bathe.h"

int bth_evolute(int step, int flag_conv, Vec * u){

  int           e,i,d,ii,gp;
  int           flag;
  int           error;
  int           input_i[10];
  int           ngp, npe;

  double        strain[VOI];
  double        disp[NPE][DIM];
  double        coor[NPE][DIM];
  double        B[VOI][NPE*DIM];
  double        der[NPE][DIM];
  double        jac[DIM][DIM];
  double        ijac[DIM][DIM];
  double        det;
  double    *** ode;
  double      * wp;

  Vec           u_local;
  
  VecGhostUpdateBegin(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*u,&u_local);

  if(step==1){
    flag = INITIALIZE;
  }else{
    flag = (flag_conv == CONV_ON) ? CALC_EVO_FI : CALC_EVO_CO;
  }

  for(e=0;e<mesh.nelemv;e++){

    ngp  = mesh.elemv[e].ngp;
    npe  = mesh.elemv[e].npe;

    for(i=0;i<npe;i++){
      for(d=0;d<DIM;d++){
        ii = mesh.elemv[e].nodel[i]*DIM+d;
        VecGetValues(u_local,1,&ii,&disp[i][d]);
        coor[i][d]=mesh.node[mesh.elemv[e].nodel[i]].coor[d];
      }
    }

    fem_calwei(npe,DIM,&wp);
    fem_calode(npe,DIM,&ode);

    for(gp=0;gp<ngp;gp++){

      fem_caljac3(coor,ode,npe,gp,jac);

      fem_invjac3(jac,ijac,&det);

      fem_calder3(ijac,npe,gp,ode,der);

      /* Strain at Gauss point */
      bth_strain(npe,disp,der,B,strain);

      
      input_i[0] = gp;
      error    = bth_property( e, input_i, NULL, strain, NULL, NULL, flag, NULL, NULL);

    }

  } /* e loop */

  return 0;

}
