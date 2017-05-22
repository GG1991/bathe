/* Materials laws*/

#include "funct.h"

int bth_property( int e, int input_i[10], double input_d[10], double strain[VOI], double c_t[VOI][VOI],     
    double stress[VOI], int flag, int output_i[10], double output_d[10]){

  /*
     Esta funcion sirve para muchas cosas segun el flag :

     - CALC_STRESS : Calcula el stress a partir del strain 
     - CALC_TANGEN : Calcula el tensor constitutivo tangente
     - CALC_EVO_CO : Evoluciona las propiedades fisicas del material segun u
     - CALC_EVO_FI : Establece las propiedades finales una vez convergido 
     - GIVE_PARAMS : Devuelve los parametros del material que son importantes
  */

  int error;

  switch(((pv_t*)mesh.elemv[e].prop)->mattyp){

    case 3:
      error = matlaw_003( e, input_i, input_d, strain, c_t, stress, flag, output_i, output_d);
      break;

    default:
      return 1;

  }

  return (error)?1:0;
}

int matlaw_003( int e, int input_i[10], double input_d[10], double strain[VOI], double c_t[VOI][VOI],     
    double stress[VOI], int flag, int output_i[10], double output_d[10]){
  
  /* Isotropic Elastic material with isotropic damage model

     Esta funcion sirve para muchas cosas segun el flag :

     - INITIALICE  : Calcula el stress a partir del strain 
     - CALC_STRESS : Calcula el stress a partir del strain 
     - CALC_TANGEN : Calcula el tensor constitutivo tangente
     - CALC_EVO_CO : Evoluciona las propiedades fisicas del material segun u
     - CALC_EVO_FI : Establece las propiedades finales una vez convergido 
     - GIVE_PARAMS : Devuelve los parametros del material que son importantes
  */

  int           i,j,k;
  int           xp;
  int           mode;
 
  double        la,mu;
  double        E,v,ft,Gf,len;
  double        d_0,d_1;
  double        r,r0;
  double        tau, A;
  double        sig_0[6];
  double        sig_1[6];
  double        str_1[6];
  double        dstra[6];
  double        C0[6][6];

  pv_t        * pv;
  pvl_t       * p_mat;
  node_list_t * pn;

  pv   = (pv_t *)mesh.elemv[e].prop;

  /*******************************/
  /* buscamos el material en la lista */

  pn   = list_mater.head;
  while(pn){
    if(((pvl_t*)pn->data)->mattyp==pv->mattyp){
      break;
    }
    pn=pn->next;
  }
  if(!pn){
    return 1;
  }
  p_mat = (pvl_t*)pn->data;

  /*******************************/ 

  E   = p_mat->params[0];
  v   = p_mat->params[1];
  ft  = p_mat->params[2]; 
  Gf  = p_mat->params[3]; 
  len = p_mat->params[4]; 
  r0  = ft/sqrt(E);

  la  = (E*v)/((1+v)*(1-2*v));
  mu  = E/(2*(1+v));

  C0[0][0]=la+2*mu ;C0[0][1]=la      ;C0[0][2]=la      ;C0[0][3]=0.0; C0[0][4]=0.0; C0[0][5]=0.0;
  C0[1][0]=la      ;C0[1][1]=la+2*mu ;C0[1][2]=la      ;C0[1][3]=0.0; C0[1][4]=0.0; C0[1][5]=0.0;
  C0[2][0]=la      ;C0[2][1]=la      ;C0[2][2]=la+2*mu ;C0[2][3]=0.0; C0[2][4]=0.0; C0[2][5]=0.0;
  C0[3][0]=0.0     ;C0[3][1]=0.0     ;C0[3][2]=0.0     ;C0[3][3]=mu ; C0[3][4]=0.0; C0[3][5]=0.0;
  C0[4][0]=0.0     ;C0[4][1]=0.0     ;C0[4][2]=0.0     ;C0[4][3]=0.0; C0[4][4]=mu ; C0[4][5]=0.0;
  C0[5][0]=0.0     ;C0[5][1]=0.0     ;C0[5][2]=0.0     ;C0[5][3]=0.0; C0[5][4]=0.0; C0[5][5]=mu ;

  if(flag == INITIALIZE){

    xp  = input_i[0];

    pv->params[1*NGP + xp] = r0;
    pv->params[0*NGP + xp] = 0.0; /* d = 0.0 */

  }else if(flag == CALC_STRESS){

    /* 
       Calculates the stress in xp 

       input_i[0] = xp;
    
    */ 

    xp  = input_i[0];

    d_0 = pv->params[0*NGP + xp];

    for(i=0;i<6;i++){
      stress[i]=0.0;
      for(j=0;j<6;j++){
        stress[i] += (1-d_0)*C0[i][j]*strain[j];
      }
    }
    return 0;

  }else if( flag == CALC_TANGEN ){

    /* 
       Calculates the tangent constitutive tensor in xp 

       input_i[0] = xp;
    */ 
    
    xp   = input_i[0];
    mode = PERTUR;
    //    mode = SECANT;

    d_0  = pv->params[0*NGP + xp];

    if(mode == PERTUR){

      r    = pv->params[1*NGP + xp];

      for(i=0;i<VOI;i++){

        /* We compute every column "i" of c_t separately */

        memset(dstra,0.0,6*sizeof(double));
        memset(sig_0,0.0,6*sizeof(double));
        memset(sig_1,0.0,6*sizeof(double));

        dstra[i] = 1.0e-12 ;

        for(j=0;j<6;j++){
          str_1[j] = strain[j] + dstra[j];
        }

        for(j=0;j<6;j++){
          for(k=0;k<6;k++){
            sig_1[j] += C0[j][k] * str_1[k];
          }
        }

        tau = 0.0;
        for(j=0;j<VOI;j++){
          tau += str_1[j] * sig_1[j];
        }
        tau = sqrt(tau);

        r   = GSL_MAX_DBL(tau,r);

        A   = 1/(Gf*E/(len*pow(ft,2))-0.5);
        if(A<=0.0){
          PetscPrintf(PETSC_COMM_WORLD,"bthmodels.c:elem length too big model 002.\n"); 
          return 1;
        }

        d_1 = 1-r0/r*exp(A*(1-r/r0));

        memset(sig_1,0.0,6*sizeof(double));

        for(j=0;j<6;j++){
          for(k=0;k<6;k++){
            sig_0[j] += (1-d_0) * C0[j][k] * strain[k];
            sig_1[j] += (1-d_1) * C0[j][k] * str_1[k];
          }
        }

        for(j=0;j<6;j++){
          c_t[j][i] =  (sig_1[j] - sig_0[j]) / dstra[i];
        }

      }

    }else if(mode == SECANT){

      for(i=0;i<6;i++){
        for(j=0;j<6;j++){
          c_t[i][j] =  (1-d_0) * C0[i][j];
        }
      }

    }
    return 0;

  }else if( flag == CALC_EVO_CO || flag == CALC_EVO_FI ){

    /* 
       Evolute material variables on every Gauss point  

       input_i[0] = xp
    */ 

    xp   = input_i[0];

    /* undamaged stress at Gauss point */
    for(i=0;i<6;i++){
      sig_0[i]=0.0;
      for(k=0;k<6;k++){
        sig_0[i] += C0[i][k] * strain[k];
      }
    }
    A   = 1/(Gf*E/(len*pow(ft,2))-0.5);
    if(A<=0.0){
      PetscPrintf(PETSC_COMM_WORLD,"bthmodels.c:elem length too big model 002.\n"); 
      return 1;
    }

    r  = pv->params[1*NGP + xp];

    tau = 0.0;
    for(j=0;j<VOI;j++){
      tau += strain[j] * sig_0[j];
    }

    tau = sqrt(tau);

    r   = GSL_MAX_DBL(tau,r);
    d_0 = 1-r0/r*exp(A*(1-r/r0));

    pv->params[0*NGP + xp] = d_0;
    if( flag == CALC_EVO_FI ){
      pv->params[1*NGP + xp] = r;
    }
    return 0;

  }else if( flag == GIVE_PARAMS ){

    /* 
       Returns parameters of the material

       input_i[0]  = xp;

       output_d[0] = tau;
       output_d[1] = r;
       output_d[2] = d;
    */ 
    xp   = input_i[0];

    /* undamaged stress at Gauss point */
    for(i=0;i<6;i++){
      sig_0[i]=0.0;
      for(k=0;k<6;k++){
        sig_0[i] += C0[i][k] * strain[k];
      }
    }

    tau = 0.0;
    for(j=0;j<VOI;j++){
      tau += strain[j] * sig_0[j];
    }
    tau = sqrt(tau);

    output_d[0] = tau;
    output_d[1] = pv->params[1*NGP + xp]; /* r */
    output_d[2] = pv->params[0*NGP + xp]; /* d */
    return 0;
  }

  return 1;

}
