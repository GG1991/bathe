
/* 
   ARC-LENGHT METHOD - DISSIPATIVE AND NON DISIPATIVE CONSTRAIN 

   Constrain condition: 
   
   

   In each time increment lambda is incremented initially "d_lam" that
   makes f_ext to produce an external work of "d_work"

 */

#include "funct.h"

int bth_arclength_3(void){


  FILE        * fp_perf;
  
  int           tot_step;
  int           sub_step;
  int           error;
  int           its;
  int           max_its;
  int           k_restart;
  int           internal_energy;
  int           flag;

  double        lambda, lambda_o;
  double        d_work, d_work_u, d_work_d;
  double        r_tol;
  double        d_lam;
  double        norm;
  double        dt;
  double        t0;
  double        tf;
  double        p_f1, p_f2, p_fa, p_fu, p_fu_o;
  double        a1 = 0.01;
  double        a2 = 0.01;

  node_list_t * pn;
  node_list_t * pl;

  Mat           Kt;

  Vec           u;       /* displacement   at time  t                           */
  Vec           u_o;     /* displacement   at time  t-1                         */
  Vec           f_ext;   /* external force at time  t                           */
  Vec           f_int;   /* internal force at time  t                           */
  Vec           R;       /* internal force at time  t                           */
  Vec           du_o;    /* delta u        at time  t-1  (converged)            */
  Vec           du_a;    /* delta u acumulated at time  t    (not converged)    */
  Vec           du_n;    /* delta u = u_i-1 - u_i-1 (not converged)             */
  Vec           du_1;    /* du_1 = - Kt^-1 (f_int(u_i-1 - lambda_i-1 * f_ext)   */
  Vec           du_2;    /* du_2 =   Kt^-1 f_ext                                */
  Vec           b;       /* rhs vector to solve for du = Kt^-1 b (include BCs)  */

  KSP           ksp;     /* linear solver context */

  PC            pc;      /* preconditioner context */

  /******************************/
  /* Alloc memory */ 
  
  MatCreateAIJ(PETSC_COMM_WORLD,mesh.nnodes*DIM,mesh.nnodes*DIM,ntot*DIM,ntot*DIM,120,NULL,120,NULL,&Kt);
  VecCreateGhost(PETSC_COMM_WORLD,mesh.nnodes*DIM,ntot*DIM,mesh.nghost*DIM,(PetscInt*)ghost,&u); 
  VecDuplicate(u,&u_o);
  VecDuplicate(u,&f_ext);
  VecDuplicate(u,&f_int);
  VecDuplicate(u,&R);
  VecDuplicate(u,&du_o);
  VecDuplicate(u,&du_a);
  VecDuplicate(u,&du_n);
  VecDuplicate(u,&du_1);
  VecDuplicate(u,&du_2);
  VecDuplicate(u,&b);
  
  /******************************/
  /* Setting solver options */

  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetType(ksp,KSPGMRES);
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCLU);
  KSPSetFromOptions(ksp);
 
  /******************************/
  /* init variables */

  d_work    = calcu.arclen_3.d_work;
  r_tol     = calcu.arclen_3.r_tol;
  max_its   = calcu.arclen_3.max_its;
  k_restart = calcu.arclen_3.k_restart;
  
  fp_perf   = fopen("arclength.dat","w");

  /* set neumann in final values, lambda is going to scale it */

  pn=calcu.time.head;
  while(pn->next){
    pn = pn->next;
  }
  if(pn){
    tf = ((tcontrol_t*)pn->data)->tf;
  }else{
    return 1;
  }

  error = bth_set_neumann(&f_ext, tf);

  calcu.t         = calcu.t0;
  t0              = calcu.t0;
  tot_step        = 1;
  sub_step        = 1;
  internal_energy = 1;
  flag            = 0;
  lambda          = 0.0;

  pn       = calcu.time.head;
  pl       = calcu.arclen_3.dlam.head;
    
  PetscPrintf(PETSC_COMM_WORLD,"\nARC-LENGTH 3\n"); 

  /******************************/
  /* temporal loop  */
  
  while(pn){
    
    PetscPrintf(PETSC_COMM_WORLD,"\nt:%lf s n:%d\n",calcu.t,tot_step); 
    PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"\nt:%lf s n:%d\n",calcu.t,tot_step); 
    
    /* non-linear loop  */

    its  = 0;

    VecSet(du_a,0.0);
    
    error = bth_set_dirichlet(&u);

    /* Evolute materials properties acording to the new u */

    error = bth_evolute( tot_step, CONV_OFF, &u);

    error = bth_intforce(&f_int,&u); 

    if(tot_step == 1){
      // tot_step = 0 lambda = 0
      d_lam   = 0.0;
      lambda += d_lam;
    }else if(tot_step == 2){
      // tot_step = 1 lambda = d_work
      d_lam   = d_work;
      lambda += d_lam;
    }else{
      flag = 1;
    }

    error   = VecWAXPY( R, -lambda, f_ext, f_int); /* R = f_int - lambda * f_ext */

    PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"lambda_0 : %e \n",lambda); 

    /******************************/
    /* Prepare RHS = b to solve */

    VecCopy(R, b);
    VecScale( b, -1.0);
    VecSetValues(b,n_dir,dir_index,dir_zeros,INSERT_VALUES);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    error = VecNorm( b, NORM_2, &norm);
    PetscPrintf(PETSC_COMM_WORLD,"AL1 it : %2d |R|=%e \n",0,norm); 
    PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"AL1 it : %2d |R|=%e \n",0,norm); 
    
    /* acà deberia entrar solo si tot_step > 1*/
    while( (norm > r_tol && its < max_its) || flag ){

      flag = 0;

      /******************************/
      /* Calculate tangent operator Kt */
      
      if(its % k_restart == 0){
        error = bth_calc_k(&Kt, &u);
      }
      KSPSetOperators(ksp,Kt,Kt);
      
      /******************************/
      /* Solve systems */ 

      KSPSolve(ksp,b,du_1);
      VecCopy(f_ext,b);
      KSPSolve(ksp,b,du_2);

      /* with du_1 and du_2 we calculate d_lam */
      VecDot(f_ext, du_1, &p_f1);
      VecDot(f_ext, du_2, &p_f2);
      VecDot(f_ext, du_a, &p_fa);
      VecDot(f_ext, u   , &p_fu);
      VecDot(f_ext, u_o , &p_fu_o);

      if(tot_step == 2){
        d_lam   = 0.0;
      }else{
        if(internal_energy == 1){
          d_lam   = (2*d_work_u - lambda   * (p_fa + p_f1))/(lambda * p_f2);
        }else{
          d_lam   = (2*d_work_d - lambda_o * (p_fu + p_f1) + lambda * p_fu_o)/(lambda_o * p_f2 - p_fu_o);
        }
      }

      VecWAXPY( du_n,d_lam,du_2, du_1);
      lambda += d_lam; 
      PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"lambda : %e \n",lambda); 


      VecAXPY( du_a, 1.0, du_n);
      VecAXPY( u, 1.0, du_n);

      /******************************/
      /* Evolute materials properties acording to the new u */

      error = bth_evolute( tot_step, CONV_OFF, &u);
      
      error = bth_intforce( &f_int, &u); 

      error = VecWAXPY( R, -lambda, f_ext, f_int);

      VecCopy(R, b);
      VecScale( b, -1.0);
      VecSetValues(b,n_dir,dir_index,dir_zeros,INSERT_VALUES);
      VecAssemblyBegin(b);
      VecAssemblyEnd(b);

      error = VecNorm( b, NORM_2, &norm);

      if(error)
        return 1;
     
      its ++;

      PetscPrintf(PETSC_COMM_WORLD,"AL1 it : %2d |R|=%e \n",its,norm); 
      PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"AL1 it : %2d |R|=%e \n",0,norm); 

    }

    if(tot_step == 1 || tot_step == 2){
      
      VecDot(f_ext, u, &p_fu);
      d_work_u = 0.5*lambda*p_fu;

    }else{

      VecDot(f_ext, u, &p_fu);
      VecDot(f_ext, u_o, &p_fu_o);

      if(internal_energy == 1){

        d_work_d = 0.5*(lambda_o * p_fu - lambda * p_fu_o);
        if(d_work_d > a1 * d_work_u){

          internal_energy = 0;
          PetscPrintf(PETSC_COMM_WORLD,"dissipative power significant.\n"); 
          d_work_d = a2 * d_work_u;
        }

      }

    }

    VecCopy(du_a, du_o);
    VecCopy(u, u_o);
    lambda_o = lambda;

    error = bth_evolute( tot_step, CONV_ON, &u);

    error = print_out(&u,&f_int,&f_ext,tot_step);

    dt = ( ((tcontrol_t*)pn->data)->tf - t0 ) / (((tcontrol_t*)pn->data)->st - 1);
    calcu.t += dt;

    if( sub_step == ((tcontrol_t*)pn->data)->st ){
      sub_step  = 0;
      t0    = ((tcontrol_t*)pn->data)->tf;
      pn    = pn->next;
      pl    = pl->next;
    }

    sub_step ++;
    tot_step ++;

  }

  return 0;

}

