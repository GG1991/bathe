
/* 
   ARC-LENGHT METHOD - CILINDRICAL PATH (Crisfield) 

   Constrain condition: 

   du_i^T * du_i^T = dl^2

   In each time increment lambda is incremented initially "d_lam" that
   can be defined differently at each time lapse with final time "tf"

 */

#include "funct.h"

int bth_arclength_2(void){


  FILE        * fp_perf;
  
  int           tot_step;
  int           sub_step;
  int           error;
  int           its;
  int           max_its;
  int           k_restart;

  double        lambda;
  double        dlen;
  double        r_tol;
  double        d_lam, d_lam_1, d_lam_2;
  double        norm;
  double        dt;
  double        t0;
  double        tf;
  double        ang_1, ang_2;
  double        C1,C2,C3;
  double        p_a2, p_12, p_aa, p_11, p_a1, p_22;

  node_list_t * pn;
  node_list_t * pl;

  Mat           Kt;

  Vec           u;       /* displacement   at time  t                           */
  Vec           f_ext;   /* external force at time  t                           */
  Vec           f_int;   /* internal force at time  t                           */
  Vec           R;       /* internal force at time  t                           */
  Vec           du_o;    /* delta u        at time  t-1  (converged)            */
  Vec           du_a;    /* delta u acumulated at time  t    (not converged)    */
  Vec           du_n;    /* delta u = u_i-1 - u_i-1 (not converged)             */
  Vec           du_n1;   /* delta u = u_i-1 - u_i-1 (not converged)             */
  Vec           du_n2;   /* delta u = u_i-1 - u_i-1 (not converged)             */
  Vec           du_1;    /* du_1 = - Kt^-1 (f_int(u_i-1 - lambda_i-1 * f_ext)   */
  Vec           du_2;    /* du_2 =   Kt^-1 f_ext                                */
  Vec           b;       /* rhs vector to solve for du = Kt^-1 b (include BCs)  */

  KSP           ksp;     /* linear solver context */

  PC            pc;      /* preconditioner context */

  /******************************/
  /* Alloc memory */ 
  
  MatCreateAIJ(PETSC_COMM_WORLD,mesh.nnodes*DIM,mesh.nnodes*DIM,ntot*DIM,ntot*DIM,120,NULL,120,NULL,&Kt);
  VecCreateGhost(PETSC_COMM_WORLD,mesh.nnodes*DIM,ntot*DIM,mesh.nghost*DIM,(PetscInt*)ghost,&u); 
  VecDuplicate(u,&f_ext);
  VecDuplicate(u,&f_int);
  VecDuplicate(u,&R);
  VecDuplicate(u,&du_n);
  VecDuplicate(u,&du_n1);
  VecDuplicate(u,&du_n2);
  VecDuplicate(u,&du_o);
  VecDuplicate(u,&du_a);
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

  dlen      = calcu.arclen_2.dlen;
  r_tol     = calcu.arclen_2.r_tol;
  max_its   = calcu.arclen_2.max_its;
  k_restart = calcu.arclen_2.k_restart;
  
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

  calcu.t  = calcu.t0;
  t0       = calcu.t0;
  tot_step = 1;
  sub_step = 1;

  lambda   = 0.0;

  pn       = calcu.time.head;
  pl       = calcu.arclen_2.dlam.head;
    
  PetscPrintf(PETSC_COMM_WORLD,"\nARC-LENGTH 2\n"); 

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

    if(tot_step > 1){
      // primer paso de tiempo lambda = 0
      d_lam   = *(double *)pl->data;
      lambda += d_lam;
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

    while( norm > r_tol && its < max_its){

      /******************************/
      /* Calculate tangent operator Kt */
      
      if(its % k_restart == 0){
        error = bth_calc_k(&Kt, &u);
      }
      KSPSetOperators(ksp,Kt,Kt);
      
      /******************************/
      /* Solve systems */ 

      if(its == 0 || tot_step == 1){

        /* Resolvemos por NR la primera iteracion
           y todas las iteraciones del primer paso 
           de tiempo 
         */

        KSPSolve(ksp,b,du_n);
        VecCopy(du_n,du_o);

      }else{

        KSPSolve(ksp,b,du_1);
        VecCopy(f_ext,b);
        KSPSolve(ksp,b,du_2);

        /* with du_1 and du_2 we calculate d_lam */
        VecDot(du_2, du_2, &p_22);
        VecDot(du_a, du_2, &p_a2);
        VecDot(du_1, du_2, &p_12);
        VecDot(du_a, du_a, &p_aa);
        VecDot(du_1, du_1, &p_11);
        VecDot(du_a, du_1, &p_a1);
//        VecDot(du_2, du_2, &p_22);
//        VecDot(du_o, du_2, &p_a2);
//        VecDot(du_1, du_2, &p_12);
//        VecDot(du_o, du_o, &p_aa);
//        VecDot(du_1, du_1, &p_11);
//        VecDot(du_o, du_1, &p_a1);

        C1 = p_22;
        C2 = 2*(p_a2 + p_12);
        C3 = p_aa + 2*p_a1 + p_11 - pow(dlen,2);

        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"p_11=%e\n",p_11); 
        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"p_22=%e\n",p_22); 
        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"p_a1=%e\n",p_a1); 
        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"p_a2=%e\n",p_a2); 
        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"p_aa=%e\n",p_aa); 
        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"C1=%e\n",C1); 
        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"C2=%e\n",C2); 
        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"C3=%e\n",C3); 
        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"pow(C2,2) - 4*C1*C3=%e\n",pow(C2,2) - 4*C1*C3); 

        if( pow(C2,2) - 4*C1*C3 >= 0){
          d_lam_1   = (-C2 + sqrt(pow(C2,2) - 4*C1*C3))/(2*C1);
          d_lam_2   = (-C2 - sqrt(pow(C2,2) - 4*C1*C3))/(2*C1);
        }else{
          PetscPrintf(PETSC_COMM_WORLD,"Imaginary roots in arclength_2.\n"); 
          PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"Imaginary roots in arclength_2.\n");
          return 1;
        }

        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"d_lam_1 = %e  & d_lam_2 = %e.\n",d_lam_1,d_lam_2);
        
        VecWAXPY( du_n1,d_lam_1,du_2, du_1);
        VecWAXPY( du_n2,d_lam_2,du_2, du_1);

        VecDot(du_n1, du_o,  &p_a1);
        VecDot(du_n2, du_o,  &p_a2);
        VecDot(du_n1, du_n1, &p_11);
        VecDot(du_n2, du_n2, &p_22);
        VecDot(du_o,  du_o , &p_aa);

        ang_1 = acos(p_a1/sqrt(p_11*p_aa));
        ang_2 = acos(p_a2/sqrt(p_22*p_aa));
//        ang_1 = p_a1/(p_11*p_aa);
//        ang_2 = p_a2/(p_22*p_aa);
        
        if(ang_1 >= 0 && ang_2 <0){
          d_lam = d_lam_1;
          VecCopy(du_n1,du_n);
          PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"ang_1 >= 0 & ang_2 < 0.\n");
        }else if(ang_2 >= 0 && ang_1 <0){
          d_lam = d_lam_2;
          VecCopy(du_n2,du_n);
          PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"ang_2 >= 0 & ang_1 < 0.\n"); 
        }else if(fabs(-C3/C2 - d_lam_1) <= fabs(-C3/C2 - d_lam_2)){
          d_lam = d_lam_1;
          VecCopy(du_n1,du_n);
          PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"|-C3/C2 - d_lam_1| <= |-C3/C2 - d_lam_2|\n"); 
        }else if(fabs(-C3/C2 - d_lam_1) > fabs(-C3/C2 - d_lam_2)){
          d_lam = d_lam_2;
          VecCopy(du_n2,du_n);
          PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"|-C3/C2 - d_lam_1| > |-C3/C2 - d_lam_2|\n"); 
        }else{
          PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"undefined state in arclength_2.\n"); 
          return 1;
        }

        lambda += d_lam; 
        PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"lambda : %e \n",lambda); 

      }

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

    VecCopy(du_a, du_o);

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

