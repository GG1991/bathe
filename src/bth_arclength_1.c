
/* 
   ARC-LENGHT METHOD - DISSIPATIVE AND NON DISIPATIVE CONSTRAIN 

   Thid method increase in one degree of freedom the complete set of equations in 
   order to impose the constrain condition. 

   ************************************************************
   Dissipative Constrain Condition: 
   
   phi_d (u,lambda) = 0.5 * ( lambda_n-1 u - lambda u_n-1)*f - d_work_d = 0

   Non Dissipative Constrain condition (internal energy): 
   
   phi_u (u,lambda) = 0.5 * ( lambda u - lambda_n-1 u_n-1)*f - d_work_u = 0

   ************************************************************

   ALGORITHM TAKEN FROM : "A new arc-length control method based on the rates 
   of the internal and dissipated energy"

   Notes: 
   * I think the matrix is not positive defined (proved with -pc_type lu)

 */

#include "bathe.h"

int bth_arclength_1(void){


  FILE        * fp_perf;
  
  int           i;
  int           tot_step;
  int           sub_step;
  int           error;
  int           its;
  int           max_its;
  int           k_restart;
  int           internal_energy;
  int         * index_k;
  int           index_i;
  int           ksp_its;
  int           reason;

  double        lambda, lambda_o;
  double        a1, a2, d_work, d_work_u, d_work_d;
  double        d_work_u_fixed, d_work_d_fixed;
  double        phi_d, phi_u, phi;
  double        energy_u, energy_d;
  double        r_tol;
  double        d_lam;
  double        norm;
  double        dt, t0, tf;
  double        p_fu, p_fu_o;
  double      * v;
  double        w;

  node_list_t * pn;

  Mat           Kt;

  Vec           u;       /* displacement   at time  t                           */
  Vec           u_o;     /* displacement   at time  t-1                         */
  Vec           f_ext;   /* external force at time  t                           */
  Vec           f_int;   /* internal force at time  t                           */
  Vec           R;       /* internal force at time  t                           */
  Vec           du;      /* delta u = u_i-1 - u_i-1 (not converged)             */
  Vec           b;       /* rhs vector to solve for du = Kt^-1 b (include BCs)  */
  Vec           xlocal;  /* local vector to get values from distributed vector  */

  KSP           ksp;     /* linear solver context */

  PetscBool     set;
  PetscBool     flg;

  /******************************/
  /* Alloc memory */ 
  if(rank == nproc-1){ 
    MatCreateAIJ(PETSC_COMM_WORLD,mesh.nnodes*DIM+1,mesh.nnodes*DIM+1,ntot*DIM+1,ntot*DIM+1,120,NULL,120,NULL,&Kt);
    VecCreateGhost(PETSC_COMM_WORLD,mesh.nnodes*DIM+1,ntot*DIM+1,mesh.nghost*DIM,(PetscInt*)ghost,&u); 
  }else{
    MatCreateAIJ(PETSC_COMM_WORLD,mesh.nnodes*DIM,mesh.nnodes*DIM,ntot*DIM+1,ntot*DIM+1,120,NULL,120,NULL,&Kt);
    VecCreateGhost(PETSC_COMM_WORLD,mesh.nnodes*DIM,ntot*DIM+1,mesh.nghost*DIM,(PetscInt*)ghost,&u); 
  }
//  ierr = MatSetOption(Kt,MAT_SYMMETRIC,PETSC_TRUE);CHKERRQ(ierr);
  MatIsSymmetricKnown(Kt,&set,&flg);
  if(set == PETSC_TRUE){
      if(flg == PETSC_TRUE){
      }
  }

  VecDuplicate(u,&u_o);
  VecDuplicate(u,&f_ext);
  VecDuplicate(u,&f_int);
  VecDuplicate(u,&R);
  VecDuplicate(u,&du);
  VecDuplicate(u,&b);

  index_k = (int*)calloc(ntot*DIM,sizeof(int));
  for(i=0;i<mesh.nnodes*DIM;i++){
    index_k[i] = i;
  }
  index_i = ntot*DIM;
  v = (double*)calloc(mesh.nnodes*DIM,sizeof(double));
  
  /******************************/
  /* Setting solver options */

  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetType(ksp,KSPCG);
  KSPSetFromOptions(ksp);
  KSPSetOperators(ksp,Kt,Kt);
 
  /******************************/
  /* init variables */

  if(!calcu.arclen_1)
      return 1;

  d_work    = calcu.arclen_1->d_work;
  r_tol     = calcu.arclen_1->r_tol;
  max_its   = calcu.arclen_1->max_its;
  k_restart = calcu.arclen_1->k_restart;
  a1        = calcu.arclen_1->a1;
  a2        = calcu.arclen_1->a2;
  
  fp_perf   = fopen("arclength_1.dat","w");

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
  pn              = calcu.time.head;
  tot_step        = 1;
  sub_step        = 1;
  internal_energy = 1;
  lambda          = 0.0;
  lambda_o        = 0.0;
  energy_u        = 0.0;
  energy_d        = 0.0;
  d_work_u_fixed  = 0.0;
  d_work_d_fixed  = 0.0;

  VecSet(u,0.0);
  VecSet(u_o,0.0);

    
  PetscPrintf(PETSC_COMM_WORLD,"\nARC-LENGTH 1\n"); 

  /******************************/
  /* temporal loop  */
  
  while(pn){
    
    PetscPrintf(PETSC_COMM_WORLD,"\nt:%lf s n:%d\n",calcu.t,tot_step); 
    
    /* non-linear loop  */

    its  = 0;

    error = bth_set_dirichlet(&u);

    /* Evolute materials properties acording to the new u */
    error = bth_evolute( tot_step, CONV_OFF, &u);

    /* calculate internal forces */
    error = bth_intforce(&f_int,&u); 

    /* Residue R = f_int - lambda * f_ext */
    error = VecWAXPY( R, -lambda, f_ext, f_int); 

    /*
    Prepare RHS = b to solve
    The value b(n) should be correctly select according 
    to the energy stage where the system is 
    */ 
    VecCopy(R, b);
    VecScale( b, -1.0);
    VecSetValues(b,n_dir,dir_index,dir_zeros,INSERT_VALUES);
    VecDot(f_ext, u,   &p_fu);
    VecDot(f_ext, u_o, &p_fu_o);
    phi_u  = 0.5*( lambda   * p_fu - lambda_o * p_fu_o) - d_work_u_fixed;
    phi_d  = 0.5*( lambda_o * p_fu - lambda   * p_fu_o) - d_work_d_fixed;
    if(rank == (nproc-1)){

      if(tot_step == 1){
        phi = 0.0;
      }else if(tot_step == 2 ){
        phi = d_work - lambda;
      }else{
        if(internal_energy == 1){
          phi = phi_u/lambda;
        }else{
          phi = phi_d/lambda_o;
        }
      }
      VecSetValues(b,1,&index_i,&phi,INSERT_VALUES);

    }
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    error = VecNorm( b, NORM_2, &norm);
    PetscPrintf(PETSC_COMM_WORLD,"AL1 it : %2d |R|=%e \n",0,norm); 

    /* acà deberia entrar solo si tot_step > 1*/
    while( norm > r_tol && its < max_its ){

      /******************************/
      /* Calculate tangent operator Kt */
      
      if(its % k_restart == 0){
        error = bth_calc_k(&Kt, &u);
      }

      VecGhostUpdateBegin(f_ext,INSERT_VALUES,SCATTER_FORWARD);
      VecGhostUpdateEnd(f_ext,INSERT_VALUES,SCATTER_FORWARD);
      VecGhostGetLocalForm(f_ext,&xlocal);

      /* Ahora armamos la columna "n" de la Kt y la fila "n" */
      VecGetValues(xlocal,mesh.nnodes*DIM,index_k,v);
      for(i=0;i<mesh.nnodes*DIM;i++){
        v[i] *= -1.0;
      }
      MatSetOption(Kt, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
      MatSetValues(Kt,mesh.nnodes*DIM,index_k,1,&index_i,v,INSERT_VALUES);

      /* v = -f (warning: sign) */
      if(rank == (nproc-1)){

        if(tot_step == 2){
          /* d phi / d u = 0 */
          for(i=0;i<mesh.nnodes*DIM;i++){
            v[i] = 0.0;
          }
        }
        
        MatSetValues(Kt,1,&index_i,mesh.nnodes*DIM,index_k,v,INSERT_VALUES);
        if(tot_step==2){
          /* d phi / d lam = 1 */
          w = 1;
        }else{
          if(internal_energy == 1){
            /* (d phi / d lam) / (-0.5 * lambda)   = - f * u / lambda */
            w = -p_fu/lambda;
          }else{
            /* (d phi / d lam) / (-0.5 * lambda_o) = - f * u_o / lambda_o */
            w =  p_fu_o/lambda_o;
          }
        }
        MatSetValues(Kt,1,&index_i,1,&index_i,&w,INSERT_VALUES);

      }
      MatAssemblyBegin(Kt,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(Kt,MAT_FINAL_ASSEMBLY);

      /* Dirichlet BC */

      MatSetOption(Kt,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
      MatZeroRowsColumns(Kt,n_dir,dir_index,1.0,NULL,NULL);
      MatAssemblyBegin(Kt,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(Kt,MAT_FINAL_ASSEMBLY);

      
      /******************************/
      /* Solve systems */ 

      KSPSolve(ksp,b,du);
      KSPGetConvergedReason(ksp,&reason);
      ierr = KSPGetIterationNumber(ksp,&ksp_its);CHKERRQ(ierr);

      if(rank == (nproc-1)){
        VecGetValues(du,1,&index_i,&d_lam);
      }
      lambda += d_lam; /* Bcast */
      VecAXPY( u, 1.0, du);

      /******************************/
      /* Evolute materials properties acording to the new u */

      error = bth_evolute( tot_step, CONV_OFF, &u);
      
      /******************************/
      /* Assembly RHS */
      error = bth_intforce( &f_int, &u); 

      error = VecWAXPY( R, -lambda, f_ext, f_int);

      VecCopy(R, b);
      VecScale( b, -1.0);
      VecSetValues(b,n_dir,dir_index,dir_zeros,INSERT_VALUES);

      VecDot(f_ext, u,   &p_fu);
      VecDot(f_ext, u_o, &p_fu_o);
      phi_u  = 0.5*( lambda   * p_fu - lambda_o * p_fu_o) - d_work_u_fixed;
      phi_d  = 0.5*( lambda_o * p_fu - lambda   * p_fu_o) - d_work_d_fixed;

      if(rank == (nproc-1)){

        if(tot_step == 1){
          phi = 0.0;
        }else if(tot_step == 2 ){
          phi = -(lambda - d_work);
        }else{
          if(internal_energy == 1){
            phi = phi_u/lambda;
          }else{
            phi = phi_d/lambda_o;
          }
        }
        VecSetValues(b,1,&index_i,&phi,INSERT_VALUES);

      }
      VecAssemblyBegin(b);
      VecAssemblyEnd(b);

      error = VecNorm( b, NORM_2, &norm);
      /******************************/

      if(error)
        return 1;
     
      its ++;

      PetscPrintf(PETSC_COMM_WORLD,"AL1 it : %2d |R|=%e KSP it: %4d reason: ",its,norm,ksp_its); 
      if(reason == KSP_CONVERGED_RTOL){
	  PetscPrintf(PETSC_COMM_WORLD,"KSP_CONVERGED_RTOL"); 
      }else if(reason == KSP_CONVERGED_ATOL){
	  PetscPrintf(PETSC_COMM_WORLD,"KSP_CONVERGED_ATOL"); 
      }else if(reason == KSP_CONVERGED_ITS){
	  PetscPrintf(PETSC_COMM_WORLD,"KSP_CONVERGED_ITS"); 
      }else{
	  PetscPrintf(PETSC_COMM_WORLD,"%3d",reason); 
      }
      PetscPrintf(PETSC_COMM_WORLD,"\n"); 

    }

    VecDot(f_ext, u,   &p_fu);
    VecDot(f_ext, u_o, &p_fu_o);
    d_work_u = 0.5*( lambda   * p_fu - lambda_o * p_fu_o);
    d_work_d = 0.5*( lambda_o * p_fu - lambda   * p_fu_o);

    if(tot_step == 2){
      d_work_u_fixed = d_work_u;
    }

    if(internal_energy == 1){

      if(d_work_d > a1 * d_work_u_fixed){
        internal_energy = 0;
        PetscPrintf(PETSC_COMM_WORLD,"dissipative power significant.\n"); 
        d_work_d_fixed = a2 * d_work_u_fixed;
      }

    }else{

      if(d_work_u < -d_work_u_fixed){

//        internal_energy = 1;
//        PetscPrintf(PETSC_COMM_WORLD,"reduction on internal energy.\n"); 
//        d_work_u_fixed *= -1.0;

      }

    }


    energy_u += d_work_u;
    energy_d += d_work_d;
    PetscFPrintf(PETSC_COMM_WORLD,fp_perf,"%e %d %d %e %e %e %e %e %e %e\n",
        calcu.t,tot_step,internal_energy,lambda,energy_u,energy_d,d_work_u,d_work_d,d_work_u_fixed,d_work_d_fixed); 

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
    }

    sub_step ++;
    tot_step ++;

  }

  fclose(fp_perf);

  return 0;

}

