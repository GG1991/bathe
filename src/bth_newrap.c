#include "bathe.h"

int bth_newrap_1(void){

  /* ARC-LENGHT METHOD - NORMAL PLANE */
  
  int           tot_step;
  int           sub_step;
  int           error;
  int           its;
  int           max_its;

  double        R_tol;
  double        norm;
  double        dt;
  double        t0;

  node_list_t * pn;

  Mat           Kt;

  Vec           u;       /* displacement   at time  t                           */
  Vec           f_ext;   /* external force at time  t                           */
  Vec           f_int;   /* internal force at time  t                           */
  Vec           R;       /* internal force at time  t                           */
  Vec           du;      /* delta u = u_t-1 - u_t (converged)                   */
  Vec           b;       /* rhs vector to solve for du = Kt^-1 b (include BCs)  */

  KSP           ksp;     /* linear solver context                               */

  PC            pc;      /* preconditioner context                              */

  /* Alloc memory */ 
  
  MatCreateAIJ(PETSC_COMM_WORLD,mesh.nnodes*DIM,mesh.nnodes*DIM,ntot*DIM,ntot*DIM,120,NULL,120,NULL,&Kt);
  VecCreateGhost(PETSC_COMM_WORLD,mesh.nnodes*DIM,ntot*DIM,mesh.nghost*DIM,(PetscInt*)ghost,&u); 
  VecDuplicate(u,&f_ext);
  VecDuplicate(u,&f_int);
  VecDuplicate(u,&R);
  VecDuplicate(u,&du);
  VecDuplicate(u,&b);
  
  /* Setting solver options */

  KSPCreate(PETSC_COMM_WORLD,&ksp);
  KSPSetType(ksp,KSPGMRES);
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCLU);
  KSPSetFromOptions(ksp);
  KSPSetOperators(ksp,Kt,Kt);
 
  /* init variables */

  R_tol   = calcu.params_d[0];
  max_its = calcu.params_i[0];

  /* set neumann in final values, lambda is going to scale it */

  calcu.t  = calcu.t0;
  t0       = calcu.t0;
  tot_step = 1;
  sub_step = 1;

  pn       = calcu.time.head;

  /* temporal loop  */
  
  while(pn){
    
    PetscPrintf(PETSC_COMM_WORLD,"t:%lf s n:%d\n",calcu.t,tot_step); 
    
    /* non-linear loop  */

    its  = 0;
    
    error = bth_set_dirichlet(&u);

    /* Evolute materials properties acording to the new u */

    error = bth_evolute( tot_step, 0, &u);

    error = bth_intforce(&f_int,&u); 

    error = bth_set_neumann(&f_ext, calcu.t);

    error = VecWAXPY( R, -1.0, f_ext, f_int); /* R = f_int - f_ext */

    /* Prepare RHS = b to solve */

    VecCopy(R, b);
    VecScale( b, -1.0);
    VecSetValues(b,n_dir,dir_index,dir_zeros,INSERT_VALUES);
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);

    error = VecNorm( b, NORM_2, &norm);
    PetscPrintf(PETSC_COMM_WORLD,"AL1 it : %2d |R|=%e \n",0,norm); 

    while( norm > R_tol && its < max_its){

      /******************************/
      /* Calculate tangent operator Kt */
      error = bth_calc_k(&Kt, &u);
      
      MatSetOption(Kt,MAT_KEEP_NONZERO_PATTERN,PETSC_TRUE);
      MatZeroRowsColumns(Kt,n_dir,dir_index,1.0,NULL,NULL);
      MatAssemblyBegin(Kt,MAT_FINAL_ASSEMBLY);
      MatAssemblyEnd(Kt,MAT_FINAL_ASSEMBLY);
      /******************************/

      KSPSolve(ksp,b,du);

      VecAXPY( u, 1.0, du);

      /* Evolute materials properties acording to the new u */

      error = bth_evolute( tot_step, 0, &u);

      error = bth_intforce( &f_int, &u); 

      /* Calculate Residue R */ 

      error = VecWAXPY( R, -1.0, f_ext, f_int);

      /* Prepare RHS = b to solve */

      VecCopy(R, b);
      VecScale( b, -1.0);
      VecSetValues(b,n_dir,dir_index,dir_zeros,INSERT_VALUES);
      VecAssemblyBegin(b);
      VecAssemblyEnd(b);

      error = VecNorm( b, NORM_2, &norm);

      if(error){
        return 1;
      }
     
      its ++;

      PetscPrintf(PETSC_COMM_WORLD,"AL1 it : %2d |R|=%e \n",its,norm); 

    }

    error = bth_evolute( tot_step, 1, &u);

    error = print_out(&u,&f_int,&f_ext,tot_step);

    dt = ( ((tcontrol_t*)pn->data)->tf - t0 ) / (((tcontrol_t*)pn->data)->st - 1);
    calcu.t += dt;

    if( sub_step >= ((tcontrol_t*)pn->data)->st ){
      sub_step  = 0;
      t0    = ((tcontrol_t*)pn->data)->tf;
      pn    = pn->next;
    }

    sub_step ++;
    tot_step ++;

  }

  return 0;

}

//VecView(u,PETSC_VIEWER_STDOUT_WORLD);
