#include "funct.h"

int bth_init(int argc,char **argv){

  /* Reads the imput file
   * Reads the mesh
   * Reads the mesh
   * Allocs mamory for K, x, b
   */

  int           error,i,d;

  node_list_t * pn;
  node_list_t * pp;
  node_list_t * pe;
  
  output_t    * po;

  SlepcInitialize(&argc,&argv,(char*)0,NULL);
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(PETSC_COMM_WORLD, &nproc);  
  calcu.exec = (nproc>1)?PARALLEL:SEQUENCIAL;
  if(argc == 1)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c:input file NS.\n\n"); 
    return 1;
  }

  //============================== 
  // PARCING INPUT FILE
  //============================== 

  list_init(&list_mater, sizeof(pvl_t), cmp_mat);
  list_init(&list_bound, sizeof(bou_t), cmp_bou);
  list_init(&list_fun1d, sizeof(f1d_t), cmp_f1d);
  list_init(&list_outpu, sizeof(output_t), NULL);

  strcpy(inputfile,argv[1]);
  PetscPrintf(PETSC_COMM_WORLD,"Parcing input file.\n");
  error=parse_input();
  if(error!=0){
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error parsing input file.\n");
    return 1;
  }

  //============================== 
  // READING MESH 
  //============================== 

  list_init(&list_nodes, sizeof(gmshN_t), gmsh_nodcmp);
  list_init(&list_ghost, sizeof(gmshN_t), gmsh_nodcmp);
  list_init(&list_elemv, sizeof(gmshE_t), cmp_int);
  list_init(&list_elems, sizeof(gmshE_t), cmp_int);
  list_init(&list_physe, sizeof(gmshP_t), cmp_int);    

  PetscPrintf(PETSC_COMM_WORLD,"Reading mesh.\n");
  error=gmsh_read(meshfile,epartfile,npartfile,rank,DIM,&list_nodes,&list_ghost,&list_elemv,&list_elems,&list_physe,&loc2gold,&loc2gnew,&npp,nproc);    

  if(error!=0){
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error reading mesh.\n"); 
    return 1;
  }

  ntot=0;
  for(i=0;i<nproc;i++)
    ntot+=npp[i];

  //==============================      
  // ERASE AND OPEN FILES 
  //==============================      

  pn=list_outpu.head;

  while(pn){

    po = (output_t*)pn->data;
    /* erase the output files if they exist */
    switch(po->kind){

      case 1:
        break;

      case 2:
        /* abrimos el archivo*/
        if(rank==0){
          po->kind_2.fp = fopen(po->kind_2.file,"w");
        }
        pp   = list_physe.head;
        while(pp){
          if( !strcmp(((gmshP_t*)pp->data)->name, po->kind_2.phys) )
            break;
          pp = pp->next;
        }

        if(!pp){
          PetscPrintf(PETSC_COMM_WORLD,"Physical entity %s NF in list_physe.\n",po->kind_2.phys);
          return 1;
        }

        if(((gmshP_t*)pp->data)->dim != 0){
          PetscPrintf(PETSC_COMM_WORLD,"Physical entity %s should refer to point in the list_physe.\n",po->kind_2.phys);
          return 1;
        }
        
        pe   = list_elems.head;

        while(pe){
          if( ((gmshE_t*)pe->data)->gmshid == ((gmshP_t*)pp->data)->gmshid ){
            po->kind_2.node = ((gmshE_t*)pe->data)->node[0] - 1;
            break;
          }
          pe = pe->next;
        }

        if(!pe){
          PetscPrintf(PETSC_COMM_WORLD,"Surface element with physe %s was not found in list_elems.\n",po->kind_2.phys);
          return 1;
        }
        break;

      default:
        break;

    }

    pn=pn->next;

  }


  //============================== 
  // PRINTING STRUCTURES
  //============================== 

  PetscPrintf(PETSC_COMM_WORLD,"Printing structures 1.\n");
  error = print_struct(1);   
  if(error!=0){
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error printing structures.\n"); 
    return 1;
  }


  //============================== 
  // CONSTRUCTING MESH
  //============================== 

  PetscPrintf(PETSC_COMM_WORLD,"Constructing mesh.\n");
  error=mesh_alloc(&list_nodes, &list_ghost, cpynode, &list_elemv, cpyelemv, &list_elems, cpyelems, &mesh);
  if(error){
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error allocating mesh.\n"); 
    return 1;
  }

  error=mesh_renum(&mesh,loc2gold,loc2gnew);
  if(error){
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error renumbering mesh nodes.\n"); 
    return 1;
  }

  PetscPrintf(PETSC_COMM_WORLD,"Allocating Matrices/Vectors.\n");
  ghost=(int*)calloc(mesh.nghost*DIM,sizeof(int));
  for(i=0;i<mesh.nghost;i++){
    for(d=0;d<DIM;d++){
      ghost[i*DIM+d]=loc2gnew[mesh.nnodes+i]*DIM+d;
    }
  }

  //==============================      
  // FINISH TO ASSEMBLY BOUNDARIES 
  //==============================      

  PetscPrintf(PETSC_COMM_WORLD,"Assemblying BCs.\n");
  error=bth_set_boundary();
  if(error)
  {
    PetscPrintf(PETSC_COMM_WORLD,"main.c: Error assembling BCs : ERROR: %d\n\n",error);
    return 1;
  }


  fem_inigau();
  if(error){
    PetscPrintf(PETSC_COMM_WORLD,"main.c:error gps init.\n\n"); 
    return 1;
  }

  //============================== 
  // PRINTING STRUCTURES
  //============================== 

  PetscPrintf(PETSC_COMM_WORLD,"Printing structures 2\n");
  error = print_struct(2);   
  if(error!=0){
    PetscPrintf(PETSC_COMM_WORLD,"main.c: Error printing structures: ERROR: %d\n",error); 
    return 1;
  }

  return 0;
}
