/*
 * output.c - This files contains the output functions
 * 
 * savevtk - Generates a vtk format file where are storage the mesh,
 * the material properties distribution and the flux
 * 
 * Autor : Guido Giuntoli
 * Last Modification: 29/12/2015
 *  
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "funct.h"
#include "global.h"
#include "list.h"
#include "gmsh.h"

#define   VTK_POINT         1
#define   VTK_LINE          3
#define   VTK_TRIANGLE      5
#define   VTK_QUADRANGLE    9
#define   VTK_TETRAHEDRON   10
#define   VTK_HEXAHEDRON    12
#define   VTK_6N_PRISM      13

int print_struct(int step){

  char          file_name[64];
  int           i,d,e,*pInt;
  gmshN_t     * p_gmsh_node;
  gmshE_t     * p_gmsh_elem;
  gmshP_t     * p_gmsh_phys;
  node_list_t * onode,*pNod;
  ps_t        * propES;
  pv_t        * pv;

  sprintf(file_name,"struct_r%d_s%d.bathe",rank,step);

  FILE *fout = fopen(file_name,"w");

  if(!fout)
  {
    printf("output.c: (output_print_structures) Error opening output file.");
    return 1;
  }
  fprintf(fout,"Structures - rank %d\n\n",rank);  

  fprintf(fout,"list_nodes : \n");
  fprintf(fout,"sizelist : %d\n",list_nodes.sizelist);
  onode = list_nodes.head;
  for(i=0;i<list_nodes.sizelist;i++){
    p_gmsh_node = (gmshN_t *)onode->data;
    fprintf(fout,"list_nodes - node : %d\n",i+1);
    fprintf(fout,"n : %d  coor :",p_gmsh_node->n);
    for(d=0;d<DIM;d++){
      fprintf(fout,"%4.3lf ",p_gmsh_node->coor[d]);
    }
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"list_ghost : \n");
  fprintf(fout,"sizelist : %d\n",list_ghost.sizelist);
  onode = list_ghost.head;
  for(i=0;i<list_ghost.sizelist;i++){
    p_gmsh_node = (gmshN_t *)onode->data;
    fprintf(fout,"list_ghost - node : %d\n",i+1);
    fprintf(fout,"n : %d  coor :",p_gmsh_node->n);
    for(d=0;d<DIM;d++)
      fprintf(fout,"%4.3lf ",p_gmsh_node->coor[d]);
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"list_elemv : \n");
  fprintf(fout,"sizelist : %d\n",list_elemv.sizelist);
  onode = list_elemv.head;
  for(i=0;i<list_elemv.sizelist;i++){
    p_gmsh_elem = (gmshE_t *)onode->data;
    fprintf(fout,"node : %d\n",i);
    fprintf(fout,"npe : %d ",p_gmsh_elem->npe);
    fprintf(fout,"gmshid : %d ",p_gmsh_elem->gmshid);
    fprintf(fout,"node : ");
    for(d=0;d<p_gmsh_elem->npe;d++)
      fprintf(fout,"%3d ",p_gmsh_elem->node[d]);
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"list_elems : \n");
  fprintf(fout,"sizelist : %d\n",list_elems.sizelist);
  onode = list_elems.head;
  for(i=0;i<list_elems.sizelist;i++){
    p_gmsh_elem = (gmshE_t *)onode->data;
    fprintf(fout,"node : %d\n",i);
    fprintf(fout,"npe : %d ",p_gmsh_elem->npe);
    fprintf(fout,"gmshid : %d ",p_gmsh_elem->gmshid);
    fprintf(fout,"elemv : %d ",p_gmsh_elem->elemv);
    fprintf(fout,"node : ");
    for(d=0;d<p_gmsh_elem->npe;d++)
      fprintf(fout,"%3d ",p_gmsh_elem->node[d]);
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"list_physe : \n");
  fprintf(fout,"sizelist : %d\n",list_physe.sizelist);
  onode = list_physe.head;
  for(i=0;i<list_physe.sizelist;i++){
    p_gmsh_phys = (gmshP_t *)onode->data;
    fprintf(fout,"node : %d\n",i);
    fprintf(fout,"name : %s ",p_gmsh_phys->name);
    fprintf(fout,"dim : %d ",p_gmsh_phys->dim);
    fprintf(fout,"gmshid : %d ",p_gmsh_phys->gmshid);
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"list_bound : \n\n");
  fprintf(fout,"sizelist : %d\n",list_bound.sizelist);
  onode = list_bound.head;
  for(i=0;i<list_bound.sizelist;i++){
    fprintf(fout,"node : %d ",i);
    fprintf(fout,"name : %s ",((bou_t *)onode->data)->nam); 
    fprintf(fout,"order: %d ",((bou_t *)onode->data)->ord); 
    fprintf(fout,"kind : %d ",((bou_t *)onode->data)->kin); 
    fprintf(fout,"dimS : %d ",((bou_t *)onode->data)->dim); 
    fprintf(fout,"nfx  : %d ",((bou_t *)onode->data)->nfx); 
    fprintf(fout,"nfy  : %d ",((bou_t *)onode->data)->nfy); 
    fprintf(fout,"nfz  : %d ",((bou_t *)onode->data)->nfz); 
    fprintf(fout,"nnodes: %d",((bou_t *)onode->data)->nl.sizelist);
    fprintf(fout,"\nnodes: ");
    pNod=((bou_t *)onode->data)->nl.head;
    while(pNod){
      fprintf(fout,"%4d ",*(int*)pNod->data);
      pNod=pNod->next;
    }
    fprintf(fout,"\nelems: ");
    pNod=((bou_t *)onode->data)->esl.head;
    while(pNod){
      fprintf(fout,"%4d ",*(int*)pNod->data);
      pNod=pNod->next;
    }
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"\nDirichlet indeces (%d) : ",n_dir);
  for(i=0;i<n_dir;i++){
    fprintf(fout,"%2d ",dir_index[i]);
  }
  fprintf(fout,"\nNeumann indeces (%d) : ",n_neu);
  for(i=0;i<n_neu;i++){
    fprintf(fout,"%2d ",neu_index[i]);
  }

  fprintf(fout,"\n list_fun1d : \n\n");
  fprintf(fout,"sizelist : %d\n",list_fun1d.sizelist);
  onode = list_fun1d.head;
  for(i=0;i<list_fun1d.sizelist;i++){
    fprintf(fout,"node : %d direction : %p\n",i,(f1d_t *)onode->data);
    fprintf(fout,"n : %d ",((f1d_t *)onode->data)->n);
    fprintf(fout,"inter : %d ",((f1d_t *)onode->data)->inter);
    fprintf(fout,"fnum : %d \n",((f1d_t *)onode->data)->fnum);
    for(d=0;d<((f1d_t *)onode->data)->n;d++){
      fprintf(fout," %lf  %lf\n",((f1d_t *)onode->data)->x[d],((f1d_t *)onode->data)->y[d]);
    }
    fprintf(fout,"\n");
    onode = onode ->next;
  }
  fprintf(fout,"\n");

  fprintf(fout,"npp : \n");
  for(i=0;i<nproc;i++)
    fprintf(fout,"rank %d : npp[%d] = %d \n",i,i,npp[i]);
  fprintf(fout,"\n");

  fprintf(fout,"mesh : \n\n");
  fprintf(fout,"coor : \n");
  fprintf(fout,"local : \n");
  for(i=0;i<mesh.nnodes;i++){
    fprintf(fout,"%3d: ",i); 
    for(d=0;d<3;d++)
      fprintf(fout,"  %+4.3lf ",mesh.node[i].coor[d]);
    onode=mesh.node[i].elemvL.head;        
    fprintf(fout,"\nelemvL: "); 
    e=0;        
    while(onode){
      pInt=(int*)onode->data;
      fprintf(fout,"%4d",*pInt); 
      onode=onode->next;
      e++; 
    }        
    while(e<8){
      fprintf(fout,"    "); 
      e++; 
    }
    onode=mesh.node[i].elemsL.head;        
    fprintf(fout,"\nelemsL: "); 
    e=0;        
    while(onode){
      pInt=(int*)onode->data;
      fprintf(fout,"%4d ",*pInt); 
      onode=onode->next;
      e++; 
    }        
    while(e<8){
      fprintf(fout,"    "); 
      e++; 
    }
    fprintf(fout,"\n");
  }
  fprintf(fout,"\n");

  fprintf(fout,"ghost : \n");
  for(i=0;i<mesh.nghost;i++){
    fprintf(fout,"%3d: ",i); 
    for(d=0;d<3;d++)
      fprintf(fout,"  %+4.3lf ",mesh.node[mesh.nnodes+i].coor[d]);
    onode=mesh.node[mesh.nnodes+i].elemvL.head;        
    fprintf(fout," elemvL: "); 
    e=0;        
    while(onode){
      pInt=(int*)onode->data;
      fprintf(fout,"%4d",*pInt); 
      onode=onode->next;
      e++; 
    }        
    while(e<8){
      fprintf(fout,"    "); 
      e++; 
    }
    onode=mesh.node[mesh.nnodes+i].elemsL.head;        
    fprintf(fout," elemsL: "); 
    e=0;        
    while(onode){
      pInt=(int*)onode->data;
      fprintf(fout,"%4d ",*pInt); 
      onode=onode->next;
      e++; 
    }        
    while(e<8){
      fprintf(fout,"    "); 
      e++; 
    }
    fprintf(fout,"\n");
  }
  fprintf(fout,"\n");

  fprintf(fout,"elemv : \n");
  for(i=0;i<mesh.nelemv;i++){
    fprintf(fout,"e=%4d npe: %4d\nnodel : ",i,mesh.elemv[i].npe);
    for(d=0;d<mesh.elemv[i].npe;d++)
      fprintf(fout,"%4d ",mesh.elemv[i].nodel[d]);
    fprintf(fout,"\nnodeg : ");
    for(d=0;d<mesh.elemv[i].npe;d++)
      fprintf(fout,"%4d ",mesh.elemv[i].nodeg[d]);
    pv = (pv_t*)mesh.elemv[i].prop;
    fprintf(fout,"\nmattyp : %d\n",pv->mattyp);
  }
  fprintf(fout,"\n");

  fprintf(fout,"elems : \n");
  for(i=0;i<mesh.nelems;i++){
    fprintf(fout,"e : %4d npe : %4d \nnodel : ",i,mesh.elems[i].npe);
    for(d=0;d<mesh.elems[i].npe;d++)
      fprintf(fout,"%4d ",mesh.elems[i].nodel[d]);
    fprintf(fout,"\nnodeg : ");
    for(d=0;d<mesh.elems[i].npe;d++)
      fprintf(fout,"%4d ",mesh.elems[i].nodeg[d]);
    propES=(ps_t*)mesh.elems[i].prop;
    fprintf(fout,"\ngmshid : %4d\n",propES->gmshid);
  }
  fprintf(fout,"\n");

  fprintf(fout,"loc2gold : \n");
  for(i=0;i<(mesh.nnodes+mesh.nghost);i++)
    fprintf(fout,"%5d %5d \n",i,loc2gold[i]);
  fprintf(fout,"\n");
  fprintf(fout,"loc2gnew : \n");
  for(i=0;i<(mesh.nnodes+mesh.nghost);i++)
    fprintf(fout,"%5d %5d \n",i,loc2gnew[i]);
  fprintf(fout,"\n");

  fclose(fout);     

  return 0;   

}

int print_vtk_forces(int step, Vec * f_int, Vec * f_ext, Vec * u){

  /* output file   */
  FILE    * vtkf;
  char      filevtk[32];

  int       n,e,d;
  int       count;
  int       index;

  double    vald;

  Vec       x_local;

  sprintf(filevtk,"force_r%d_s%d.vtk",rank,step);
  vtkf = fopen(filevtk, "w");

  /* MESH DATA */

  fprintf(vtkf, "# vtk DataFile Version 2.0\n");
  fprintf(vtkf, "Bathe\n");
  fprintf(vtkf, "ASCII\n");
  fprintf(vtkf, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(vtkf, "POINTS %d double\n", (mesh.nnodes+mesh.nghost));

  for (n=0;n<(mesh.nnodes+mesh.nghost);n++){
    for(d=0;d<3;d++)
      fprintf(vtkf, "%lf ", mesh.node[n].coor[d]);
    fprintf(vtkf, "\n");
  }

  count=0;
  for(e=0;e<mesh.nelemv;e++){
    count += mesh.elemv[e].npe + 1;
  }
  fprintf(vtkf, "CELLS %d %d\n", mesh.nelemv, count);
  for (e=0;e<mesh.nelemv;e++){
    fprintf(vtkf, "%d ", mesh.elemv[e].npe);
    for (n=0;n<mesh.elemv[e].npe;n++){
      fprintf(vtkf, "%d ", mesh.elemv[e].nodel[n]);
    }
    fprintf(vtkf, "\n");
  }

  fprintf(vtkf, "CELL_TYPES %i\n", mesh.nelemv);
  for (e=0;e<mesh.nelemv;e++){
    fprintf(vtkf, "%d\n",vtkcode(DIM,mesh.elemv[e].npe));  
  }


  /************************************************************/
  /* DATA IN NODES                                            */
  /************************************************************/

  fprintf(vtkf, "POINT_DATA %i\n",(mesh.nnodes+mesh.nghost));

  /* u */
  VecGhostUpdateBegin(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*u,&x_local);
  fprintf(vtkf,"VECTORS u FLOAT\n");
  for (n=0;n<(mesh.nnodes+mesh.nghost);n++){
    for (d=0;d<3;d++){
      index=n*DIM+d;
      VecGetValues(x_local,1,&index,&vald);
      fprintf(vtkf,"%lf ",vald);
    }
    fprintf(vtkf,"\n");
  }        

  /* f_int */
  VecGhostUpdateBegin(*f_int,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*f_int,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*f_int,&x_local);
  fprintf(vtkf,"VECTORS f_int FLOAT\n");
  for (n=0;n<(mesh.nnodes+mesh.nghost);n++){
    for (d=0;d<3;d++){
      index=n*DIM+d;
      VecGetValues(x_local,1,&index,&vald);
      fprintf(vtkf,"%lf ",vald);
    }
    fprintf(vtkf,"\n");
  }        

  /* f_int */
  VecGhostUpdateBegin(*f_ext,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*f_ext,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*f_ext,&x_local);
  fprintf(vtkf,"VECTORS f_ext FLOAT\n");
  for (n=0;n<(mesh.nnodes+mesh.nghost);n++){
    for (d=0;d<3;d++){
      index=n*DIM+d;
      VecGetValues(x_local,1,&index,&vald);
      fprintf(vtkf,"%lf ",vald);
    }
    fprintf(vtkf,"\n");
  }        

  fclose(vtkf);

  return 0;
}

int print_vtk(int step, Vec *u){

  int       n,e,d;
  int       count;
  int       index;

  double    vald;

  char      filevtk[32];
  char      ending[16];
  FILE    * vtkf;
  pv_t    * pv;
  Vec       xlocal;

  strncpy(filevtk,inputfile,strlen(inputfile)-6);
  filevtk[strlen(inputfile)-6]='\0';
  sprintf(ending,"_r%d_s%d.vtk",rank,step);
  strcat(filevtk,ending);
  vtkf = fopen(filevtk, "w");

  /************************************************************/
  /*Mesh geometry data                                        */
  /************************************************************/    
  fprintf(vtkf, "# vtk DataFile Version 2.0\n");
  fprintf(vtkf, "Bathe\n");
  fprintf(vtkf, "ASCII\n");
  fprintf(vtkf, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(vtkf, "POINTS %d double\n", (mesh.nnodes+mesh.nghost));

  for (n=0;n<(mesh.nnodes+mesh.nghost);n++){
    for(d=0;d<3;d++)
      fprintf(vtkf, "%lf ", mesh.node[n].coor[d]);
    fprintf(vtkf, "\n");
  }

  count=0;
  for(e=0;e<mesh.nelemv;e++)
    count += mesh.elemv[e].npe + 1;
  fprintf(vtkf, "CELLS %d %d\n", mesh.nelemv, count);
  for (e=0;e<mesh.nelemv;e++){
    fprintf(vtkf, "%d ", mesh.elemv[e].npe);
    for (n=0;n<mesh.elemv[e].npe;n++)
      fprintf(vtkf, "%d ", mesh.elemv[e].nodel[n]);
    fprintf(vtkf, "\n");
  }

  fprintf(vtkf, "CELL_TYPES %i\n", mesh.nelemv);
  for (e=0;e<mesh.nelemv;e++)
    fprintf(vtkf, "%d\n",vtkcode(DIM,mesh.elemv[e].npe));  


  /************************************************************/
  /* DATA IN NODES                                            */
  /************************************************************/
  fprintf(vtkf, "POINT_DATA %i\n",(mesh.nnodes+mesh.nghost));

  /* Ownership */
  fprintf(vtkf, "SCALARS ownership FLOAT\n");
  fprintf(vtkf, "LOOKUP_TABLE default\n");
  for(n=0;n<mesh.nnodes;n++)
    fprintf(vtkf,"%lf\n",1.0);
  for(n=0;n<mesh.nghost;n++)
    fprintf(vtkf,"%lf\n",0.0);

  /* Displacement */
  VecGhostUpdateBegin(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*u,&xlocal);
  fprintf(vtkf,"VECTORS displac FLOAT\n");
  for (n=0;n<(mesh.nnodes+mesh.nghost);n++){
    for (d=0;d<3;d++){
      index=n*DIM+d;
      VecGetValues(xlocal,1,&index,&vald);
      fprintf(vtkf,"%lf ",vald);
    }
    fprintf(vtkf,"\n");
  }        

  /************************************************************/
  /* DATA IN CELLS                                            */
  /************************************************************/      
  fprintf(vtkf, "CELL_DATA %i\n",mesh.nelemv);

  fprintf(vtkf, "SCALARS damage FLOAT\n");
  fprintf(vtkf, "LOOKUP_TABLE default\n");
  for(e=0;e<mesh.nelemv;e++){
    pv=(pv_t*)mesh.elemv[e].prop;
    fprintf(vtkf, "%lf\n",pv->params[0]);
  }

  fprintf(vtkf, "TENSORS strain FLOAT\n");
  for(e=0;e<mesh.nelemv;e++){
    pv=(pv_t*)mesh.elemv[e].prop;
    fprintf(vtkf, "%lf %lf %lf\n", pv->strain[0],pv->strain[4],pv->strain[5]);
    fprintf(vtkf, "%lf %lf %lf\n", pv->strain[4],pv->strain[1],pv->strain[3]);
    fprintf(vtkf, "%lf %lf %lf\n", pv->strain[5],pv->strain[3],pv->strain[2]);
    fprintf(vtkf, "\n"); 
  }

  fprintf(vtkf, "TENSORS stress FLOAT\n");
  for(e=0;e<mesh.nelemv;e++){
    pv=(pv_t*)mesh.elemv[e].prop;
    fprintf(vtkf, "%lf %lf %lf\n", pv->stress[0],pv->stress[4],pv->stress[5]);
    fprintf(vtkf, "%lf %lf %lf\n", pv->stress[4],pv->stress[1],pv->stress[3]);
    fprintf(vtkf, "%lf %lf %lf\n", pv->stress[5],pv->stress[3],pv->stress[2]);
    fprintf(vtkf, "\n"); 
  }

  fclose(vtkf);

  return 0;

}


int print_vtk_damage(int step, Vec *u){

  FILE        * vtkf;
  char          filevtk[] = "damage";
  char          ending[16];
  int           i,n,e,d,gp;
  int           ii;
  int           count;
  int           ngp;
  int           npe;
  int           input_i[10];
  int           error;
  double        disp[NPE][DIM];
  double        output_d[10];
  double        tau_ave, r_ave, d_ave;
  double        vol;
  double        strain[VOI];
  double        coor[NPE][DIM];
  double        B[VOI][NPE*DIM];
  double        der[NPE][DIM];
  double        jac[DIM][DIM];
  double        ijac[DIM][DIM];
  double        det;
  double    *** ode;
  double      * wp;
  pv_t        * pv;
  Vec           u_local;


  sprintf(ending,"_r%d_s%d.vtk",rank,step);
  strcat(filevtk,ending);

  vtkf = fopen(filevtk, "w");

  /************************************************************/
  /*Mesh geometry data                                        */
  /************************************************************/    

  fprintf(vtkf, "# vtk DataFile Version 2.0\n");
  fprintf(vtkf, "Bathe\n");
  fprintf(vtkf, "ASCII\n");
  fprintf(vtkf, "DATASET UNSTRUCTURED_GRID\n");
  fprintf(vtkf, "POINTS %d double\n", (mesh.nnodes+mesh.nghost));

  for (n=0;n<(mesh.nnodes+mesh.nghost);n++){
    for(d=0;d<3;d++)
      fprintf(vtkf, "%lf ", mesh.node[n].coor[d]);
    fprintf(vtkf, "\n");
  }

  count=0;
  for(e=0;e<mesh.nelemv;e++)
    count += mesh.elemv[e].npe + 1;
  fprintf(vtkf, "CELLS %d %d\n", mesh.nelemv, count);
  for (e=0;e<mesh.nelemv;e++){
    fprintf(vtkf, "%d ", mesh.elemv[e].npe);
    for (n=0;n<mesh.elemv[e].npe;n++)
      fprintf(vtkf, "%d ", mesh.elemv[e].nodel[n]);
    fprintf(vtkf, "\n");
  }

  fprintf(vtkf, "CELL_TYPES %i\n", mesh.nelemv);
  for (e=0;e<mesh.nelemv;e++)
    fprintf(vtkf, "%d\n",vtkcode(DIM,mesh.elemv[e].npe));  

  /************************************************************/
  /* DATA IN CELLS                                            */
  /************************************************************/      

  VecGhostUpdateBegin(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*u,&u_local);

  fprintf(vtkf, "CELL_DATA %i\n",mesh.nelemv);

  fprintf(vtkf, "SCALARS tau_ave FLOAT\n");
  fprintf(vtkf, "LOOKUP_TABLE default\n");

  for(e=0;e<mesh.nelemv;e++){

    npe = mesh.elemv[e].npe;
    ngp = mesh.elemv[e].ngp;
    pv=(pv_t*)mesh.elemv[e].prop;

    for(i=0;i<npe;i++){
      for(d=0;d<DIM;d++){
        ii = mesh.elemv[e].nodel[i]*DIM+d;
        VecGetValues(u_local,1,&ii,&disp[i][d]);
      }
    }
    tau_ave = vol = 0.0;

    fem_calwei(npe,DIM,&wp);
    fem_calode(npe,DIM,&ode);

    for(gp=0;gp<ngp;gp++){

      fem_caljac3(coor,ode,npe,gp,jac);

      fem_invjac3(jac,ijac,&det);
      if(det<0.0) return 1;

      fem_calder3(ijac,npe,gp,ode,der);

      /* el strain lo calculamos para juntar tau que no se guarda en memoria */
      bth_strain(npe,disp,der,B,strain); 

      error = bth_property( e, input_i, NULL, NULL, NULL, NULL, GIVE_PARAMS, NULL, output_d);
      if(error) return 1;

      tau_ave += output_d[0]*wp[gp]*det;
      vol     += wp[gp]*det;
     
    }
    tau_ave /= vol;

    if(pv->mattyp == 3){
      fprintf(vtkf, "%lf\n",tau_ave);
    }else{
      fprintf(vtkf, "%lf\n",0.0);
    }

  }

  /**************************************************/
  fprintf(vtkf, "SCALARS r_ave FLOAT\n");
  fprintf(vtkf, "LOOKUP_TABLE default\n");

  for(e=0;e<mesh.nelemv;e++){

    npe = mesh.elemv[e].npe;
    pv=(pv_t*)mesh.elemv[e].prop;

    for(i=0;i<npe;i++){
      for(d=0;d<DIM;d++){
        ii = mesh.elemv[e].nodel[i]*DIM+d;
        VecGetValues(u_local,1,&ii,&disp[i][d]);
      }
    }

    fem_calwei(npe,DIM,&wp);
    fem_calode(npe,DIM,&ode);

    for(gp=0;gp<ngp;gp++){

      fem_caljac3(coor,ode,npe,gp,jac);

      fem_invjac3(jac,ijac,&det);
      if(det<0.0) return 1;

      fem_calder3(ijac,npe,gp,ode,der);

      error = bth_property( e, input_i, NULL, NULL, NULL, NULL, GIVE_PARAMS, NULL, output_d);
      if(error) return 1;

      r_ave += output_d[1]*wp[gp]*det;
      vol     += wp[gp]*det;
     
    }
    r_ave /= vol;

    if(pv->mattyp == 3){
      fprintf(vtkf, "%lf\n",r_ave);
    }else{
      fprintf(vtkf, "%lf\n",0.0);
    }

  }
    
  /**************************************************/
  fprintf(vtkf, "SCALARS d_ave FLOAT\n");
  fprintf(vtkf, "LOOKUP_TABLE default\n");

  for(e=0;e<mesh.nelemv;e++){

    npe = mesh.elemv[e].npe;

    for(i=0;i<npe;i++){
      for(d=0;d<DIM;d++){
        ii = mesh.elemv[e].nodel[i]*DIM+d;
        VecGetValues(u_local,1,&ii,&disp[i][d]);
      }
    }

    pv=(pv_t*)mesh.elemv[e].prop;

    fem_calwei(npe,DIM,&wp);
    fem_calode(npe,DIM,&ode);

    for(gp=0;gp<ngp;gp++){

      fem_caljac3(coor,ode,npe,gp,jac);

      fem_invjac3(jac,ijac,&det);
      if(det<0.0) return 1;

      fem_calder3(ijac,npe,gp,ode,der);

      error = bth_property( e, input_i, NULL, NULL, NULL, NULL, GIVE_PARAMS, NULL, output_d);
      if(error) return 1;

      d_ave += output_d[2]*wp[gp]*det;
      vol     += wp[gp]*det;

    }
    r_ave /= vol;

    if(pv->mattyp == 3){
      fprintf(vtkf, "%lf\n",d_ave);
    }else{
      fprintf(vtkf, "%lf\n",0.0);
    }
  }

  return 0;
}

int print_out(Vec *u, Vec *f_int, Vec *f_ext, int step){

  int           i;
  int           indeces[10];
  int           step_o;
  int           error;
  double        values[10];
  Vec           u_local;
  Vec           f_local;
  node_list_t * pNod;
  output_t    * po;
        
  VecGhostUpdateBegin(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*u,&u_local);

  VecGhostUpdateBegin(*f_int,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*f_int,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*f_int,&f_local);

  error = 0;

  pNod = list_outpu.head;

  while(pNod){

    po=(output_t *)pNod->data;

    switch(po->kind){

      case 1:

//        double        vloc[2], vglo[2];
//
//        bthdisf(po->phys, po->norm, vloc, u);
//        MPI_Reduce(vloc, vglo, 2, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
//
//        if(!rank){
//          fl=fopen(po->file,"a");
//        }
//
//        fprintf(fl,"%lf %lf \n",vglo[0],vglo[1]);
//
//        fclose(fl);
        break;

      case 2:

        /* 

           plots f_int(n) vs u(n) being n a node specified by phys 
           the index is saved on po->params_i[0]  

           $Output
             file f_vs_u.dat
             kind 2
             phys "MOV_B"
           $EndOutput

        */

        for(i=0;i<3;i++){
          indeces[0+i] = po->kind_2.node * DIM + i;
        }

        VecGetValues(u_local,3,indeces,&values[0]);
        VecGetValues(f_local,3,indeces,&values[3]);
        fprintf(po->kind_2.fp,"%e %e %e %e %e %e\n",values[0],values[1],values[2],values[3],values[4],values[5]);
        break;

      case 3:

        /* 
           plots averange damage on elements on a vtk file, the total damage is printed on a file 
           an is the total on those elements with physical entity physe 

           $Output
             kind 3
             step 10
             phys "MAT"
           $EndOutput
           
           output : damage_s#_r#.vtk    d, tau on elements
         */

        step_o = po->kind_3.step;

        if(step % step_o == 0){
          error = print_vtk_damage(step, u);
        }
        break;

      case 4:

        /* 
           plots f_int, f_ext and u vectors in a vtk format

           $Output
             kind 4
             step 10
           $EndOutput
           
           output : ffu_s#_r#.vtk    d, tau on elements
         */

        step_o = po->kind_4.step;

        if(step % step_o == 0){
          error = print_vtk_forces( step, f_int, f_ext, u);
        }
        break;

      default:
        return 1;

    }

    if(error){
      return 1;
    }

    pNod=pNod->next;
  }

  return 0;
}

int bthdisf(char *phys, double *norm, double val[2], Vec *u){

  double        strmat[3][3],disp[3],dispa,forc[3],aux,area,areat;
  int           d1,d2,locind;
  node_list_t * pNodB,*pNodE;
  Vec           xlocal;

  VecGhostUpdateBegin(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostUpdateEnd(*u,INSERT_VALUES,SCATTER_FORWARD);
  VecGhostGetLocalForm(*u,&xlocal);

  memset(val,0.0,2*sizeof(double));
  pNodB=list_bound.head;
  while(pNodB){
    if(!strcmp(((bou_t*)pNodB->data)->nam,phys)) {
      areat=0.0;
      pNodE=((bou_t*)pNodB->data)->esl.head;
      while(pNodE){

        mesh_calcarea(&mesh,&mesh.elems[*(int*)pNodE->data],DIM,&area);
        areat+=area;

        memset(disp,0.0,3*sizeof(double));
        for(d1=0;d1<mesh.elems[*(int*)pNodE->data].npe;d1++){
          for(d2=0;d2<DIM;d2++){
            locind=mesh.elems[*(int*)pNodE->data].nodel[d1]*DIM+d2;
            VecGetValues(xlocal,1,&locind,&dispa);
            disp[d2]+=dispa;
          }                                      
        }
        aux=0.0;
        for(d2=0;d2<DIM;d2++)
          aux+=disp[d2]*norm[d2];
        aux*=area;
        val[0]+=aux/mesh.elems[*(int*)pNodE->data].npe;

        aux=0.0;
        voigt2mat(((pv_t*)mesh.elemv[((ps_t*)mesh.elems[*(int*)pNodE->data].prop)->elemv].prop)->stress,strmat);
        for(d1=0;d1<3;d1++){
          forc[d1]=0.0;
          for(d2=0;d2<3;d2++)
            forc[d1]+=strmat[d1][d2]*norm[d2];
          aux+=forc[d1]*norm[d1];
        }
        val[1]+=aux*area;
        pNodE=pNodE->next;
      }
      val[0]/=areat;
    }
    pNodB=pNodB->next;
  }
  return 0;
}

int printMatrixR(char *name,double *A, int m, int n){

  int    i,j;
  FILE * f = fopen(name,"w");

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      fprintf(f,"%lf ",A[i*n+j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);   
  return 0;
}

int printMatrixRC(char *name,double **A, int m, int n){

  int    i,j;

  FILE * f = fopen(name,"w");

  for(i=0;i<m;i++){
    for(j=0;j<n;j++){
      fprintf(f,"%lf ",A[i][j]);
    }
    fprintf(f,"\n");
  }
  fclose(f);   
  return 0;
}

int printMat(char *name,Mat *A,int n){

  int      i,j;
  double   val;
  FILE   * f = fopen(name,"w");

  for(i=Istart;i<Iend;i++){
    for(j=0;j<n;j++){
      MatGetValues(*A,1,&i,1,&j,&val);
      fprintf(f,"%lf ",val);
    }
    fprintf(f,"\n");
  }
  fclose(f);   
  return 0;
}

int printVector(char *name,double *vec, int m){

  int j;

  FILE *f=fopen(name,"w");

  for(j=0;j<m;j++){
    fprintf(f,"%lf %lf\n ",vec[j],vec[j]);
  }

  fclose(f);   
  return 0;
}

int printVec(char *name,Vec *vec, int m){

  /* Prints the values of a Vec, it can include ghost values*/

  int      j;
  double   val;
  FILE   * f = fopen(name,"w");
  Vec      xlocal;

  VecGhostGetLocalForm(*vec,&xlocal);

  for(j=0;j<m;j++){
    VecGetValues(xlocal,1,&j,&val);
    fprintf(f,"%lf %lf\n ",val,val);
  }

  fclose(f);   
  return 0;
}

int vtkcode(int dim,int npe){

  switch(dim){
    case 1:
      switch(npe){
        case 2 :
          return VTK_LINE;
        default:
          return -1;
      }
    case 2:
      switch(npe){
        case 3 :
          return VTK_TRIANGLE;
        case 4 :
          return VTK_QUADRANGLE;
        default:
          return -1;
      }
    case 3:
      switch(npe){
        case 4 :
          return VTK_TETRAHEDRON;
        case 6 :
          return VTK_6N_PRISM;  
        case 8 :
          return VTK_HEXAHEDRON;  
        default:
          return -1;
      }
    default:
      return -1;
  }
}
