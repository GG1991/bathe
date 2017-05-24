/* 
 * parser.c : parser rutine for BATHE inputs
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "global.h"
#include "bathe.h"
#include "list.h"
#include "types.h"
#include "mesh.h"

#define NBUF 128

int parse_input(void){

  if(access(inputfile,F_OK) == -1){
    PetscPrintf(PETSC_COMM_WORLD,"parser.c: file %s not found.\n",inputfile);
    return 1;
  }

  if(parse_mesh())return 1; 
  if(parse_mats())return 2;
  if(parse_mode())return 3;
  if(parse_func())return 4; 
  if(parse_boun())return 5; 
  if(parse_outp())return 5; 
  return 0;
}

/*******************************/

int parse_mesh(void){
  
  FILE * file;
  char * data,buf[NBUF];
  int    fl=0,com=0,ln=0;

  file= fopen(inputfile,"r");

  while(fgets(buf,NBUF,file))
  {
    ln++;
    data=strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Mesh")){
        fl=1;
      }else if(!strcmp(data,"$EndMesh")){
        if(com==7){ 
          return 0;
        }else if(nproc==1 && com>=1){
          return 0;
        }else if(nproc>1&&(!(com&4)||!(com&8)) ){
          PetscPrintf(PETSC_COMM_WORLD,"parser.c:part file NF.\n");
          return 1;
        }else{
          PetscPrintf(PETSC_COMM_WORLD,"parser.c:$Mesh sect BF.\n");
          return 1;    
        }
      }
      if(fl==2){
        if(strcmp(data,"mesh_file") == 0){    
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:meshfile line %d.\n",ln); 
            return 1;
          }
          strcpy(meshfile,data);   
          if(access(meshfile,F_OK) == -1){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:%s NF.\n",meshfile); 
            return 1;
          }
          com=com|1;  
        }else if(strcmp(data,"parfile_e") == 0){    
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:partfile line %d.\n",ln);
            return 1;
          }
          strcpy(epartfile,data);
          if(access(epartfile,F_OK) == -1){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:%s NF.\n",epartfile); 
            return 1;
          }
          com=com|2;
        }else if(strcmp(data,"parfile_n") == 0){    
          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:partfile line %d.\n",ln);
            return 1;
          }
          strcpy(npartfile,data);
          if(access(npartfile,F_OK) == -1){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:%s NF.\n",npartfile); 
            return 1;
          }
          com=com|4;     
        }else if(strcmp(data,"#") != 0){    
          PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d of %s\n",ln,inputfile); 
          return 1;
        }    
      }
      if(fl==1)
        fl=2;
    }
  }
  fclose(file);
  return 1;
}

/*******************************/

int parse_mats(void){

  FILE *file= fopen(inputfile,"r");   
  char *data,buf[NBUF],bufcpy[NBUF];
  int fl=0,ln=0;
  pvl_t mat;

  while(fgets(buf,NBUF,file))
  {
    ln++;
    strcpy(bufcpy,buf);
    data=strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Materials")){
        fl=1;
      }else if(!strcmp(data,"$EndMaterials")){
        if(list_mater.sizelist) 
          return 0;
        PetscPrintf(PETSC_COMM_WORLD,"parser.c:no materials specified.\n");
        return 1;     
      }
      if(fl==2)
      {
        if(data[0]!='#')
        {
          if(parse_material(bufcpy,&mat))
            return 1;
          if(list_insert_se(&list_mater,(void*)&mat)<0)
            return 1;
        }
      }
      if(fl==1)
        fl=2;
    }
  }
  fclose(file);
  return 1;
}

/*******************************/

int parse_mode(void){

  FILE   * file;   

  char   * data,buf[NBUF];

  int      fl,com,ln,auxi;

  double   auxd;
  double   dlam;

  list_t   list_st, list_tf, list_dt;        

  list_init(&calcu.time, sizeof(tcontrol_t), cmp_time);
  list_init(&list_tf   , sizeof(double)    , cmp_int );
  list_init(&list_dt   , sizeof(double)    , cmp_int );
  list_init(&list_st   , sizeof(int)       , cmp_int );
  
  fl  = 0;
  com = 0;
  ln  = 0;
  
  file= fopen(inputfile,"r");

  while(fgets(buf,NBUF,file)){

    ln++;
    data=strtok(buf," \n");

    if(data){

      if(!strcmp(data,"$Mode")){

        fl=1;

      }else if(!strcmp(data,"$EndMode")){


        node_list_t * pnt = list_tf.head;
        node_list_t * pnd = list_dt.head;
        node_list_t * pns = list_st.head;

        tcontrol_t time;
        
        if((com&8)==8){
          
          while(pnt && pnd){
            time.tf = *(double*)pnt->data;
            time.dt = *(double*)pnd->data;
            time.st = -1.0;
            if(list_insert_se(&calcu.time,(void*)&time)){
              return 1;
            }
            pnt=pnt->next;
            pnd=pnd->next;
          }

        }else if((com&16)==16){
          
          while(pnt && pns){
            time.tf = *(double*)pnt->data;
            time.dt = -1.0;
            time.st = *(int*)pns->data;
            if(list_insert_se(&calcu.time,(void*)&time)){
              return 1;
            }
            pnt=pnt->next;
            pns=pns->next;
          }

        }

        fclose(file);

        return 0;    

      }

      if(fl==2){

        if(strcmp(data,"geoappr") == 0){    

          data = strtok(NULL," \n");    

          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
          }
          
          if(!strcmp(data,"SD")){    
            calcu.geomodel=SD;
          }else{    
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:Invalid option at line %d.\n",ln);
            return 1;
          }

          com=com|1;

        }else if(strcmp(data,"timedep") == 0){

          data = strtok(NULL," \n");    

          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
          }

          if(!strcmp(data,"QSTATIC")){    
            calcu.timedep=QS;
          }else if(!strcmp(data,"DYNAMIC")){    
            calcu.timedep=QS;
          }else{    
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:Invalid option at line %d.\n",ln);
            return 1;
          }

          com=com|2;

        }else if(strcmp(data,"t0") == 0){

          data = strtok(NULL," \n");    

          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }

          calcu.t0=atof(data);

        }else if(strcmp(data,"tf") == 0){

          data = strtok(NULL," \n");

          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }

          while(data){
            auxd=atof(data);
            list_insertlast(&list_tf,(void*)&auxd);
            data=strtok(NULL," \n");
          }

          com=com|4;

        }else if(strcmp(data,"dt") == 0){

          data = strtok(NULL," \n");

          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }

          while(data){
            auxd=atof(data);
            list_insertlast(&list_dt,(void*)&auxd);
            data=strtok(NULL," \n");
          }

          com=com|8;

        }else if(strcmp(data,"steps") == 0){

          data = strtok(NULL," \n");

          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }

          while(data){
            auxi=atoi(data);
            list_insertlast(&list_st,(void*)&auxi);
            data=strtok(NULL," \n");
          }

          com=com|16;

        }else if(strcmp(data,"vtkstep") == 0){

          data = strtok(NULL," \n");

          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }

          calcu.vtkstep = atoi(data);

        }else if(strcmp(data,"steps") == 0){

          data = strtok(NULL," \n");

          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }

          while(data){
            auxi=atoi(data);
            list_insertlast(&list_st,(void*)&auxi);
            data=strtok(NULL," \n");
          }

          com=com|16;

        }else if(strcmp(data,"algorithm") == 0){

          data = strtok(NULL," \n");    

          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }

          if(!strcmp(data,"NEWTON_RAPHSON")){    

            calcu.algorithm = NR1;
            
            /* Reading R_tol parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"R_tol expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"R_tol")){
              PetscPrintf(PETSC_COMM_WORLD,"R_tol expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"R_tol value expected line %d.\n",ln);
              return 1;
            }
            calcu.params_d[0]=atof(data);

            /* Reading max_its parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"max_its")){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its value expected line %d.\n",ln);
              return 1;
            }
            calcu.params_i[0]=atoi(data);
            
          }else if(!strcmp(data,"ARC_LENGTH_1")){    

            /*
               algorithm    ARC_LENGTH_1
               dlen         0.0
               r_tol        1.0e-5
               max_its      100
               d_lam        0.0001
               k_restart    1
             */

            calcu.algorithm = AL1;
            list_init( &calcu.arclen_1.dlam, sizeof(double), NULL );

            /* Reading dlen parameter */ 
            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"dlen expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"dlen")){
              PetscPrintf(PETSC_COMM_WORLD,"dlen expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"dlen value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_1.dlen=atof(data);

            /* Reading r_tol parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"r_tol")){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_1.r_tol=atof(data);

            /* Reading max_its parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"max_its")){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_1.max_its=atoi(data);

            /* Reading d_lam parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"d_lam expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"d_lam")){
              PetscPrintf(PETSC_COMM_WORLD,"d_lam expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    

            while(data){

              if(!data){
                PetscPrintf(PETSC_COMM_WORLD,"d_lam value expected line %d.\n",ln);
                return 1;
              }
              dlam = atof(data);
              list_insertlast(&calcu.arclen_1.dlam,(void*)&dlam);

              data = strtok(NULL," \n");    

            }
            if(calcu.arclen_1.dlam.sizelist != list_tf.sizelist){
                PetscPrintf(PETSC_COMM_WORLD,"d_lam list should be the same size as time lapses line %d.\n",ln);
                return 1;
            }

            /* Reading k_restart parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"k_restart")){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_1.k_restart=atof(data);

          }else if(!strcmp(data,"ARC_LENGTH_2")){    

            /*
               algorithm    ARC_LENGTH_2
               dlen         0.0
               r_tol        1.0e-5
               max_its      100
               d_lam        0.0001
               k_restart    1
             */

            calcu.algorithm = AL2;
            list_init( &calcu.arclen_2.dlam, sizeof(double), NULL );

            /* Reading dlen parameter */ 
            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"dlen expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"dlen")){
              PetscPrintf(PETSC_COMM_WORLD,"dlen expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"dlen value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_2.dlen=atof(data);

            /* Reading r_tol parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"r_tol")){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_2.r_tol=atof(data);

            /* Reading max_its parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"max_its")){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_2.max_its=atoi(data);

            /* Reading d_lam parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"d_lam expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"d_lam")){
              PetscPrintf(PETSC_COMM_WORLD,"d_lam expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    

            while(data){

              if(!data){
                PetscPrintf(PETSC_COMM_WORLD,"d_lam value expected line %d.\n",ln);
                return 1;
              }
              dlam = atof(data);
              list_insertlast(&calcu.arclen_2.dlam,(void*)&dlam);

              data = strtok(NULL," \n");    

            }
            if(calcu.arclen_2.dlam.sizelist != list_tf.sizelist){
                PetscPrintf(PETSC_COMM_WORLD,"d_lam list should be the same size as time lapses line %d.\n",ln);
                return 1;
            }

            /* Reading k_restart parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"k_restart")){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_2.k_restart=atof(data);

          }else if(!strcmp(data,"ARC_LENGTH_3")){    

            /*
               algorithm    ARC_LENGTH_3
               d_work       1.0e-5
               r_tol        1.0e-5
               max_its      100
               d_lam        0.0001
               k_restart    1
             */

            calcu.algorithm = AL3;
            list_init( &calcu.arclen_3.dlam, sizeof(double), NULL );

            /* Reading dlen parameter */ 
            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"d_work expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"d_work")){
              PetscPrintf(PETSC_COMM_WORLD,"d_work expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"d_work value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_3.d_work=atof(data);

            /* Reading r_tol parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"r_tol")){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_3.r_tol=atof(data);

            /* Reading max_its parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"max_its")){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_3.max_its=atoi(data);

            /* Reading d_lam parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"d_lam expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"d_lam")){
              PetscPrintf(PETSC_COMM_WORLD,"d_lam expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    

            while(data){

              if(!data){
                PetscPrintf(PETSC_COMM_WORLD,"d_lam value expected line %d.\n",ln);
                return 1;
              }
              dlam = atof(data);
              list_insertlast(&calcu.arclen_3.dlam,(void*)&dlam);

              data = strtok(NULL," \n");    

            }
            if(calcu.arclen_3.dlam.sizelist != list_tf.sizelist){
                PetscPrintf(PETSC_COMM_WORLD,"d_lam list should be the same size as time lapses line %d.\n",ln);
                return 1;
            }

            /* Reading k_restart parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"k_restart")){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_3.k_restart=atof(data);

          }else if(!strcmp(data,"ARC_LENGTH_4")){    

            /*
               algorithm    ARC_LENGTH_4
               d_work       1.0e-5 0.001 0.6
               r_tol        1.0e-5
               max_its      100
               k_restart    1
             */

            calcu.algorithm = AL4;

            /* Reading d_work parameter */ 
            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"d_work expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"d_work")){
              PetscPrintf(PETSC_COMM_WORLD,"d_work expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"d_work value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_4.d_work=atof(data);

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"d_work value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_4.a1=atof(data);
            data = strtok(NULL," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"d_work value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_4.a2=atof(data);

            /* Reading r_tol parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"r_tol")){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"r_tol value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_4.r_tol=atof(data);

            /* Reading max_its parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"max_its")){
              PetscPrintf(PETSC_COMM_WORLD,"max_its expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"max_its value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_4.max_its=atoi(data);

            /* Reading k_restart parameter */ 

            if(!fgets(buf,NBUF,file)){
              return 1;
            }
            ln ++;

            data = strtok(buf," \n");    

            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart expected line %d.\n",ln);
              return 1;
            }
            if(strcmp(data,"k_restart")){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart expected line %d.\n",ln);
              return 1;
            }

            data = strtok(NULL," \n");    
            if(!data){
              PetscPrintf(PETSC_COMM_WORLD,"k_restart value expected line %d.\n",ln);
              return 1;
            }
            calcu.arclen_4.k_restart=atof(data);

          }else{    
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:Invalid option at line %d.\n",ln);
            return 1;
          }

        }else if(strcmp(data,"#") != 0){    
          PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
          return 1;
        }
      }

      if(fl==1)
        fl=2;

    }
  }

  return 1;    

}

/*******************************/

int parse_boun(void){

  FILE *file= fopen(inputfile,"r");   

  char *data,buf[NBUF],bufcpy[NBUF];

  int   fl=0,ln=0;

  bou_t bou;

  while(fgets(buf,NBUF,file)){
    ln++;
    strcpy(bufcpy,buf);
    data=strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Boundary")){
        fl=1;
      }else if(!strcmp(data,"$EndBoundary")){
        if(list_bound.sizelist){ 
          return 0;
        }
        PetscPrintf(PETSC_COMM_WORLD,"parser.c:no boundaries in %s.\n",inputfile);
        return 1;    
      }
      if(fl==2){
        if(data[0]!='#')
        {
          if(parse_boundary(bufcpy,&bou))
            return 1;
          if(list_insert_se(&list_bound,(void*)&bou)<0)
            return 1;
        }
      }
      if(fl==1)
        fl=2;
    }
  }
  fclose(file);
  return 1;

}

/*******************************/

int parse_func(void){

  FILE *file= fopen(inputfile,"r");   
  char *data,buf[NBUF];
  int fl=0,com=0,ln=0,i,fs;
  double *xy;
  f1d_t * f1d;
  list_t list_xy;

  list_init(&list_xy, 2*sizeof(double), cmp_dou);

  while(fgets(buf,NBUF,file)){
    ln++;
    data=strtok(buf," \n");
    if(data){
      if(!strcmp(data,"$Function")){
        fl=1;
        fs=0;
      }else if(!strcmp(data,"$EndFunction")){
        if(com==15){ 
          com=0;
          fl=0;
        }else if(fl==0){
          PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
          return 1;
        }else{
          PetscPrintf(PETSC_COMM_WORLD,"parser.c:incomplete $Function section.\n");
          return 1;
        }
      }
      if(fl==2){
        if(strcmp(data,"funcknd") == 0){
          data=strtok(NULL," \n");
          if(!data)
            return 1;
          if(!strcmp(data,"1D")){
            f1d=(f1d_t*)calloc(1,sizeof(f1d_t));
            com=com|1;
          }
        }else if(strcmp(data,"funcnum") == 0){
          if(!f1d)
            return 1;
          data=strtok(NULL," \n");
          if(!data)
            return 1;
          f1d->fnum=atoi(data);
          com=com|2;
        }else if(strcmp(data,"funcint") == 0){
          if(!f1d)
            return 1;
          data=strtok(NULL," \n");
          if(!data)
            return 1;
          if(!strcmp(data,"INTER1")){
            f1d->inter=atoi(data);
          }    
          com=com|4;
        }else if(strcmp(data,"start") == 0){
          fs=1;
          if( (com & 1 ) != 1){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:<funckind> should be < start>.\n");  
            return 1;
          }
        }else if(strcmp(data,"end") == 0){
          if(fs==0){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          if(!list_xy.sizelist){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          f1d->n=list_xy.sizelist;
          f1d->x=(double*)calloc(list_xy.sizelist,sizeof(double));
          f1d->y=(double*)calloc(list_xy.sizelist,sizeof(double));
          for(i=0;i<f1d->n;i++){
            f1d->x[i]=*(((double*)(list_xy.head->data))+0);
            f1d->y[i]=*(((double*)(list_xy.head->data))+1);
            list_delfirst(&list_xy);
          }
          if(list_insert_se(&list_fun1d,(void*)f1d)<0)
            return 1;
          com=com|8;
        }else if(fs){
          xy=(double*)calloc(2,sizeof(double));
          xy[0]=atof(data);
          data=strtok(NULL," \n");
          if(!data)
            return 1;
          xy[1]=atof(data);
          if(list_insert_se(&list_xy,(void*)xy)<0)
            return 1;
        }
      }
      if(fl==1)
        fl=2;
    }
  }
  fclose(file);
  return 0;    
}       

/*******************************/

int parse_material(char *buff, pvl_t *mat){

  char * data;
  int    i;

  if(!strlen(buff))
    return 1;

  data = strtok(buff," \n");
  if(!data)
    return 1;

  if(data[0]!='\"'||data[strlen(data)-1]!='\"')
    return 1;

  strcpy(mat->name,data);

  data = strtok(NULL," \n");
  if(!data)
    return 1;

  mat->mattyp=atoi(data);

  switch(mat->mattyp){

    case 001:
      mat->nparfix = 2; /* E nu */
      mat->nparvar = 0;
      break;

    case 002:
      mat->nparfix = 5; /* E nu ft Gf */
      mat->nparvar = 2*NGP; /* d(NGP) r(NGP) */
      break;

    case 003:
      mat->nparfix = 5;     /* E nu ft Gf len  */
      mat->nparvar = 2*NGP; /* d(NGP) r(NGP)   */
      break;

    default:
      return 1;

  }

  /* leemos los parametros fix*/
  for(i=0;i<mat->nparfix;i++){

    data = strtok(NULL," \n");
    if(!data) 
      return 1;

    mat->params[i]=atof(data);
  }

  return 0;
}

/******************************/

int parse_outp(void){

  FILE * file;

  char * data,buf[NBUF];

  int    fl, ln;

  output_t output;

  file = fopen(inputfile,"r");

  fl   = 0;
  ln   = 0;

  while(fgets(buf,NBUF,file)){

    ln++;
    data=strtok(buf," \n");

    if(data){

      if(!strcmp(data,"$Output")){

        if(fl) /* quiere decir que otro no se cerró*/
          return 1;

        fl=1;

      }else if(!strcmp(data,"$EndOutput")){

        list_insertlast(&list_outpu,(void*)&output);
        fl=0;

      }

      if(fl==2){

        if(strcmp(data,"kind") == 0){

          data = strtok(NULL," \n");    
          if(!data){
            PetscPrintf(PETSC_COMM_WORLD,"parser.c:BF line %d.\n",ln);
            return 1;
          }
          output.kind = atoi(data);
            
          switch(output.kind){
            
            case 1:
              break;

            case 2:

              /* kind = 2  u(node) vs f_int(node)  node -> phys */

              if(!fgets(buf,NBUF,file)) 
                return 1;
              ln ++;

              if( get_char(buf,"file",output.kind_2.file)) 
                return 1;

              if(!fgets(buf,NBUF,file)) 
                return 1;
              ln ++;

              if( get_char(buf,"phys",output.kind_2.phys)) 
                return 1;

              break;

            case 3:

              /* kind = 3  VTK of damage variables */

              if(!fgets(buf,NBUF,file)) 
                return 1;
              ln ++;

              if( get_int(buf,"step",&output.kind_3.step)) 
                return 1;

              break;

            case 4:

              /* kind = 4  VTK of f_int, f_ext and u */

              if(!fgets(buf,NBUF,file)) 
                return 1;
              ln ++;

              if( get_int(buf,"step",&output.kind_4.step)) 
                return 1;
              break;

            default:
              break;

          }
          
        }else{
          PetscPrintf(PETSC_COMM_WORLD,"parser.c:kind spected line %d.\n",ln);
          return 1;
        }


      }

      if(fl==1){
        fl=2;
      }

    }

  }

  fclose(file);

  return 0;

}

int get_int(char *buf, const char *name,int *a){

  char *data;

  data = strtok(buf," \n");    

  if(!data){
    PetscPrintf(PETSC_COMM_WORLD,"%s expected.\n",name);
    return 1;
  }

  if(strcmp(data,name)){
    PetscPrintf(PETSC_COMM_WORLD,"%s expected.\n",name);
    return 1;
  }

  data = strtok(NULL," \n");    
  if(!data){
    PetscPrintf(PETSC_COMM_WORLD,"%s value expected.\n",name);
    return 1;
  }
  *a = atoi(data);

  return 0;
}

int get_char(char *buf, const char *name,char *a){

  char *data;

  data = strtok(buf," \n");    

  if(!data){
    PetscPrintf(PETSC_COMM_WORLD,"%s expected.\n",name);
    return 1;
  }

  if(strcmp(data,name)){
    PetscPrintf(PETSC_COMM_WORLD,"%s expected.\n",name);
    return 1;
  }

  data = strtok(NULL," \n");    
  if(!data){
    PetscPrintf(PETSC_COMM_WORLD,"%s value expected.\n",name);
    return 1;
  }
  strcpy(a,data);

  return 0;
}

/******************************/

int cmp_mat(void *a, void *b){

  return (strcmp(((pvl_t*)a)->name,((pvl_t*)b)->name)==0)?0:-1;
}

/******************************/

int parse_boundary(char *buff, bou_t *bou){

  char *data;    
  node_list_t *pNode;
  if(!strlen(buff))
    return 1;
  data = strtok(buff," \n");
  if(!data)
    return 1;
  if(data[0]!='\"'||data[strlen(data)-1]!='\"')
    return 1;
  strcpy(bou->nam,data);
  data = strtok(NULL," \n");
  if(!data)
    return 1;
  bou->ord=atoi(data);
  data = strtok(NULL," \n");
  if(!data)
    return 1;
  if(strBin2Dec(data,&bou->kin))
    return 1;

  /*fx*/ 
  data = strtok(NULL," \n");
  if(!data)
    return 1;
  bou->nfx=atoi(data);
  pNode=list_fun1d.head;
  while(pNode){
    if(((f1d_t*)pNode->data)->fnum==bou->nfx)
      break;
    pNode=pNode->next;
  }
  if(!pNode){
    PetscPrintf(PETSC_COMM_WORLD,"parser.c:Function number %d NF.\n",bou->nfx);
    return 1;   
  }
  bou->fx=pNode->data;


  /*fy*/
  data = strtok(NULL," \n");
  if(!data)
    return 1;
  bou->nfy=atoi(data);
  pNode=list_fun1d.head;
  while(pNode){
    if(((f1d_t*)pNode->data)->fnum==bou->nfy)
      break;
    pNode=pNode->next;
  }
  if(!pNode){
    PetscPrintf(PETSC_COMM_WORLD,"parser.c:Function number %d NF.\n",bou->nfy);
    return 1;   
  }
  bou->fy=pNode->data;


  /*fz*/
  data = strtok(NULL," \n");
  if(!data)
    return 1;
  bou->nfz=atoi(data);
  pNode=list_fun1d.head;
  while(pNode){
    if(((f1d_t*)pNode->data)->fnum==bou->nfz)
      break;
    pNode=pNode->next;
  }
  if(!pNode){
    PetscPrintf(PETSC_COMM_WORLD,"parser.c:Function number %d NF.\n",bou->nfz);
    return 1;   
  }
  bou->fz=pNode->data;


  list_init(&bou->nl,sizeof(int),cmp_int);
  list_init(&bou->evl,sizeof(int),cmp_int);
  list_init(&bou->esl,sizeof(int),cmp_int);

  return 0;
}

/******************************/
int cmp_bou(void *a, void *b){

  if( ((bou_t*)a)->ord>((bou_t*)b)->ord){
    return 1;
  }else if( ((bou_t*)a)->ord==((bou_t*)b)->ord){
    return 0;
  }else{
    return -1;
  }
}

/******************************/

int cmp_time(void *a, void *b){

  if( ((tcontrol_t*)a)->tf>((tcontrol_t*)b)->tf){
    return 1;
  }else if( ((tcontrol_t*)a)->tf==((tcontrol_t*)b)->tf){
    return 0;
  }else{
    return -1;
  }
}
