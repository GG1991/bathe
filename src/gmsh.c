#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "list.h"
#include "gmsh.h"

#define MAXC 128

int gmsh_read(char *file,char *efile,char *nfile,int rank,int dim, list_t *list_nodes, list_t *list_ghost, list_t
*list_elemv, list_t *list_elems, list_t *list_phyce, int **loc2gold, int **loc2gnew, int **npp, int nproc){
       
    if(gmsh_readnodes(file,nproc,nfile,rank,list_nodes))                     return 1;   
    if(gmsh_readelemv(file,nproc,efile,rank,dim,list_elemv))                 return 2;
    if(gmsh_readelems(file,dim,list_elemv,list_elems))                       return 3;
    if(gmsh_readphys(file,list_elemv,list_elems,list_phyce))                 return 4;
    if(gmsh_detghosts(list_nodes,list_elemv,list_ghost))                     return 5;
    if(gmsh_readghosts(file,list_ghost))                                     return 6;
    if(gmsh_reenumerate(nfile,rank,list_nodes,list_ghost,loc2gold,loc2gnew,npp,nproc)) return 7;
    return 0;    
}

enum {SN,NN};

int gmsh_readnodes(char *file,int nproc,char *nfile,int rank,list_t *list_nodes){

    FILE *fm,*fn;
    int  sel=NN,p,skip=0,d;
    char buff1[MAXC],buff2[MAXC];
    char *data;
    gmshN_t node;

    fm = fopen(file,"r");
    if(!fm)
        return 1;
    if(nproc>1){
        fn = fopen(nfile,"r");
        if(!fn)
            return 1;
    }
    while(fgets(buff1,MAXC,fm)!=NULL){
        data=strtok(buff1," \n");
        if(!data)
            return 1;
        if(strcmp(data,"$Nodes")==0){
            if(fgets(buff1,MAXC,fm) == NULL)
                return 1;
            sel = SN;
            skip= 1;
        }else if(strcmp(data,"$EndNodes")==0){
            if(nproc>1)
                fclose(fn);
            fclose(fm);
            return 0;   
        }
        if(sel==SN && !skip){
            if(nproc>1){
                if(!fgets(buff2,MAXC,fn)) 
                    return 1;
                sscanf(buff2,"%d",&p);
                if(p<0) 
                    return 1;
            }else{
                p=rank;    
            }
            if(p==rank){
                node.n=atoi(data);
                for(d=0;d<3;d++){
                    data=strtok(NULL," \n");
                    node.coor[d]=atof(data); 
                }
                if(list_insertlast(list_nodes,(void*)&node)!=0)
                    return 1;
            }

        }
        if(skip)skip=0;

    }
    return 0;   
}

int gmsh_readelemv(char *file,int nproc,char *efile,int rank,int dim,list_t *list_elemv){

    FILE *fm,*fe;
    int  sel=NN,p,skip=0,i,ntags;
    char buff1[MAXC],buff2[MAXC];   
    char *data;
    gmshE_t  elemv;

    fm = fopen(file,"r");
    if(!fm)
        return 1;
    if(nproc>1){
        fe = fopen(efile,"r");
        if(!fe)
            return 1;
    }
    while(fgets(buff1,MAXC,fm)!=NULL){
        data=strtok(buff1," \n");
        if(strcmp(data,"$Elements")==0){
            if(fgets(buff1,MAXC,fm) == NULL)
                return 1;
            sel = SN;
            skip= 1;
        }else if(strcmp(data,"$EndElements")==0){
            if(nproc>1)fclose(fe);
            fclose(fm);return 0;    
        }
        if(sel==SN && !skip){
            data = strtok(NULL," \n");
            if(!data)
                return 1;
            if(gmsh_isvol(atoi(data),dim)<0)
                return 1;
            if(gmsh_isvol(atoi(data),dim)){
                if(nproc>1){
                    if(fgets(buff2,MAXC,fe)==NULL)
                        return 1;
                    sscanf(buff2,"%d",&p);
                    if(p<0) 
                        return 1;
                }else{
                    p=rank;
                    if(p==rank){

                    }
                }
                if(p==rank){
                    if(gmsh_npe(atoi(data))<0)
                        return 1;
                    elemv.npe = gmsh_npe(atoi(data));
                    data = strtok(NULL," \n");
                    if(!data)
                        return 1;
                    ntags = atoi(data);
                    data = strtok(NULL," \n");
                    if(!data)
                        return 1;
                    elemv.gmshid = atoi(data);
                    for(i=0;i<ntags-1;i++){
                        data = strtok(NULL," \n");
                        if(!data)
                            return 1;
                    }
                    for(i=0;i<elemv.npe;i++){
                        data = strtok(NULL," \n");
                        if(!data)
                            return 1;
                        elemv.node[i]=atoi(data);
                    }
                    if(list_insertlast(list_elemv,(void *)&elemv))
                        return 1;
                }
            }    
        }
        if(skip)
            skip=0;
    }

    return 0;   
}

int gmsh_readelems(char *file,int dim,list_t *list_elemv,list_t *list_elems){

    FILE *fm;
    int  sel=NN,skip=0,i,ntags;
    char buff1[MAXC];   
    char *data;
    gmshE_t  elems;

    fm = fopen(file,"r");
    if(!fm)
        return 1;
    while(fgets(buff1,MAXC,fm)!=NULL){
        data=strtok(buff1," \n");
        if(strcmp(data,"$Elements")==0){
            if(fgets(buff1,MAXC,fm) == NULL)
                return 1;
            sel = SN;
            skip= 1;
        }else if(strcmp(data,"$EndElements")==0){
            fclose(fm);return 0;    
        }
        if(sel==SN && !skip){
            data = strtok(NULL," \n");
            if(!data)
                return 1;
            if(gmsh_isvol(atoi(data),dim)<0)
                return 1;
            if(!gmsh_isvol(atoi(data),dim)){
                if(gmsh_npe(atoi(data))<0)
                    return 1;
                elems.npe = gmsh_npe(atoi(data));
                data = strtok(NULL," \n");
                if(!data)
                    return 1;
                ntags = atoi(data);
                data = strtok(NULL," \n");
                if(!data)
                    return 1;
                elems.gmshid = atoi(data);
                for(i=0;i<ntags-1;i++){
                    data = strtok(NULL," \n");
                    if(!data)
                        return 1;
                }
                for(i=0;i<elems.npe;i++){
                    data = strtok(NULL," \n");
                    if(!data)
                        return 1;
                    elems.node[i]=atoi(data);
                }
//                if(gmsh_elems_belongs(list_elemv,&elems,&elems.elemv)){
                if(list_insertlast(list_elems,(void*)&elems))
                  return 1;
//                }
            }
        }
        if(skip)
            skip=0;
    }    
    return 0;
}

int gmsh_readphys(char *file,list_t *list_elemv,list_t *list_elems,list_t *list_phys){

    FILE *fm;
    int  sel=NN,skip=0;
    char buff1[MAXC];   
    char *data;
    gmshP_t  phys;

    if(!file || !list_elemv || !list_elems || !list_phys)
        return 1;
    fm = fopen(file,"r");
    if(!fm)
        return 1;
    while(fgets(buff1,MAXC,fm)){
        data=strtok(buff1," \n");
        if(!data)
            return 1;
        if(strcmp(data,"$PhysicalNames")==0){
            if(fgets(buff1,MAXC,fm) == NULL)
                return 1;
            sel = SN;
            skip= 1;
        }else if(strcmp(data,"$EndPhysicalNames")==0){
            fclose(fm);
            return 0;    
        }
        if(sel==SN && !skip){
            phys.dim = atoi(data);
            data = strtok(NULL," \n");
            if(!data)
                return 1;
            phys.gmshid = atoi(data);
            data = strtok(NULL," \n");
            if(!data)
                return 1;
            strcpy(phys.name, data);
            //if(gmsh_phys_belongs(list_elemv,list_elems,&phys)){
                if(list_insertlast(list_phys,(void*)&phys))
                    return 1;
            //}

        }
        if(skip)skip=0;
    }
    return 0;
}

int gmsh_detghosts(list_t *list_nodes,list_t *list_elemv,list_t *list_ghost){

    /* Fills the list_ghost with ghosts nodes */
    int e,n,i,fl;
    if(!list_nodes || !list_elemv || !list_ghost)
        return 1;

    node_list_t *onode;
    node_list_t *anode;
    gmshE_t *elemv;
    gmshN_t *nodes;
    gmshN_t ghost;
    for(i=0;i<3;i++)
        ghost.coor[i]=-1;
    onode=list_elemv->head;
    for(e=0;e<list_elemv->sizelist;e++){
        elemv=(gmshE_t *)onode->data;
        for(n=0;n<elemv->npe;n++){
            fl=0;
            anode=list_nodes->head;
            for(i=0;i<list_nodes->sizelist;i++){
                nodes=(gmshN_t *)anode->data;
                if(nodes->n==elemv->node[n]){
                    fl=1;    
                    break;   
                }
                anode=anode->next;
            }
            if(!fl){
                ghost.n=elemv->node[n];
                list_insert_se(list_ghost,(void*)&ghost);
            }
        }
        onode=onode->next;
    }

    return 0;   
}

int gmsh_readghosts(char *file,list_t *list_ghost){
    
    FILE *fm;
    int  sel=NN,skip=0,d;
    char buff1[MAXC],*data;
    node_list_t *onode;
    gmshN_t *ghost;

    if(!file || !list_ghost)
        return 1;
    if(!list_ghost->sizelist)
        return 0;
    fm = fopen(file,"r");
    if(!fm)
        return 1;
    onode=list_ghost->head;
    ghost=(gmshN_t *)onode->data;
    while(fgets(buff1,MAXC,fm)){
        data=strtok(buff1," \n");
        if(!data)
            return 1;
        if(strcmp(data,"$Nodes")==0){
            if(!fgets(buff1,MAXC,fm))
                return 1;
            sel = SN;
            skip= 1;
        }else if(strcmp(data,"$EndNodes")==0){
            fclose(fm);return 0;    
        }
        if(sel==SN && !skip){
            if(atoi(data)==ghost->n){
                for(d=0;d<3;d++){
                    data=strtok(NULL," \n");
                    if(!data)
                        return 1;
                    ghost->coor[d]=atof(data);
                }
                onode=onode->next;
                if(onode)
                    ghost=(gmshN_t *)onode->data;
            }
        }
	if(skip)skip=0;
    }
    return 0;
}

int gmsh_reenumerate(char *nfile,int rank,list_t *list_nodes,list_t *list_ghost, int **loc2gold, int **loc2gnew, int **npp,int nproc){

    FILE *fn;
    int  p,i,n,na,nl,ng,nt,k,*nppa;
    char buff[MAXC];
    char *data;
    node_list_t *onode;
    gmshN_t *ghost;

    (*loc2gold) = (int*)calloc((list_nodes->sizelist+list_ghost->sizelist),sizeof(int));
    (*loc2gnew) = (int*)calloc((list_nodes->sizelist+list_ghost->sizelist),sizeof(int));
    (*npp) = (int*)calloc(nproc,sizeof(int));
    nppa= (int*)calloc(nproc,sizeof(int));
    if(nproc==1){
        if(list_ghost->sizelist)
            return 1;
        for(i=0;i<list_nodes->sizelist;i++){
            (*loc2gold)[i]=i;
            (*loc2gnew)[i]=i;
        }
        (*npp)[0]=list_nodes->sizelist;
        return 0;
    }

    fn = fopen(nfile,"r");
    if(!fn)
        return 1;
    while(fgets(buff,MAXC,fn)!=NULL){
        data=strtok(buff," \n");
        if(!data)
            return 1;
        p=atoi(data);
        (*npp)[p]++;
    }

    n=0;
    for(p=0;p<rank;p++)
        n+=(*npp)[p];

    rewind(fn);    

    ng=0;
    nl=0;
    nt=0;
    while(fgets(buff,MAXC,fn)!=NULL){
        data=strtok(buff," \n");
        if(!data)
            return 1;
        p=atoi(data);
        if(p==rank){
            (*loc2gold)[nl]=nt;
            (*loc2gnew)[nl]=n+nl;
            nl++;
        }else{
            onode=list_ghost->head;
            while(onode){
                ghost=(gmshN_t *)onode->data;
                if(ghost->n==(nt+1)){
                    break;
                }
                onode=onode->next;
            }
            if(onode){
                na=0;
                for(k=0;k<p;k++)
                    na+=(*npp)[k];
                (*loc2gold)[(*npp)[rank]+ng]=nt; 
                (*loc2gnew)[(*npp)[rank]+ng]=na+nppa[p]; 
                ng++;
            }
        }
        nppa[p]++;
        nt++;
    }
    return 0;
}

int gmsh_elems_belongs(list_t *list_elemv,gmshE_t *elems,int *nelemv){

    gmshE_t  *elemv;
    node_list_t  *node; 
    int nn,i,j,ne=0;    
    node = list_elemv->head;

    *nelemv=-1;
    while(node){
        elemv = (gmshE_t *)node->data;
        nn=0;
        for(i=0;i<elems->npe;i++){
            for(j=0;j<elemv->npe;j++){
                if(elemv->node[j]==elems->node[i]){ 
                    nn++; 
                    *nelemv=ne;
                    break;
                }
            }
        }
        if(nn==elems->npe)
            return 1;
        node=node->next;
        ne++;
    }    
    return 0;
}

int gmsh_phys_belongs(list_t *list_elemv,list_t *list_elems,gmshP_t *phys){

    gmshE_t  *elem;
    node_list_t  *node; 
    int i;

    for(i=0;i<2;i++){
        if(i==0){
            node = list_elemv->head;
        }else if(i==1){
            node = list_elems->head;
        }
        while(node){
            elem = (gmshE_t *)node->data;
            if(elem->gmshid==phys->gmshid) 
                return 1;
            node=node->next;
        }
    }
    return 0;    
}

int gmsh_isvol(int code,int dim){
    
    switch(dim){
	case 1:
	    return (code == 1)?1:0;
	case 2:
	    return (code == 2 || code == 3)?1:0;
	case 3:
	    return (code == 4 || code == 5 || code == 6)?1:0;
	default:
	    return -1;
    }        
}

int gmsh_npe(int code){
    
    switch(code){
	case 1:  
	    return 2; 
	case 2:    
	    return 3;
	case 3:    
	    return 4;
	case 4:    
	    return 4;
	case 5:    
	    return 8;
	case 6:   
	    return 6;
	case 15:   
	    return 1;
	default:    
	    return -1;
    }
}

int gmsh_nodcmp(void *a, void *b){
    gmshN_t *na,*nb;
    na=(gmshN_t *)a;
    nb=(gmshN_t *)b;
    if ( na->n > nb->n){
        return 1;
    }else if( na->n == nb->n){
        return 0;
    }else{
	return -1;
    }
}
