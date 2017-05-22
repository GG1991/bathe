#include <petscksp.h>
#include "list.h"
#include "types.h"
#include "mesh.h"

#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#define DIM 3
#define NPE 8
#define NGP 8
#define VOI 6

#define NR_TOL 1.0e-6
#define NR_MAX_IT 100
#define NR_RESTART_K 100

enum { INITIALIZE, CALC_STRESS, CALC_TANGEN, CALC_EVO_CO, CALC_EVO_FI, GIVE_PARAMS};

#define PERTUR 1
#define SECANT 2

enum {QS, TR};            /*Quasi Static, transient*/
enum {SD, TL , UL};       /*Small deformations, Total Lagrangean, Updated Lagrangean*/
enum {K1,K2};             /*Calculation of elemental matrix by this K modes*/
enum {NR1,AL1,AL2,AL3,AL4};       /* NR=Newton-Raphson, AL!=Arc-Lenght sphere distance */
enum {SEQUENCIAL, PARALLEL};
enum {CONV_ON, CONV_OFF};

PetscViewer        viewer;
KSPConvergedReason reason;

int      rank;
int      nproc;   
int      nummat;                 
        
int    * loc2gold, *loc2gnew;
int      n_dir;
int      n_neu;
int    * npp;
int      ntot;
int    * dir_index;
int      its;
int    * ghost;

double * dir_value;
double * dir_zeros;
int    * neu_index;
double * neu_value;
                    
int      memory;

list_t list_nodes;
list_t list_ghost;
list_t list_elemv;
list_t list_elems;
list_t list_physe;
list_t list_mater;
list_t list_bound;
list_t list_fun1d;
list_t list_outpu;

mesh_t mesh;

int     nke,nre;

int     Istart,Iend;
double  kspnorm;

char inputfile[32];
char meshfile[32];
char epartfile[32];
char npartfile[32];
                    
calcu_t  calcu;

#endif
