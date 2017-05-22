
#ifndef _TYPES_H_
#define _TYPES_H_

#include "list.h"
#include "mesh.h"
#include "fun.h"


typedef struct bit{
    
    unsigned int val : 1;
    
}bit_t;

typedef struct _ps_t{
    
    int gmshid;  /* ID that correspond to evewry element in the gmshfile */
    int elemv;   /* The elemv local numeration at which belongs */
    
}ps_t;

typedef struct _bou_t{

    char   nam[16];
    int    kin;
    int    ord;
    int    dim;    
    int    nfx;
    int    nfy;
    int    nfz;
    f1d_t  *fx; 
    f1d_t  *fy;
    f1d_t  *fz;
    list_t nl;    
    list_t evl;    
    list_t esl;    

}bou_t;


typedef struct _pv_t{
    
    int    mattyp;
    double *params;
    double strain[6];
    double stress[6];
    
}pv_t;

typedef struct _pvl_t{
    
    char     name[16];
    int      mattyp;
    int      nparam;
    int      nparfix;
    int      nparvar;
    double   params[10];
    
}pvl_t;

/*******************************/

typedef struct _kind_1_t
{
  char     phys[16];
}
kind_1_t;

typedef struct _kind_2_t
{
  FILE     * fp;
  char       phys[16];
  char       file[16];
  int        node;
}
kind_2_t;

typedef struct _kind_3_t
{
  int      node;
  int      step;
}
kind_3_t;

typedef struct _kind_4_t
{
  int      node;
  int      step;
}
kind_4_t;

typedef struct _output_t{

  int        kind;

  kind_1_t   kind_1;
  kind_2_t   kind_2;
  kind_3_t   kind_3;
  kind_4_t   kind_4;

}output_t;
                                            
/*******************************/

typedef struct _tcontrol_t{
    
    double tf;
    double dt;
    double st;

}tcontrol_t;

typedef struct _arclength_1_t{

    int     max_its;
    int     k_restart;
    double  r_tol;
    double  dlen;
    list_t  dlam;       /* dlam increments used for arc-length method */ 

}arclength_1_t;

typedef struct _arclength_2_t{

    int     max_its;
    int     k_restart;
    double  r_tol;
    double  dlen;
    list_t  dlam;       /* dlam increments used for arc-length method */ 

}arclength_2_t;

typedef struct _arclength_3_t{

    int     max_its;
    int     k_restart;
    double  r_tol;
    double  d_work;
    list_t  dlam;       /* dlam increments used for arc-length method */ 

}arclength_3_t;

typedef struct _arclength_4_t{

    int     max_its;
    int     k_restart;
    double  r_tol;
    double  d_work;
    double  a1;
    double  a2;

}arclength_4_t;

typedef struct{
    
    int             geomodel;
    int             timedep;
    int             vtkstep;
    int             algorithm;
    int             kmode;  
    int             mode;
    int             exec;
    int             params_i[5];
    double          params_d[5];
    
    double          t0;
    double          t;
    list_t          time;

    arclength_1_t   arclen_1;
    arclength_2_t   arclen_2;
    arclength_3_t   arclen_3;
    arclength_4_t   arclen_4;

}calcu_t;

/*******************************/


#endif
