/* structure types and constants
 * 
 */

#ifndef _FUNCT_H_
#define _FUNCT_H_

#include <slepceps.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "mesh.h"
#include "list.h"
#include "global.h"
#include "types.h"
#include "gmsh.h"
#include "fem.h"
#include "utils.h"
#include "fun.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_minmax.h>

//parser.c
int parse_input(void);
int parse_mesh(void);
int parse_mats(void);
int parse_mode(void);
int parse_boun(void);
int parse_func(void);
int parse_outp(void);
int parse_material(char *buff, pvl_t *mat);
int cmp_mat(void *a, void *b);
int parse_boundary(char *bufcpy, bou_t *bou);
int cmp_bou(void *a, void *b);
int cmp_time(void *a, void *b);
int get_int(char *buf, const char *name, int *a);
int get_char(char *buf, const char *name, char *a);

//lst2msh.c
int cpynode  (node_list_t *node_nl, node_t *node);
int cpyelemv (node_list_t *elem_nl, elem_t *elemv);
int cpyelems (node_list_t *elem_nl, elem_t *elems);

// bth_boundary.c
int bth_set_boundary(void);

// bth_force.c
int bth_extforce(Vec *f_ext);
int bth_intforce(Vec *f_int, Vec *u);
int bth_set_dirichlet(Vec * u);
int bth_set_neumann(Vec * f, double t);

/**************************************************/
// Iterative methods
//
//  bth_arclenght_x.c
int bth_arclength_1(void);
int bth_arclength_2(void);
int bth_arclength_3(void);
int bth_arclength_4(void); /* Arclength based on energy (internal & dissipated)*/

//  bth_newrap_x.c
int bth_newrap_1(void);    /* Newton-Raphson */
/**************************************************/

//bth_calc_k.c
int bth_calc_k(Mat * K, Vec *u);

/* bth_evolute.c */
int bth_evolute(int step,int flag_conv, Vec * u);

//bth_property.c
int bth_property( int e, int input_i[10], double input_d[10], double strain[VOI], double c_t[VOI][VOI],     
    double stress[VOI], int flag, int output_i[10], double output_d[10]);

int matlaw_003( int e, int input_i[10], double input_d[10], double strain[VOI], double c_t[VOI][VOI],     
    double stress[VOI], int flag, int output_i[10], double output_d[10]);

//bth_elemental.c
int bth_elem_intforce(int e,int *npe, int ix[NPE*DIM], double vec_e[NPE*DIM], Vec * xlocal);
int bth_elem_k(int e, int *npe, int ix[NPE*DIM], double K_e[NPE*DIM], Vec * u_local);

//bthstra.c
int bth_strain(int npe,double disp[NPE][DIM],double derivs[NPE][DIM],
    double B[VOI][NPE*DIM],double strain[VOI]);

int cmp_nod(void *a, void *b);

//procedur.c
int bth_init(int argc,char **argv);
int bth_fini(void);

//output.c
int print_struct(int step);
int print_vtk(int step, Vec * u);
int print_vtk_forces(int step, Vec * f_int, Vec * f_out, Vec * u);
int print_out(Vec * u, Vec * f_int, Vec * f_ext, int step);
int bthdisf(char *phys, double *norm, double val[2], Vec *u);
int printMatrixR(char *name,double *A, int m, int n);
int printMatrixRC(char *name,double **A, int m, int n);
int printVector(char *name,double *vec, int m);
int printMat(char *name,Mat *A, int n);
int printVec(char *name,Vec *vec, int m);
int vtkcode(int dim,int npe);

#endif
