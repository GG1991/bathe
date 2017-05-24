/* 

Bathe main routine

This code is used to solve the displacements field inside a solid body.

Main properties:

-> Non linear solver         : Newton Raphson (for boundary condition completly imposed)
-> Non linear control method : Arc-length energy based 
-> Material models           : Elastic Isotropic with Isotropic damage 

Approximations:

-> Small deformation approach

Input:

In the input we specified basically:

-> The mesh
-> The materials in the mesh
-> The kind of simulation (for example: boundary completly applied or automatic control load)
-> Boundary conditions
-> The outputs that user want to receive 

Output:

The output are clasified with "kinds" specified by users.

-> VTK plots of internal & external forces
-> VTK plots of displacement fields
-> VTK plots of damage variables
-> force vs. displacement on nodes in ASCII file

*/

#include "bathe.h"

int main(int argc,char **argv){

  //
  // Input parsing
  // Main variables initialization
  // Mesh reading and partition (future) 
  // Output files opening
  // 
  if(bth_init(argc,argv))
    goto ERROR;

  //
  // Selection of the calculation procedure
  //
  switch(calcu.algorithm){

    case NR1:
      if(bth_newrap_1())
        goto ERROR;
      break;

    case AL1:
      if(bth_arclength_1())
        goto ERROR;
      break;

    default:
      goto ERROR;

  }

  //
  // Code end
  // Deallocation of memory
  // Output files closing
  //
ERROR:
  if(bth_fini());

  return 0;
}
