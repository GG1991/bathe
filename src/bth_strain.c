#include "bathe.h"

int bth_strain(int npe, double disp[NPE][DIM], double derivs[NPE][DIM], double B[VOI][NPE*DIM], 
    double strain[VOI]){

  int i,j,d;

  for(i=0;i<npe;i++){

    B[0][i*DIM+0] = derivs[i][0]; 
    B[0][i*DIM+1] = 0.0         ;
    B[0][i*DIM+2] = 0.0         ; 

    B[1][i*DIM+0] = 0.0         ;
    B[1][i*DIM+1] = derivs[i][1];
    B[1][i*DIM+2] = 0.0         ; 
    
    B[2][i*DIM+0] = 0.0         ;
    B[2][i*DIM+1] = 0.0         ;
    B[2][i*DIM+2] = derivs[i][2]; 

    B[3][i*DIM+0] = derivs[i][1];
    B[3][i*DIM+1] = derivs[i][0];
    B[3][i*DIM+2] = 0.0         ; 

    B[4][i*DIM+0] = 0.0         ;
    B[4][i*DIM+1] = derivs[i][2];
    B[4][i*DIM+2] = derivs[i][1];

    B[5][i*DIM+0] = derivs[i][2];
    B[5][i*DIM+1] = 0.0         ;
    B[5][i*DIM+2] = derivs[i][0]; 

  }

  for(j=0;j<VOI;j++){
    strain[j]=0.0; 
    for(i=0;i<npe;i++){
      for(d=0;d<DIM;d++){
        strain[j] += B[j][i*DIM+d]*disp[i][d];
      }
    }
  }

  return 0;

}
