/* Bathe main program*/

#include "funct.h"

int main(int argc,char **argv){

  if(bth_init(argc,argv))
    goto ERROR;

  switch(calcu.algorithm){

    case NR1:
      if(bth_newrap_1())
        goto ERROR;
      break;

    case AL1:
      if(bth_arclength_1())
        goto ERROR;
      break;

    case AL2:
      if(bth_arclength_2())
        goto ERROR;
      break;

    case AL3:
      if(bth_arclength_3())
        goto ERROR;
      break;
    
    case AL4:
      if(bth_arclength_4())
        goto ERROR;
      break;

    default:
      goto ERROR;

  }

ERROR:
  if(bth_fini());

  return 0;
}
