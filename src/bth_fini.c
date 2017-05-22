/* frees all the memory */

#include "funct.h"

int bth_fini(void){

  node_list_t * pn;
  output_t    * po;

  pn=list_outpu.head;

  while(pn){
    
    po = (output_t*)pn;

    if(po->kind==2) fclose(po->kind_2.fp);

    pn=pn->next;

  }

  SlepcFinalize();
  return 0;
}
