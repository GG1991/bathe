#include "utils.h"

int cmp_int(void *a, void *b){
    if ( *(int*)a > *(int*)b){
	return 1;
    }else if(*(int*)a == *(int*)b){
	return 0;
    }else{
	return -1;
    }
}

int cmp_dou(void *a, void *b){
    if ( *((double *)a) > *((double *)b) ){
	return 1;
    }else if( *((double *)a) == *((double *)b) ){
	return 0;
    }else{
	return -1;
    }
}

int strBin2Dec(char *str, int *dec){

    /* Converts the string in "str" that suppose to have a continuue
     * sequence of "0" and "1" to its decimal representation.
     */

    int i;
    *dec=0; 
    for(i=strlen(str)-1;i>=0;i--){
        if(str[i]=='0' || str[i]=='1'){
            if(str[i]=='1')
                *dec+=(int)pow(2,strlen(str)-1-i);
        }else{
            return 1;
        }
    }
    return 0;
}

int voigt2mat(double strvoi[6],double strmat[3][3]){

    strmat[0][0]=strvoi[0];strmat[0][1]=strvoi[4];strmat[0][2]=strvoi[5];
    strmat[1][0]=strvoi[4];strmat[1][1]=strvoi[1];strmat[1][2]=strvoi[3];
    strmat[2][0]=strvoi[5];strmat[2][1]=strvoi[3];strmat[2][2]=strvoi[2];
    return 0;
}

double normvec(double * vec, int n){

  double norm;
  
  int    i;

  norm = 0.0;

  for(i=0;i<n;i++){

    norm += pow(vec[i],2);

  }
  norm = sqrt(norm);

  return norm;
 
}
