/* ***************************************************

C code for finding all possible combinations with p 
elements from a set with n elements (where p<=n).

Based on Fortran code found somewhere on the web.

S�ren H�jsgaard

** ***************************************************/

#include <Rdefines.h>
#include <R.h>
#include "combnPrimC.h"


void combnC(int *nsel, int *ncand, int *nset, int *ans){

  int ii, ofs, jj, kk;
  int *list;
  //  int ans_len = *nset * *nsel;

  list = (int *) R_alloc(*nsel, sizeof(int));
  for (ii=0;ii<*nsel;ii++)
    list[ii] = 0;

  for (ii=0; ii<*nset; ii++){
    sla_combn__(nsel, ncand, list, &jj); 
    ofs = ii * *nsel;
    for (kk=0;kk<*nsel;kk++)
      ans[ofs+kk] = list[kk];
  }
}
