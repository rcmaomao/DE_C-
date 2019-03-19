# include "de.h"

float evaluate(int D, float tmp[], long *nfeval) {
  // polynomial fitting problem
  int i, j;
  int const M = 60;
  float px, x = -1, dx = M, result = 0;

  (*nfeval)++;
  dx = 2/dx;
  for (i=0; i<=M; i++) {
    px = tmp[0];
    for (j=1; j<D; j++) {
      px = x*px + tmp[j];
    }
    if (px<-1 || px>1) {
      result += (1 - px)*(1 - px);
    }
    x+=dx;
  }
  px = tmp[0];
  for (j=1; j<D; j++) {
    px = (float) 1.2*px+tmp[j];
  }
  px = (float) (px - 72.661);
  if (px<0) {
    result += px*px;
  }
  px = tmp[0];
  for (j=1; j<D; j++) {
    px = (float) (-1.2*px + tmp[j]);
  }
  px = (float) (px-72.661);
  if (px<0) {
    result += px*px;
  }
  // printf("%2d %15.4f\n", i, result);
  return result;
}