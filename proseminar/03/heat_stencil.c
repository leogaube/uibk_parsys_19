#include <stdio.h>
#include <stdlib.h>

#include "heat_stencil.h"


Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }
