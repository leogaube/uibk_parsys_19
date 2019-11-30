#include <stdio.h>
#include <stdlib.h>

#include "heat_stencil.h"


Vector createVector(int N) {
  // create data and index vector
  return malloc(sizeof(value_t) * N);
}

void releaseVector(Vector m) { free(m); }

void printTemperature(Vector m, int nx, int ny, int nz)
{
#ifdef SLICE_3D
  if (nz != 1){
    nz=1;
    printf("########### ONLY DISPLAYING SLICE 0 of 3D room ##############\n");
  }
#endif

  const char *colors = " .-:=+*^X#%@";
  const int numColors = 12;

  // boundaries for temperature (for simplicity hard-coded)
  const value_t max = 273 + 30;
  const value_t min = 273 + 0;

  // set the 'render' resolution
  int Wx = MIN(RESOLUTION, nx);
  int Wy = MIN(RESOLUTION, ny);
  int Wz = MIN(RESOLUTION, nz);

  // step size in each dimension
  int sWx = nx / Wx;
  int sWy = ny / Wy;
  int sWz = nz / Wz;

  // room
  // actual room
  for (int k = 0; k < Wz; k++)
  {
    if (Wz != 1)
      printf("########### SLICE %d ##############\n", k);
    for (int j = 0; j < Wy; j++)
    {
      // left wall
      printf("X");
      for (int i = 0; i < Wx; i++)
      {
        // get max temperature in this tile
        value_t max_t = 0;
        for (int z = sWz * k; z < sWz * k + sWz; z++)
        {
          for (int y = sWy * j; y < sWy * j + sWy; y++)
          {
            for (int x = sWx * i; x < sWx * i + sWx; x++)
            {
              max_t = (max_t < m[IDX_3D(x, y, z, nx, ny)]) ? m[IDX_3D(x, y, z, nx, ny)] : max_t;
            }
          }
        }
        value_t temp = max_t;

        // pick the 'color'
        int c = ((temp - min) / (max - min)) * numColors;
        c = (c >= numColors) ? numColors - 1 : ((c < 0) ? 0 : c);

        // print the average temperature
        // if numbers are desired use printf("%2.2f\t",temp-273);
        printf("%c", colors[c]);
      }
      // right wall
      printf("X\n");
    }
  }
  printf("\n##################################\n");
}
