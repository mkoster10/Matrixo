#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>

#include "mk.h"
mk_mat *mk_mat_new(unsigned int num_rows, unsigned int num_cols) {
  if (num_rows == 0) {
    return NULL;
  }
  if (num_cols == 0) {
    return NULL;
  }
  mk_mat *m = calloc(1, sizeof(*m));
  m->num_rows = num_rows;
  m->num_cols = num_cols;
  m->is_square = (num_rows == num_cols) ? 1 : 0;
  m->mk_data = calloc(m->num_rows, sizeof(*m->mk_data));
  int i;
  for(i = 0; i < m->num_rows; ++i) {
    m->mk_data[i] = calloc(m->num_cols, sizeof(**m->mk_data));
  }
  return m;
}


void mk_mat_free(mk_mat *matrix) {
  int i;
  for(i = 0; i < matrix->num_rows; ++i) {
    free(matrix->mk_data[i]);
  }
  free(matrix->mk_data);
  free(matrix);
}

mk_mat *nml_mat_rnd(unsigned int num_rows, unsigned int num_cols, double min, double max) {
  mk_mat *r = nml_mat_new(num_rows, num_cols);
  int i, j;
  for(i = 0; i < num_rows; i++) {
    for(j = 0; j < num_cols; j++) {
      r->mk_data[i][j] = nml_rand_interval(min, max);
    }
  }
  return r;
}

