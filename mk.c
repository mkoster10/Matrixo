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

//remove memory assigent
void mk_mat_free(mk_mat *matrix) {
  int i;
  for(i = 0; i < matrix->num_rows; ++i) {
    free(matrix->mk_data[i]);
  }
  free(matrix->mk_data);
  free(matrix);
}
//create a matrix with random entries
mk_mat *mk_mat_rnd(unsigned int num_rows, unsigned int num_cols, double min, double max) {
  mk_mat *r = mk_mat_new(num_rows, num_cols);
  int i, j;
  for(i = 0; i < num_rows; i++) {
    for(j = 0; j < num_cols; j++) {
      r->mk_data[i][j] = mk_rand_interval(min, max);
    }
  }
  return r;
}
#define	RAND_MAX	0x7fffffff
double mk_rand_interval(double min, double max) {
  double d;
  d = (double) rand() / ((double) RAND_MAX + 1);
  return (min + d * (max - min));
}

mk_mat *mk_mat_cp(mk_mat *m) {
  mk_mat *r  = mk_mat_new(m->num_rows, m->num_cols);
  int i,j;
  for(i = 0; i < r->num_rows; i++) {
    for(j = 0; j < r->num_cols; j++) {
      r->mk_data[i][j] = m->mk_data[i][j];
    }
  }
  return r;
}
//creating a square matrix where c = r 
mk_mat *mk_mat_sqr(unsigned int size) {
  return mk_mat_new(size, size);
}

//returns an identity matrix
mk_mat *mk_mat_eye(unsigned int size) {
  mk_mat *r = mk_mat_new(size, size);
  int i;
  for(i = 0; i < r->num_rows; i++) {
    r->mk_data[i][i] = 1.0;
  }
  return r;
}


//FUNCTIONS CONCERNED WITH EQUALITY


// Checks if two matrices have the same dimesions
int mk_mat_eqdim(mk_mat *m1, mk_mat *m2) {
  return (m1->num_cols == m2->num_cols) &&
          (m1->num_rows == m2->num_rows);
}

// Checks if two matrices have the same dimensions, and the elements
// are all equal to each other with a given tolerance;
// For exact equality use tolerance = 0.0
int mk_mat_eq(mk_mat *m1, mk_mat *m2, double tolerance) {
  if (!mk_mat_eqdim(m1, m2)) {
    return 0;
  }
  int i, j;
  for(i = 0; i < m1->num_rows; i++) {
    for(j = 0; j < m1->num_cols; j++) {
      if (fabs(m1->mk_data[i][j] - m2->mk_data[i][j]) > tolerance) {
        return 0;
      }
    }
  }
  return 1;
}


void mk_mat_print(mk_mat *matrix) {
  mk_mat_printf(matrix, "%lf\t\t");
}
// Prints the matrix on the stdout (with a custom formatting for elements)
void mk_mat_printf(mk_mat *matrix, const char *d_fmt) {
  int i, j;
  fprintf(stdout, "\n");
  for(i = 0; i < matrix->num_rows; ++i) {
    for(j = 0; j < matrix->num_cols; ++j) {
      fprintf(stdout, d_fmt, matrix->mk_data[i][j]);
    }
    fprintf(stdout, "\n");
  }
  fprintf(stdout, "\n");
}
// Accessing/modifying matrix elements
double mk_mat_get(mk_mat *matrix, unsigned int i, unsigned int j){
      return matrix->mk_data[i][j];

}

void mk_mat_set(mk_mat *matrix, unsigned int i, unsigned int j, double value) {
      matrix->mk_data[i][j] = value;

}
mk_mat *mk_mat_col_get(mk_mat *m, unsigned int col) {
  if (col >= m->num_cols) {
    return NULL;
  }
  mk_mat *r = mk_mat_new(m->num_rows, 1);
  int j;
  for(j = 0; j < r->num_rows; j++) {
    r->mk_data[j][0] = m->mk_data[j][col];
  }
  return r;
} 

mk_mat *mk_mat_row_get(mk_mat *m, unsigned int row) {
    if (row >=m->num_rows) {
        return NULL;
    }
    mk_mat *r = mk_mat_new(1,m->num_cols);
    memcpy(r->mk_data[0], m->mk_data[row], m->num_cols * sizeof(*r->mk_data[0]));
  return r;
}

// Sets all elements of a matrix to a given value
void mk_mat_all_set(mk_mat *matrix, double value) {
  int i, j;
  for(i = 0; i < matrix->num_rows; i++) {
    for(j = 0; j < matrix->num_cols; j++) {
      matrix->mk_data[i][j] = value;
    }
  }
}

// Sets all elements of the matrix to given value
int mk_mat_diag_set(mk_mat *m, double value) {
  if (!m->is_square) {
    return 0;
  }
  int i;
  for(i = 0; i < m->num_rows; i++) {
    m->mk_data[i][i] = value;
  }
  return 1;
} 

int mk_mat_row_mult_r(mk_mat *m, unsigned int row, double num) {
  if (row>= m->num_rows) {
    return 0;
  }
  int i;
  for(i=0; i < m->num_cols; i++) {
    m->mk_data[row][i] *= num;
  }
  return 1;
}

mk_mat *mk_mat_row_mult(mk_mat *m, unsigned int row, double num) {
  mk_mat *r = mk_mat_cp(m);
  if (!mk_mat_row_mult_r(r, row, num)) {
    mk_mat_free(r);
    return NULL;
  }
  return r;
} 


mk_mat *mk_mat_row_addrow(mk_mat *m, unsigned int where, 
unsigned int row, double multiplier) {
  mk_mat *r = mk_mat_cp(m);
  if (!mk_mat_row_addrow_r(m, where, row, multiplier)) {
    mk_mat_free(r);
    return NULL;
  }
  return r;
}

int mk_mat_row_addrow_r(mk_mat *m, unsigned int where, 
unsigned int row, double multiplier) {
  if (where >= m->num_rows || row >= m->num_rows) {
    return 0;
  }
  int i = 0;
  for(i = 0; i < m->num_cols; i++) {
    m->mk_data[where][i] += multiplier * m->mk_data[row][i];
  }
  return 1;
} 


mk_mat *mk_mat_smult(mk_mat *m, double num) {
  mk_mat *r = mk_mat_cp(m);
  mk_mat_smult_r(r, num);
  return r;
}

int mk_mat_smult_r(mk_mat *m, double num) {
  int i, j;
  for(i = 0; i < m->num_rows; i++) {
    for(j = 0; j < m->num_cols; j++) {
      m->mk_data[i][j] *= num;
    }
  }
  return 1;
} 

mk_mat *mk_mat_col_rem(mk_mat *m, unsigned int column) {
  if(column >= m->num_cols) {
    return NULL;
  }
  mk_mat *r = mk_mat_new(m->num_rows, m->num_cols-1);
  int i, j, k;
  for(i = 0; i < m->num_rows; i++) {
    for(j = 0, k=0; j < m->num_cols; j++) {
      if (column!=j) {
        r->mk_data[i][k++] = m->mk_data[i][j];
      }
    }
  }
  return r;
}

mk_mat *mk_mat_row_rem(mk_mat *m, unsigned int row) {
  if (row >= m->num_rows) {
    return NULL;
  }
  mk_mat *r = mk_mat_new(m->num_rows-1, m->num_cols);
  int i, j, k;
  for(i = 0, k = 0; i < m->num_rows; i++) {
    if (row!=i) {
      for(j = 0; j < m->num_cols; j++) {
        r->mk_data[k][j] = m->mk_data[i][j];
      }
      k++;
    }
  }
  return r;
}

int mk_mat_row_swap_r(mk_mat *m, unsigned int row1, unsigned int row2) {
  if (row1 >= m->num_rows || row2 >= m->num_rows) {
    return 0;
  }
  double *tmp = m->mk_data[row2];
  m->mk_data[row2] = m->mk_data[row1];
  m->mk_data[row1] = tmp;
  return 1;
} 

mk_mat *mk_mat_row_swap(mk_mat *m, unsigned int row1, unsigned int row2) {
  mk_mat *r = mk_mat_cp(m);
  if (!mk_mat_row_swap_r(r, row1, row2)) {
    mk_mat_free(r);
    return NULL;
  }
  return r;
} 

int mk_mat_col_swap_r(mk_mat *m, unsigned int col1, unsigned int col2) {
  if (col1 >= m->num_cols || col2 >= m->num_rows) {
    return 0;
  }
  double tmp;
  int j;
  for(j = 0; j < m->num_rows; j++) {
    tmp = m->mk_data[j][col1];
    m->mk_data[j][col1] = m->mk_data[j][col2];
    m->mk_data[j][col2] = tmp;
  }
  return 1;
} 

mk_mat *mk_mat_col_swap(mk_mat *m, unsigned int col1, unsigned int col2) {
  mk_mat *r = mk_mat_cp(m);
  if (!mk_mat_col_swap_r(r, col1, col2)) {
    mk_mat_free(r);
    return NULL;
  }
  return r;
} 
mk_mat *mk_mat_add(mk_mat *m1, mk_mat *m2) {
  mk_mat *r = mk_mat_cp(m1);
  if (!mk_mat_add_r(r, m2)) {
    mk_mat_free(r);
    return NULL;
  }
  return r;
}

int mk_mat_add_r(mk_mat *m1, mk_mat *m2) {
  if (!mk_mat_eqdim(m1, m2)) {
    return 0;
  }
  int i, j;
  for(i = 0; i < m1->num_rows; i++) {
    for(j = 0; j < m2->num_rows; j++) {
      m1->mk_data[i][j] += m2->mk_data[i][j];
    }
  }
  return 1;
}

mk_mat *mk_mat_sub(mk_mat *m1, mk_mat *m2) {
  mk_mat *r = mk_mat_cp(m2);
  if (!mk_mat_sub_r(r, m2)) {
    mk_mat_free(r);
    return NULL;
  }
  return r;
}

int mk_mat_sub_r(mk_mat *m1, mk_mat *m2) {
  if (!mk_mat_eqdim(m1, m2)) {
    return 0;
  }
  int i, j;
  for(i = 0; i < m1->num_rows; i++) {
    for(j = 0; j < m1->num_cols; j++) {
      m1->mk_data[i][j] -= m2->mk_data[i][j];
    }
  }
  return 1;
}

mk_mat *mk_mat_dot(mk_mat *m1, mk_mat *m2) {
  if (!(m1->num_cols == m2->num_rows)) {
    return NULL;
  }
  int i, j, k;
  mk_mat *r = mk_mat_new(m1->num_rows, m2->num_cols);
  for(i = 0; i < r->num_rows; i++) {
    for(j = 0; j < r->num_cols; j++) {
      for(k = 0; k < m1->num_cols; k++) {
        r->mk_data[i][j] += m1->mk_data[i][k] * m2->mk_data[k][j];
      }
    }
  }
  return r;
}

mk_mat *mk_mat_transp(mk_mat *m) {
  int i, j;
  mk_mat *r = mk_mat_new(m->num_cols, m->num_rows);
  for(i = 0; i < r->num_rows; i++) {
    for(j = 0; j < r->num_cols; j++) {
      r->mk_data[i][j] = m->mk_data[j][i];
    }
  }
  return r;
}

double mk_mat_trace(mk_mat* m) {
  int i;
  double trace = 0.0;
  for(i = 0; i < m->num_rows; i++) {
    trace += m->mk_data[i][i];
  }
  return trace;
}
// *****************************************************************************
//
// Row Echelon
//
// *****************************************************************************

// Finds the first non-zero element on the col column, under the row row.
// This is used to determine the pivot in gauss Elimination
// If not pivot is found, returns -1
int _mk_mat_pivotidx(mk_mat *m, unsigned int col, unsigned int row) {
  // No validations are made, this is an API Method
  int i;
  for(i = row; i < m->num_rows; i++) {
    if (fabs(m->mk_data[i][col]) > mk_MIN_COEF) {
      return i;
    }
  }
  return -1;
}

// Find the max element from the column "col" under the row "row"
// This is needed to pivot in Gauss-Jordan elimination
// If pivot is not found, return -1
int _mk_mat_pivotmaxidx(mk_mat *m, unsigned int col, unsigned int row) {
  int i, maxi;
  double micol;
  double max = fabs(m->mk_data[row][col]);
  maxi = row;
  for(i = row; i < m->num_rows; i++) {
    micol = fabs(m->mk_data[i][col]);
    if (micol>max) {
      max = micol;
      maxi = i;
    }
  }
  return (max < mk_MIN_COEF) ? -1 : maxi;
}

// Retrieves the matrix in Row Echelon form using Gauss Elimination
mk_mat *mk_mat_ref(mk_mat *m) {
  mk_mat *r = mk_mat_cp(m);
  int i, j, k, pivot;
  j = 0, i = 0;
  while(j < r->num_cols && i < r->num_cols) {
    // Find the pivot - the first non-zero entry in the first column of the matrix
    pivot = _mk_mat_pivotidx(r, j, i);
    if (pivot<0) {
      // All elements on the column are zeros
      // We move to the next column without doing anything
      j++;
      continue;
    }
    // We interchange rows moving the pivot to the first row that doesn't have
    // already a pivot in place
    if (pivot!=i) {
      mk_mat_row_swap_r(r, i, pivot);
    }
    // Multiply each element in the pivot row by the inverse of the pivot
    mk_mat_row_mult_r(r, i, 1/r->mk_data[i][j]);
    // We add multiplies of the pivot so every element on the column equals 0
    for(k = i+1; k < r->num_rows; k++) {
      if (fabs(r->mk_data[k][j]) > mk_MIN_COEF) {
        mk_mat_row_addrow_r(r, k, i, -(r->mk_data[k][j]));
      } 
    }
    i++;
    j++;
  }
  return r;
}

// Retrieves the matrix in Reduced Row Echelon using Guass-Jordan Elimination
mk_mat *mk_mat_rref(mk_mat *m) {
  mk_mat* r = mk_mat_cp(m);
  int i,j,k,pivot;
  i = 0;
  j = 0;
  while(j < r->num_cols && i < r->num_rows) {
    // We find the pivot, the maximum row id (fabs) in the column
    pivot = _mk_mat_pivotmaxidx(r, j, i);
    if (pivot<0) {
      // No pivot, we change columns
      j++;
      continue;
    }
    // We interchange rows to out the pivot row into the 
    // desired position
    if (pivot!=i) {
      mk_mat_row_swap_r(r, i, pivot);
    }
    // We create 1 in the pivot position
    mk_mat_row_mult_r(r, i, 1/r->mk_data[i][j]);
     // We put zeros on the colum with the pivot
    for(k = 0; k < r->num_rows; k++) {
      if (!(k==i)) {
        mk_mat_row_addrow_r(r, k, i, -(r->mk_data[k][j]));
      }
    }
    i++;
    j++;
  }
  return r;
}

// *****************************************************************************
//
// LUP Decomposition
//
// *****************************************************************************

// Finds the maxid on the column (starting from k -> num_rows)
// This method is used for pivoting in LUP decomposition
int _mk_mat_absmaxr(mk_mat *m, unsigned int k) {
  // Find max id on the column;
  int i;
  double max = m->mk_data[k][k];
  int maxIdx = k;
  for(i = k+1; i < m->num_rows; i++) {
    if (fabs(m->mk_data[i][k]) > max) {
      max = fabs(m->mk_data[i][k]);
      maxIdx = i;
    }
  }
  return maxIdx;
}

mk_mat_lup *mk_mat_lup_new(mk_mat *L, mk_mat *U, mk_mat *P, unsigned int num_permutations) {
  mk_mat_lup *r = malloc(sizeof(*r));
  NP_CHECK(r);
  r->L = L;
  r->U = U;
  r->P = P;
  r->num_permutations = num_permutations;
  return r;
}
void mk_mat_lup_free(mk_mat_lup* lu) {
  mk_mat_free(lu->P);
  mk_mat_free(lu->L);
  mk_mat_free(lu->U);
  free(lu);
}

void mk_mat_lup_print(mk_mat_lup *lu) {
  mk_mat_print(lu->L);
  mk_mat_print(lu->U);
  mk_mat_print(lu->P);
}

void mk_mat_lup_printf(mk_mat_lup *lu, const char *fmt) {
  mk_mat_printf(lu->L, fmt);
  mk_mat_printf(lu->U, fmt);
  mk_mat_printf(lu->P, fmt);
}

mk_mat_lup *mk_mat_lup_solve(mk_mat *m) {
  if (!m->is_square) {
    return NULL;
  }
  mk_mat *L = mk_mat_new(m->num_rows, m->num_rows);
  mk_mat *U = mk_mat_cp(m);
  mk_mat *P = mk_mat_eye(m->num_rows);

  int j,i, pivot;
  unsigned int num_permutations = 0;
  double mult;

  for(j = 0; j < U->num_cols; j++) {
    // Retrieves the row with the biggest element for column (j)
    pivot = _mk_mat_absmaxr(U, j);
    if (fabs(U->mk_data[pivot][j]) < mk_MIN_COEF) {
      return NULL;
    }
    if (pivot!=j) {
      // Pivots LU and P accordingly to the rule
      mk_mat_row_swap_r(U, j, pivot);
      mk_mat_row_swap_r(L, j, pivot);
      mk_mat_row_swap_r(P, j, pivot);
      // Keep the number of permutations to easily calculate the
      // determinant sign afterwards
      num_permutations++;
    }
    for(i = j+1; i < U->num_rows; i++) {
      mult = U->mk_data[i][j] / U->mk_data[j][j];
      // Building the U upper rows
      mk_mat_row_addrow_r(U, i, j, -mult);
      // Store the multiplier in L
      L->mk_data[i][j] = mult;
    }
  }
  mk_mat_diag_set(L, 1.0);

  return mk_mat_lup_new(L, U, P, num_permutations);
}

// After the LU(P) factorisation the determinant can be easily calculated
// by multiplying the main diagonal of matrix U with the sign.
// the sign is -1 if the number of permutations is odd
// the sign is +1 if the number of permutations is even
double mk_mat_det(mk_mat_lup* lup) {
  int k;
  int sign = (lup->num_permutations%2==0) ? 1 : -1;
  mk_mat *U = lup->U;
  double product = 1.0;
  for(k = 0; k < U->num_rows; k++) {
    product *= U->mk_data[k][k];
  }
  return product * sign;
}

// Returns LU matrix from a LUP structure
mk_mat *mk_mat_lu_get(mk_mat_lup* lup) {
  mk_mat *r = mk_mat_cp(lup->U);
  // Copy L (without first diagonal in result)
  int i, j;
  for(i = 1; i < lup->L->num_rows; i++) {
    for(j = 0; j < i; j++) {
      r->mk_data[i][j] = lup->L->mk_data[i][j];
    }
  }
  return r;
}

// *****************************************************************************
//
// Solving linear systems of equations
//
// *****************************************************************************

// Forward substitution algorithm
// Solves the linear system L * x = b
//
// L is lower triangular matrix of size NxN
// B is column matrix of size Nx1
// x is the solution column matrix of size Nx1
//
// Note: In case L is not a lower triangular matrix, the algorithm will try to
// select only the lower triangular part of the matrix L and solve the system
// with it.
//
// Note: In case any of the diagonal elements (L[i][i]) are 0 the system cannot
// be solved
//
// Note: This function is usually used with an L matrix from a LU decomposition
mk_mat *mk_ls_solvefwd(mk_mat *L, mk_mat *b) {
  mk_mat* x = mk_mat_new(L->num_cols, 1);
  int i,j;
  double tmp;
  for(i = 0; i < L->num_cols; i++) {
    tmp = b->mk_data[i][0];
    for(j = 0; j < i ; j++) {
      tmp -= L->mk_data[i][j] * x->mk_data[j][0];
    }
    x->mk_data[i][0] = tmp / L->mk_data[i][i];
  }
  return x;
}


// Back substition algorithm
// Solves the linear system U *x = b
//
// U is an upper triangular matrix of size NxN
// B is a column matrix of size Nx1
// x is the solution column matrix of size Nx1
//
// Note in case U is not an upper triangular matrix, the algorithm will try to
// select only the upper triangular part of the matrix U and solve the system
// with it
//
// Note: In case any of the diagonal elements (U[i][i]) are 0 the system cannot
// be solved
mk_mat *mk_ls_solvebck(mk_mat *U, mk_mat *b) {
  mk_mat *x = mk_mat_new(U->num_cols, 1);
  int i = U->num_cols, j;
  double tmp;
  while(i-->0) {
    tmp = b->mk_data[i][0];
    for(j = i; j < U->num_cols; j++) {
      tmp -= U->mk_data[i][j] * x->mk_data[j][0];
    }
    x->mk_data[i][0] = tmp / U->mk_data[i][i];
  }
  return x;
}

// A[n][n] is a square matrix
// m contains matrices L, U, P for A[n][n] so that P*A = L*U
//
// The linear system is:
// A*x=b  =>  P*A*x = P*b  =>  L*U*x = P*b  =>
// (where b is a matrix[n][1], and x is a matrix[n][1])
//
// if y = U*x , we solve two systems:
//    L * y = P b (forward substition)
//    U * x = y (backward substition)
//
// We obtain and return x
mk_mat *mk_ls_solve(mk_mat_lup *lu, mk_mat* b) {
  if (lu->U->num_rows != b->num_rows || b->num_cols != 1) {
      return NULL;
  }
  mk_mat *Pb = mk_mat_dot(lu->P, b);

  // We solve L*y = P*b using forward substition
  mk_mat *y = mk_ls_solvefwd(lu->L, Pb);

  // We solve U*x=y
  mk_mat *x = mk_ls_solvebck(lu->U, y);

  mk_mat_free(y);
  mk_mat_free(Pb);
  return x;
}

// Calculates the inverse of a matrix
mk_mat *mk_mat_inv(mk_mat_lup *lup) {
  unsigned n = lup->L->num_cols;
  mk_mat *r = mk_mat_sqr(n);
  mk_mat *I = mk_mat_eye(lup->U->num_rows);
  mk_mat *invx;
  mk_mat *Ix;
  int i,j;
  for(j =0; j < n; j++) {
    Ix = mk_mat_col_get(I, j);
    invx = mk_mat_solve(lup, Ix);
    for(i = 0; i < invx->num_rows; i++) {
      r->mk_data[i][j] = invx->mk_data[i][0];
    }
    mk_mat_free(invx);
    mk_mat_free(Ix);
  }
  mk_mat_free(I);
  return r;
}


// *****************************************************************************
//
// QR Decomposition
//
// *****************************************************************************

// Useful for QR decomposition
// Represents the (dot) product of two vectors:
// vector1 = m1col column from m1
// vector2 = m2col column from m2
double mk_vect_dot(mk_mat *m1, unsigned int m1col, mk_mat *m2, unsigned m2col) {
  int i;
  double dot = 0.0;
  for(i = 0; i < m1->num_rows; i++) {
    dot += m1->mk_data[i][m1col] * m2->mk_data[i][m2col];
  }
  return dot;
}

// Calculates the l2 norm for a colum in the matrix
double mk_mat_col_l2norm(mk_mat *m, unsigned int col) {
  double doublesum = 0.0;
  int i;
  for(i = 0; i < m->num_rows; i++) {
    doublesum += (m->mk_data[i][col]*m->mk_data[i][col]);
  }
  return sqrt(doublesum);
}

// Calculates the l2norm for each column
// Keeps results into 1 row matrix
mk_mat *mk_mat_l2norm(mk_mat *m) {
  int i, j;
  mk_mat *r = mk_mat_new(1, m->num_cols);
  double square_sum;
  for(j = 0; j < m->num_cols; j++) {
    square_sum = 0.0;
    for(i = 0; i < m->num_rows; i++) {
      square_sum+=m->mk_data[i][j]*m->mk_data[i][j];
    }
    r->mk_data[0][j] = sqrt(square_sum);
  }
  return r;
}

mk_mat *mk_mat_normalize(mk_mat *m) {
  mk_mat *r = mk_mat_cp(m);
  if (!mk_mat_mat_normalize_r(r)) {
    mk_mat_free(r);
    return NULL;
  }
  return r;
}

int mk_mat_normalize_r(mk_mat *m) {
  mk_mat *l2norms = mk_mat_l2norm(m);
  int j;
  for(j = 0; j < m->num_cols; j++) {
    if (l2norms->mk_data[0][j] < mk_MIN_COEF) {
      mk_mat_free(l2norms);
      return 0;
    }
    mk_mat_col_mult_r(m, j, 1/l2norms->mk_data[0][j]);
  }
  mk_mat_free(l2norms);
  return 1;
}

mk_mat_qr *mk_mat_qr_new() {
  mk_mat_qr *qr = malloc(sizeof(*qr));
  NP_CHECK(qr);
  return qr;
}

void mk_mat_qr_free(mk_mat_qr *qr) {
  mk_mat_free(qr->Q);
  mk_mat_free(qr->R);
  free(qr);
}

// M = QR
mk_mat_qr *mk_mat_qr_solve(mk_mat *m) {

  mk_mat_qr *qr = mk_mat_qr_new();
  mk_mat *Q = mk_mat_cp(m);
  mk_mat *R = mk_mat_new(m->num_rows, m->num_cols);

  int j, k;
  double l2norm;
  double rkj;
  mk_mat *aj;
  mk_mat *qk;
  for(j=0; j < m->num_cols; j++) {    
    rkj = 0.0;
    aj = mk_mat_col_get(m, j);
    for(k = 0; k < j; k++) {
       rkj = mk_vect_dot(m, j, Q, k);
       R->mk_data[k][j] = rkj;
       qk = mk_mat_col_get(Q, k);
       mk_mat_col_mult_r(qk, 0, rkj);
       mk_mat_sub_r(aj, qk);
       mk_mat_free(qk);
    }
    for(k = 0; k < Q->num_rows; k++) {
      Q->mk_data[k][j] = aj->mk_data[k][0];
    }
    l2norm = mk_mat_col_l2norm(Q, j);
    mk_mat_col_mult_r(Q, j, 1/l2norm);
    R->mk_data[j][j] = l2norm;
    mk_mat_free(aj);
  }
  qr->Q = Q;
  qr->R = R;
  return qr;
}
