#ifdef __cplusplus
extern "C" {
#endif

#ifndef mk_H_
#define mk_H_

#define mk_MIN_COEF 0.000000000000001


// *****************************************************************************
//
// Library structures
//
// *****************************************************************************
typedef struct mk_mat_s {
  unsigned int num_rows;
  unsigned int num_cols;
  double **mk_data;
  int is_square;
} mk_mat;

typedef struct mk_mat_lup_s {
  mk_mat *L;
  mk_mat *U;
  mk_mat *P;
  unsigned int num_permutations;
} mk_mat_lup;

typedef struct mk_mat_qr_s {
  mk_mat *Q;
  mk_mat *R;
} mk_mat_qr;




// *****************************************************************************
//
// Constructing and destroying a matrix struct
//
// *****************************************************************************
mk_mat *mk_mat_new(unsigned int num_rows, unsigned int num_cols);
mk_mat *mk_mat_rnd(unsigned int num_rows, unsigned int num_cols, double min, double max);
mk_mat *mk_mat_sqr(unsigned int size);
mk_mat *mk_mat_eye(unsigned int size);
mk_mat *mk_mat_cp(mk_mat *m);

// *****************************************************************************
//
// Matrix Equality
//
// *****************************************************************************
int mk_mat_eqdim(mk_mat *m1, mk_mat *m2);
int mk_mat_eq(mk_mat *m1, mk_mat *m2, double tolerance);

// *****************************************************************************
//
// Accessing and modifying matrix elements
//
// *****************************************************************************
double mk_mat_get(mk_mat *matrix, unsigned int i, unsigned int j);
void mk_mat_set(mk_mat *matrix, unsigned int i, unsigned int j, double value);
mk_mat *mk_mat_col_get(mk_mat *m, unsigned int col);
mk_mat *mk_mat_col_mult(mk_mat *m, unsigned int col, double num);
int mk_mat_col_mult_r(mk_mat *m, unsigned int col, double num);
mk_mat *mk_mat_row_get(mk_mat *m, unsigned int row);
mk_mat *mk_mat_row_mult(mk_mat *m, unsigned int row, double num);
int mk_mat_row_mult_r(mk_mat *m, unsigned int row, double num);
void mk_mat_all_set(mk_mat *matrix, double value);
int mk_mat_diag_set(mk_mat *matrix, double value);
mk_mat *mk_mat_row_addrow(mk_mat *m, unsigned int where, unsigned int row, double multiplier);
int mk_mat_row_addrow_r(mk_mat *m, unsigned int where, unsigned int row, double multiplier);
mk_mat *mk_mat_smult(mk_mat *m, double num);
int mk_mat_smult_r(mk_mat *m, double num);
// *****************************************************************************
//
// Modifying the matrix structure
//
// *****************************************************************************
mk_mat *mk_mat_col_rem(mk_mat *m, unsigned int column);
mk_mat *mk_mat_row_rem(mk_mat *m, unsigned int row);
mk_mat *mk_mat_row_swap(mk_mat *m, unsigned int row1, unsigned int row2);
int mk_mat_row_swap_r(mk_mat *m, unsigned int row1, unsigned int row2);
mk_mat *mk_mat_col_swap(mk_mat *m, unsigned int col1, unsigned int col2);
int mk_mat_col_swap_r(mk_mat *m, unsigned int col1, unsigned int col2);

// *****************************************************************************
//
// Matrix Operations
//
// *****************************************************************************
mk_mat *mk_mat_add(mk_mat *m1, mk_mat *m2);
int mk_mat_add_r(mk_mat *m1, mk_mat *m2);
mk_mat *mk_mat_sub(mk_mat *m1, mk_mat *m2);
int mk_mat_sub_r(mk_mat *m1, mk_mat *m2);
mk_mat *mk_mat_dot(mk_mat *m1, mk_mat *m2);
mk_mat *mk_mat_transp(mk_mat *m);
double mk_mat_trace(mk_mat* m);
// *****************************************************************************
//
// Row Echelon
//
// *****************************************************************************
mk_mat *mk_mat_ref(mk_mat *m);
mk_mat *mk_mat_rref(mk_mat *m);

// *****************************************************************************
//
// LUP Decomposition
//
// *****************************************************************************

mk_mat_lup *mk_mat_lup_new(mk_mat *L, mk_mat *U, mk_mat *P, unsigned int num_permutations);
mk_mat_lup *mk_mat_lup_solve(mk_mat *m);
void mk_mat_lup_free(mk_mat_lup* lu);
void mk_mat_lup_print(mk_mat_lup *lu);
void mk_mat_lup_printf(mk_mat_lup *lu, const char *fmt);
double mk_mat_det(mk_mat_lup* lup);
mk_mat *mk_mat_lu_get(mk_mat_lup* lup);
mk_mat *mk_mat_inv(mk_mat_lup *m);
//*****************************************************************************
//
// Solving linear systems of equations
//
// *****************************************************************************

mk_mat *mk_ls_solvefwd(mk_mat *low_triang, mk_mat *b);
mk_mat *mk_ls_solvebck(mk_mat *upper_triang, mk_mat *b);
mk_mat *mk_ls_solve(mk_mat_lup *lup, mk_mat* b);

// *****************************************************************************
//
// QR Decomposition
//
// *****************************************************************************

double mk_vect_dot(mk_mat *m1, unsigned int m1col, mk_mat *m2, unsigned m2col);
mk_mat *mk_mat_l2norm(mk_mat *m);
double mk_mat_col_l2norm(mk_mat *m1, unsigned int j);
mk_mat *mk_mat_normalize(mk_mat *m);
int mk_mat_normalize_r(mk_mat *m);
mk_mat_qr *mk_mat_qr_new();
void mk_mat_qr_free(mk_mat_qr *qr);
mk_mat_qr * mk_mat_qr_solve(mk_mat *m);

#endif
#ifdef __cplusplus
}
#endif


