#ifdef __cplusplus
extern "C" {
#endif

#ifndef mk_H_
#define mk_H_



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
#endif
#ifdef __cplusplus
}
#endif
