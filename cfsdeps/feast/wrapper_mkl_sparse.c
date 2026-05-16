/* wrapper_mkl_sparse.c - C shim so Fortran can call Inspector-Executor API */
#include "mkl_spblas.h"

/* Fortran calls these with trailing underscore - we forward to C API */
void mkl_sparse_d_create_csr_(
    sparse_matrix_t *A, sparse_index_base_t *indexing,
    int *rows, int *cols,
    int *rows_start, int *rows_end, int *col_indx,
    double *values, sparse_status_t *stat)
{
    *stat = mkl_sparse_d_create_csr(A, *indexing, *rows, *cols,
                                     rows_start, rows_end, col_indx, values);
}

void mkl_sparse_d_mm_(
    sparse_operation_t *operation, double *alpha,
    sparse_matrix_t *A, struct matrix_descr *descr,
    sparse_layout_t *layout, double *x, int *columns,
    int *ldx, double *beta, double *y, int *ldy,
    sparse_status_t *stat)
{
    *stat = mkl_sparse_d_mm(*operation, *alpha, *A, *descr,
                             *layout, x, *columns, *ldx, *beta, y, *ldy);
}

void mkl_sparse_destroy_(sparse_matrix_t *A, sparse_status_t *stat)
{
    *stat = mkl_sparse_destroy(*A);
}

/* Complex double - wzcsrmm, wzhcsrmm */
void mkl_sparse_z_create_csr_(
    void **A, int *indexing,
    int *rows, int *cols,
    int *rows_start, int *rows_end, int *col_indx,
    void *values, int *stat)
{
    *stat = mkl_sparse_z_create_csr(
        (sparse_matrix_t*)A, *indexing, *rows, *cols,
        rows_start, rows_end, col_indx, (MKL_Complex16*)values);
}

void mkl_sparse_z_mm_(
    int *operation, void *alpha,
    void **A, int *descr,
    int *layout, void *x, int *columns,
    int *ldx, void *beta, void *y, int *ldy,
    int *stat)
{
    struct matrix_descr d;
    d.type = descr[0]; d.mode = descr[1]; d.diag = descr[2];
    MKL_Complex16 al = *(MKL_Complex16*)alpha;
    MKL_Complex16 be = *(MKL_Complex16*)beta;
    *stat = mkl_sparse_z_mm(*operation, al, *(sparse_matrix_t*)A, d,
                             *layout, (MKL_Complex16*)x, *columns,
                             *ldx, be, (MKL_Complex16*)y, *ldy);
}

/* Real single - wscsrmm */
void mkl_sparse_s_create_csr_(
    void **A, int *indexing,
    int *rows, int *cols,
    int *rows_start, int *rows_end, int *col_indx,
    float *values, int *stat)
{
    *stat = mkl_sparse_s_create_csr(
        (sparse_matrix_t*)A, *indexing, *rows, *cols,
        rows_start, rows_end, col_indx, values);
}

void mkl_sparse_s_mm_(
    int *operation, float *alpha,
    void **A, int *descr,
    int *layout, float *x, int *columns,
    int *ldx, float *beta, float *y, int *ldy,
    int *stat)
{
    struct matrix_descr d;
    d.type = descr[0]; d.mode = descr[1]; d.diag = descr[2];
    *stat = mkl_sparse_s_mm(*operation, *alpha, *(sparse_matrix_t*)A, d,
                             *layout, x, *columns, *ldx, *beta, y, *ldy);
}

/* Complex single - wccsrmm, wchcsrmm */
void mkl_sparse_c_create_csr_(
    void **A, int *indexing,
    int *rows, int *cols,
    int *rows_start, int *rows_end, int *col_indx,
    void *values, int *stat)
{
    *stat = mkl_sparse_c_create_csr(
        (sparse_matrix_t*)A, *indexing, *rows, *cols,
        rows_start, rows_end, col_indx, (MKL_Complex8*)values);
}

void mkl_sparse_c_mm_(
    int *operation, void *alpha,
    void **A, int *descr,
    int *layout, void *x, int *columns,
    int *ldx, void *beta, void *y, int *ldy,
    int *stat)
{
    struct matrix_descr d;
    d.type = descr[0]; d.mode = descr[1]; d.diag = descr[2];
    MKL_Complex8 al = *(MKL_Complex8*)alpha;
    MKL_Complex8 be = *(MKL_Complex8*)beta;
    *stat = mkl_sparse_c_mm(*operation, al, *(sparse_matrix_t*)A, d,
                             *layout, (MKL_Complex8*)x, *columns,
                             *ldx, be, (MKL_Complex8*)y, *ldy);
}
