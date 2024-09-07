#ifndef __NUM_LINALG_LIB__
#define __NUM_LINALG_LIB__

#include <mkl.h>
#include <stdlib.h>

/**
 * @brief Real, double matrix
 * 
 * @param mat_data: array that contains data organized column major
 * @param rows: number of rows
 * @param columns: number of columns
 */
typedef struct rd_mat {
    double *mat_data;
    MKL_INT rows;
    MKL_INT columns;
} rd_mat_t;

/**
 * @brief subscript struct
 * 
 * @param i: row number
 * @param j: column number
 */
typedef struct sub {
    MKL_INT i;
    MKL_INT j;
} sub_t;

/**
 * @brief Converts subscript to index
 * 
 * @param rows: # rows in matrix
 * @param columns: # col in matrix
 * @param sub: subscript
 * @return MKL_INT: index
 */
MKL_INT sub2ind(MKL_INT rows, MKL_INT columns, sub_t sub);

/**
 * @brief Converts subscript to index
 * 
 * @param rows: # rows in matrix
 * @param columns: # col in matrix
 * @param idx: index
 * @return sub_t: subscript
 */
sub_t ind2sub(MKL_INT rows, MKL_INT columns, MKL_INT idx);

/**
 * @brief Implementation of matlab linspace using double
 * 
 * @param start 
 * @param end 
 * @param n 
 * @param mat_addr matrix we fill
 */
void rd_linspace(double start, double end, MKL_INT n, rd_mat_t mat_addr);

void rd_meshgrid(rd_mat_t x, rd_mat_t y, rd_mat_t *X, rd_mat_t *Y);
#endif