#ifndef __NUM_LINALG_LIB__
#define __NUM_LINALG_LIB__

#include <cblas.h>
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
    size_t rows;
    size_t columns;
} rd_mat_t;

/**
 * @brief subscript struct
 * 
 * @param i: row number
 * @param j: column number
 */
typedef struct sub {
    size_t i;
    size_t j;
} sub_t;

/**
 * @brief Converts subscript to index
 * 
 * @param rows: # rows in matrix
 * @param columns: # col in matrix
 * @param sub: subscript
 * @return size_t: index
 */
size_t sub2ind(size_t rows, size_t columns, sub_t sub);

/**
 * @brief Converts subscript to index
 * 
 * @param rows: # rows in matrix
 * @param columns: # col in matrix
 * @param idx: index
 * @return sub_t: subscript
 */
sub_t ind2sub(size_t rows, size_t columns, size_t idx);

/**
 * @brief Creates vector of all ones
 * 
 * @param ones_vector matrix address to store
 */
void rd_ones(rd_mat_t ones_mat);

/**
 * @brief Implementation of matlab linspace using double
 * 
 * @param start 
 * @param end 
 * @param n 
 * @param mat_addr matrix we fill
 */
void rd_linspace(double start, double end, size_t n, rd_mat_t mat_addr);

#endif