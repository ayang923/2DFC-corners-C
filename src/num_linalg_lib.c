#include <stdlib.h>
#include <mkl.h>
#include <assert.h>

#include "num_linalg_lib.h"

MKL_INT sub2ind(MKL_INT rows, MKL_INT columns, sub_t sub) {
    assert(sub.i < rows && sub.j < columns);
    return sub.j*rows + sub.i;
}

sub_t ind2sub(MKL_INT rows, MKL_INT columns, MKL_INT idx) {
    assert(idx < columns*rows);
    return (sub_t) {idx%rows, idx/rows};
}

void rd_linspace(double start, double end, MKL_INT n, rd_mat_t mat_addr) {
    assert(mat_addr.rows == n && mat_addr.columns == 1 && start < end);

    double h = (end-start) / (n-1);
    for (MKL_INT i = 0; i < n; i++) {
        mat_addr.mat_data[i] = start + i*h;
    }
}

void rd_meshgrid(rd_mat_t x, rd_mat_t y, rd_mat_t *X, rd_mat_t *Y) {
    // assumes x and y are column vectors

    X->rows = y.rows;
    X->columns = x.rows;
    Y->rows = y.rows;
    Y->columns = x.rows;

    for (int i = 0; i < y.rows; i++) {
        for (int j = 0; j < x.rows; j++) {
            MKL_INT idx = sub2ind(y.rows, x.rows, (sub_t) {i, j});
            X->mat_data[idx] = x.mat_data[j];
            Y->mat_data[idx] = y.mat_data[i];
        }
    }
}


