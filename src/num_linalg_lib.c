#include <stdlib.h>
#include <cblas.h>
#include <assert.h>

#include "num_linalg_lib.h"

size_t sub2ind(size_t rows, size_t columns, sub_t sub) {
    assert(sub.i < rows && sub.j < columns);
    return sub.j*rows + sub.i;
}

sub_t ind2sub(size_t rows, size_t columns, size_t idx) {
    assert(idx < columns*rows);
    return (sub_t) {idx%rows, idx/rows};
}

void rd_ones(rd_mat_t ones_mat) {
    for (size_t i = 0; i < ones_mat.rows * ones_mat.columns; i++) {
        ones_mat.mat_data[i] = 1.0;
    }
}

void rd_linspace(double start, double end, size_t n, rd_mat_t mat_addr) {
    assert(mat_addr.rows == n && mat_addr.columns == 1 && start < end);

    double h = (end-start) / (n-1);
    for (size_t i = 0; i < n; i++) {
        mat_addr.mat_data[i] = start + i*h;
    }
}


