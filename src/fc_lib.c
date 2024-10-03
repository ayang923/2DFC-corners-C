#include <stdlib.h>
#include <stdio.h>

#include "fc_lib.h"
#include "num_linalg_lib.h"

int read_fc_matrix(MKL_INT d, MKL_INT C, MKL_INT n_r, char* A_filename, char* Q_filename, rd_mat_t *A, rd_mat_t *Q) {
    FILE *A_file = fopen(A_filename, "r");
    FILE *Q_file = fopen(Q_filename, "r");

    if (A_file == NULL || Q_file == NULL) {
        printf("Read failed");
        return 0;
    }

    Q->rows = d;
    Q->columns = d;

    A->rows = C*n_r;
    A->columns = d;

    char c;
    for (int i = 0; i < Q->rows; i++) {
        for (int j = 0; j < Q->columns; j++) {
            MKL_INT idx = sub2ind(Q->rows, Q->columns, (sub_t) {i, j});
            fscanf(Q_file, "%lf%c", Q->mat_data+idx, &c); // Read a double from the file
        }
    }

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
            MKL_INT idx = sub2ind(A->rows, A->columns, (sub_t) {i, j});
            fscanf(A_file, "%lf%c", A->mat_data+idx, &c); // Read a double from the file
        }
    }

    return 1;
}

MKL_INT fc_A_numel(MKL_INT d, MKL_INT C, MKL_INT n_r) {
    return d*C *n_r;
}

MKL_INT fc_Q_numel(MKL_INT d) {
    return d*d;
}

void fcont_gram_blend_S(rd_mat_t fx, MKL_INT d, rd_mat_t A, rd_mat_t Q, rd_mat_t *fcont) {
    // constructs matching point set
    double f_matching_data[d*fx.columns];

    for (int j = 0; j < fx.columns; j++) {
        for (int i = 0; i < d; i++) {
            MKL_INT fx_idx = sub2ind(fx.rows, fx.columns, (sub_t) {i, j});
            MKL_INT f_matching_idx = sub2ind(fx.rows, fx.columns, (sub_t) {i, j});
            f_matching_data[f_matching_idx] = fx_idx;
        }
    }

    fcont->rows = A.rows;
    fcont->columns = fx.columns;

    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, A.rows, fx.columns, d, 1.0, Q.mat_data, A.rows, A.mat_data, fx.columns, 0, fcont->mat_data, A.rows);
}