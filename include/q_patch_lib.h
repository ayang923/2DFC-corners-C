#ifndef __Q_PATCH_LIB_H__
#define __Q_PATCH_LIB_H__

#include <stddef.h>

typedef struct q_patch {
    size_t n_xi;
    size_t n_eta;
    double xi_start;
    double xi_end;
    double eta_start;
    double eta_end;
} q_patch_t;

void q_patch_init(q_patch_t *q_patch, size_t n_xi, size_t n_eta, double xi_start, double xi_end, double eta_start, double eta_end);

#endif