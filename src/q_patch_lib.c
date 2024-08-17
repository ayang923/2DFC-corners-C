#include "q_patch_lib.h"
#include <stddef.h>

void q_patch_init(q_patch_t *q_patch, size_t n_xi, size_t n_eta, double xi_start, double xi_end, double eta_start, double eta_end) {
    q_patch->n_xi = n_xi;
    q_patch->n_eta = n_eta;
    q_patch->xi_start = xi_start;
    q_patch->xi_end = xi_end;
    q_patch->eta_start = eta_start;
    q_patch->eta_end = eta_end;
}