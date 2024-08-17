#include <stdio.h>
#include "q_patch_lib.h"

int main() {
    printf("hello world!\n");
    q_patch_t q_patch;

    q_patch_init(&q_patch, 20, 20, 0, 1, 0, 1);

    printf("n_xi: %d\n", q_patch.n_xi);
    printf("n_eta: %d\n", q_patch.n_eta);
    printf("xi_start: %f\n", q_patch.xi_start);
    printf("xi_end: %f\n", q_patch.xi_end);
    printf("eta_start: %f\n", q_patch.eta_start);
    printf("eta_end: %f\n", q_patch.eta_end);
}