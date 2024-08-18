#ifndef __Q_PATCH_LIB_H__
#define __Q_PATCH_LIB_H__

#include <stdlib.h>
#include "num_linalg_lib.h"

/**
 * @brief M_p function handle type
 * 
 * @param xi: vector of xi values
 * @param eta: vector eta values
 * @param x: address of where to store x values
 * @param y: address of where to store y values
 */
typedef void (*M_p_t) (rd_mat_t xi, rd_mat_t eta, rd_mat_t x, rd_mat_t y);

/**
 * @brief phi handle type
 * 
 * @param xi: vector of xi values
 * @param eta: vector of eta values
 * @param phi_vals: address of where to store phi values
 */
typedef void (*phi_t) (rd_mat_t xi, rd_mat_t eta, rd_mat_t phi_vals);

/**
 * @brief phi_1D handle type
 * 
 * @param x: vector of 1D input
 * @param phi_1D_vals: address of where to store output
 */
typedef void (*phi_1D_t) (rd_mat_t x, rd_mat_t phi_1D_vals);

/**
 * @brief Jacobian handle type
 * 
 * @param v: 2 element vector containing xi and eta value
 * @param J_vals: address of where to store J_vals
 */
typedef void (*J_t) (rd_mat_t v, rd_mat_t J_vals);

/**
 * @brief q_patch type struct
 * 
 */
typedef struct q_patch {
    M_p_t M_p;
    J_t J;

    size_t n_xi;
    size_t n_eta;
    double xi_start;
    double xi_end;
    double eta_start;
    double eta_end;

    rd_mat_t f_XY;
    phi_1D_t phi_1D;
    phi_t phi;

    double eps_xi_eta;
    double eps_xy;
} q_patch_t;


void q_patch_init(q_patch_t *q_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, size_t n_xi, size_t n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t f_XY, phi_t phi);

phi_1D_t return_phi_1D(void);

#endif