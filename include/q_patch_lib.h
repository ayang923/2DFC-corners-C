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
 * @param extra_param: allows us to put in a struct if needed for extra parameters
 */
typedef void (*M_p_t) (rd_mat_t xi, rd_mat_t eta, rd_mat_t x, rd_mat_t y, void* extra_param);

/**
 * @brief phi handle type
 * 
 * @param xi: vector of xi values
 * @param eta: vector of eta values
 * @param phi_vals: address of where to store phi values
 * @param extra_param: allows us to put in a struct if needed for extra parameters
 */
typedef void (*phi_t) (rd_mat_t xi, rd_mat_t eta, rd_mat_t phi_vals, void* extra_param);

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
 * @param extra_param: allows us to put in a struct if needed for extra parameters
 */
typedef void (*J_t) (rd_mat_t v, rd_mat_t J_vals, void* extra_param);

typedef struct phi_param {
    double R_xi;
    double R_eta;
    double xi_0;
    double eta_0;
} phi_param_t;

/**
 * @brief q_patch type struct
 * 
 */
typedef struct q_patch {
    M_p_t M_p;
    J_t J;

    MKL_INT n_xi;
    MKL_INT n_eta;
    double xi_start;
    double xi_end;
    double eta_start;
    double eta_end;

    rd_mat_t f_XY;
    phi_1D_t phi_1D;
    phi_param_t phi_param;

    double eps_xi_eta;
    double eps_xy;
} q_patch_t;

/**
 * @brief Initialize preallocated q_patch_t
 * 
 * @param q_patch 
 * @param M_p 
 * @param J 
 * @param eps_xi_eta 
 * @param eps_xy 
 * @param n_xi 
 * @param n_eta 
 * @param xi_start 
 * @param xi_end 
 * @param eta_start 
 * @param eta_end 
 * @param f_XY 
 * @param phi 
 */
void q_patch_init(q_patch_t *q_patch, M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t f_XY, phi_param_t phi_param);

/**
 * @brief wrapper function for evaluating M_p
 * 
 */
void evaulate_M_p(q_patch_t *q_patch, rd_mat_t xi, rd_mat_t eta, rd_mat_t x, rd_mat_t y);

/**
 * @brief wrapper function for evaluating J
 * 
 * @param q_patch 
 * @param v 
 * @param J_vals 
 */
void evaulate_J(q_patch_t *q_patch, rd_mat_t v, rd_mat_t J_vals);

phi_1D_t return_phi_1D(void);

#endif