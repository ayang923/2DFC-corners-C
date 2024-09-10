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
typedef void (*M_p_handle_t) (rd_mat_t xi, rd_mat_t eta, rd_mat_t *x, rd_mat_t *y, void* extra_param);

/**
 * @brief Info required to evaluate M_p type function
 */
typedef struct M_p {
    M_p_handle_t M_p_handle;
    void* extra_param;
} M_p_t;

/**
 * @brief Jacobian handle type
 * 
 * @param v: 2 element vector containing xi and eta value
 * @param J_vals: address of where to store J_vals
 * @param extra_param: allows us to put in a struct if needed for extra parameters
 */
typedef void (*J_handle_t) (rd_mat_t v, rd_mat_t *J_vals, void* extra_param);

typedef struct J {
    J_handle_t J_handle;
    void* extra_param;
} J_t;

/**
 * @brief phi handle type
 * 
 * @param xi: vector of xi values
 * @param eta: vector of eta values
 * @param phi_vals: address of where to store phi values
 * @param extra_param: allows us to put in a struct if needed for extra parameters
 */
typedef void (*phi_handle_t) (rd_mat_t xi, rd_mat_t eta, rd_mat_t *phi_vals, void* extra_param);

/**
 * @brief phi_1D handle type
 * 
 * @param x: vector of 1D input
 * @param phi_1D_vals: address of where to store output
 */
typedef void (*phi_1D_t) (rd_mat_t x, rd_mat_t *phi_1D_vals);

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

    rd_mat_t *f_XY;
    phi_1D_t phi_1D;
    void* phi_param;

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
q_patch_t q_patch_init(M_p_t M_p, J_t J, double eps_xi_eta, double eps_xy, MKL_INT n_xi, MKL_INT n_eta, double xi_start, double xi_end, double eta_start, double eta_end, rd_mat_t *f_XY, void* phi_param);

MKL_INT q_patch_grid_num_el(q_patch_t *q_patch);

/**
 * @brief wrapper function for evaluating M_p
 * 
 */
void q_patch_evaulate_M_p(q_patch_t *q_patch, rd_mat_t xi, rd_mat_t eta, rd_mat_t *x, rd_mat_t *y);

/**
 * @brief wrapper function for evaluating J
 * 
 * @param q_patch 
 * @param v 
 * @param J_vals 
 */
void q_patch_evaulate_J(q_patch_t *q_patch, rd_mat_t v, rd_mat_t *J_vals);

void q_patch_xi_mesh(q_patch_t *q_patch, rd_mat_t *xi_mesh_vals);

void q_patch_eta_mesh(q_patch_t *q_patch, rd_mat_t *eta_mesh_vals);

void q_patch_xi_eta_mesh(q_patch_t *q_patch, rd_mat_t *XI_vals, rd_mat_t *ETA_vals);

void q_patch_convert_to_XY(q_patch_t *q_patch, rd_mat_t XI, rd_mat_t ETA, rd_mat_t *X_vals, rd_mat_t *Y_vals);

void q_patch_xy_mesh(q_patch_t *q_patch, rd_mat_t *X_vals, rd_mat_t *Y_vals);

phi_1D_t return_phi_1D(void);

#endif