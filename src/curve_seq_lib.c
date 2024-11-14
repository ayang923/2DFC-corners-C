#include <stdlib.h>
#include <mkl.h>
#include <math.h>

#include "curve_seq_lib.h"
#include "num_linalg_lib.h"
#include "s_patch_lib.h"

void M_p_S_general_handle(rd_mat_t xi, rd_mat_t eta, double H, rd_mat_t *x, rd_mat_t *y, void* extra_param) {
    rd_mat_shape(x, xi.rows, xi.columns);
    rd_mat_shape(y, xi.rows, xi.columns);
    
    // using p for conciseness
    M_p_S_general_param_t *p = (M_p_S_general_param_t*) extra_param;

    for (MKL_INT i = 0; i < xi.rows*xi.columns; i++) {
        double xi_val = xi.mat_data[i];
        double eta_val = eta.mat_data[i];

        double xi_tilde = p->xi_diff*xi_val + p->xi_0;
        double nu_norm = sqrt(pow(p->l_1_prime(xi_tilde), 2) + pow(p->l_2_prime(xi_tilde), 2));

        x->mat_data[i] = p->l_1(xi_tilde) - eta_val*H*p->l_2_prime(xi_tilde)/nu_norm;
        y->mat_data[i] = p->l_2(xi_tilde) + eta_val*H*p->l_1_prime(xi_tilde)/nu_norm;
    }
}

void J_S_general_handle(rd_mat_t v, double H, rd_mat_t *J_vals, void* extra_param) {
    rd_mat_shape(J_vals, 2, 2);

    J_S_general_param_t *p = (J_S_general_param_t*) extra_param;

    double xi = v.mat_data[0];
    double eta = v.mat_data[1];

    double xi_tilde = p->xi_diff*xi + p->xi_0;
    double nu_norm = sqrt(pow(p->l_1_prime(xi_tilde), 2) + pow(p->l_2_prime(xi_tilde), 2));

    J_vals->mat_data[0] = p->xi_diff * (p->l_1_prime(xi_tilde)-eta*H*(p->l_2_dprime(xi_tilde)*pow(nu_norm, 2)-p->l_2_prime(xi_tilde)*(p->l_2_dprime(xi_tilde)*p->l_2_prime(xi_tilde)+p->l_1_dprime(xi_tilde)*p->l_1_prime(xi_tilde)))/pow(nu_norm, 3));
    J_vals->mat_data[1] = p->xi_diff * (p->l_2_prime(xi_tilde)+eta*H*(p->l_1_dprime(xi_tilde)*pow(nu_norm, 2)-p->l_1_prime(xi_tilde)*(p->l_2_dprime(xi_tilde)*p->l_2_prime(xi_tilde)+p->l_1_dprime(xi_tilde)*p->l_1_prime(xi_tilde)))/pow(nu_norm, 3));
    J_vals->mat_data[2] = -H*p->l_2_prime(xi_tilde) / nu_norm;
    J_vals->mat_data[3] = H*p->l_1_prime(xi_tilde) / nu_norm;
}

void curve_seq_init(curve_seq_t *curve_seq) {
    curve_seq->first_curve = NULL;
    curve_seq->last_curve = NULL;
    curve_seq->n_curves = 0;
}

void curve_init(curve_t *curve, scalar_func_t l_1, scalar_func_t l_2, scalar_func_t l_1_prime, scalar_func_t l_2_prime, scalar_func_t l_1_dprime, scalar_func_t l_2_dprime, MKL_INT n, double frac_n_C_0, double frac_n_C_1, double frac_n_S_0, double frac_n_S_1, double h_norm, curve_t *next_curve) {
    curve->l_1 = l_1;
    curve->l_2 = l_2;
    curve->l_1_prime = l_1_prime;
    curve->l_2_prime = l_2_prime;
    curve->l_1_dprime = l_1_dprime;
    curve->l_2_dprime = l_2_dprime;

    curve->n = n;

    if(frac_n_C_0 == 0) {
        curve->n_C_0 = ceil(1.0/10.0*curve->n);
    } else {
        curve->n_C_0 = ceil(frac_n_C_0*curve->n);
    }

    if(frac_n_C_1 == 0) {
        curve->n_C_1 = ceil(1.0/10.0*curve->n);
    } else {
        curve->n_C_1 = ceil(frac_n_C_1*curve->n);
    }

    if(frac_n_S_0 == 0) {
        curve->n_S_0 = ceil(2.0/3.0*curve->n_C_0);
    } else {
        curve->n_S_0 = ceil(frac_n_S_0*curve->n_C_0);
    }

    if(frac_n_S_1 == 0) {
        curve->n_S_1 = ceil(2.0/3.0*curve->n_C_1);
    } else {
        curve->n_S_1 = ceil(frac_n_S_1*curve->n_C_1);
    }

    curve->h_tan = 1.0/(curve->n-1);
    curve->h_norm = h_norm;

    if(next_curve == NULL) {
        curve->next_curve = curve;
    } else {
        curve->next_curve = next_curve;
    }
}

void curve_seq_add_curve(curve_seq_t *curve_seq, curve_t *curve, scalar_func_t l_1, scalar_func_t l_2, scalar_func_t l_1_prime, scalar_func_t l_2_prime, scalar_func_t l_1_dprime, scalar_func_t l_2_dprime, MKL_INT n, double frac_n_C_0, double frac_n_C_1, double frac_n_S_0, double frac_n_S_1, double h_norm) {
    curve_init(curve, l_1, l_2, l_1_prime, l_2_prime, l_1_dprime, l_2_dprime, n, frac_n_C_0, frac_n_C_1, frac_n_S_0, frac_n_S_1, h_norm, curve_seq->first_curve);
    if(curve_seq->first_curve == NULL) {
        curve_seq->first_curve = curve;
        curve_seq->last_curve = curve;
    }

    curve_seq->last_curve->next_curve = curve;
    curve_seq->last_curve = curve;            
    curve_seq->n_curves += 1;
}

MKL_INT curve_S_patch_num_el(curve_t *curve, double d) {
    MKL_INT s_patch_n_xi = curve->n - (curve->n_C_1-curve->n_S_1) - (curve->n_C_0-curve->n_S_0);
    return s_patch_n_xi * d;
}

void curve_construct_S_patch(curve_t *curve, s_patch_t *s_patch, rd_mat_t *s_patch_f_XY, scalar_func_2D_t f, MKL_INT d, double eps_xi_eta, double eps_xy) {
    double xi_diff = 1-(curve->n_C_1-curve->n_S_1)*curve->h_tan - (curve->n_C_0-curve->n_S_0)*curve->h_tan;
    double xi_0 = (curve->n_C_0 - curve->n_S_0) * curve->h_tan;

    curve->M_p_S_general_param = (M_p_S_general_param_t) {xi_diff, xi_0, curve->l_1, curve->l_2, curve->l_1_prime, curve->l_2_prime};
    M_p_general_t M_p_general = {(M_p_general_handle_t) M_p_S_general_handle, (void*) &(curve->M_p_S_general_param)};

    curve->J_S_general_param = (J_S_general_param_t) {xi_diff, xi_0, curve->l_1, curve->l_2, curve->l_1_prime, curve->l_2_prime, curve->l_1_dprime, curve->l_2_dprime};
    J_general_t J_general = {(J_general_handle_t) J_S_general_handle, (void*) &(curve->J_S_general_param)};

    MKL_INT s_patch_n_xi = curve->n - (curve->n_C_1-curve->n_S_1) - (curve->n_C_0-curve->n_S_0);
    s_patch_init(s_patch, M_p_general, J_general, curve->h_norm, eps_xi_eta, eps_xy, s_patch_n_xi, d, s_patch_f_XY);

    q_patch_evaluate_f(&(s_patch->q_patch), (scalar_func_2D_t) f);
}

