#ifndef __R_CARTESIAN_MESH_LIB__
#define __R_CARTESIAN_MESH_LIB__

#include "num_linalg_lib.h"
#include <mkl.h>
#include <string.h>

void inpolygon_mesh(rd_mat_t R_X, rd_mat_t R_Y, rd_mat_t boundary_X, rd_mat_t boundary_Y, ri_mat_t *in_msk);

typedef struct r_cartesian_mesh_obj {
    double x_start;
    double x_end;
    double y_start;
    double y_end;
    double h;
    MKL_INT n_x;
    MKL_INT n_y;

    ri_mat_t *in_interior;
    ri_mat_t *interior_idxs;
    rd_mat_t *f_R;
    
    rd_mat_t *fc_coeffs;
} r_cartesian_mesh_obj_t;

void r_cartesian_mesh_init(r_cartesian_mesh_obj_t *r_cartesian_mesh_obj, double x_start, double x_end, double y_start, double y_end, double h, rd_mat_t boundary_X, rd_mat_t boundary_Y, ri_mat_t *in_interior, ri_mat_t *interior_idxs, rd_mat_t *f_XY);


#endif