#include <stdlib.h>
#include <math.h>

#include "r_cartesian_mesh_lib.h"
#include "num_linalg_lib.h"

void r_cartesian_mesh_init(r_cartesian_mesh_obj_t *r_cartesian_mesh_obj, double x_start, double x_end, double y_start, double y_end, double h, rd_mat_t boundary_X, rd_mat_t boundary_Y, ri_mat_t *in_interior, ri_mat_t *interior_idxs, rd_mat_t *f_XY) {
    r_cartesian_mesh_obj->x_start = x_start;
    r_cartesian_mesh_obj->y_start = y_start;

    r_cartesian_mesh_obj->x_end = ceil((x_end-x_start)/h)*h + x_start;
    r_cartesian_mesh_obj->y_end = ceil((y_end-y_start)/h)*h + y_start;

    r_cartesian_mesh_obj->h = h;
    
    r_cartesian_mesh_obj->n_x = round((r_cartesian_mesh_obj->x_end-r_cartesian_mesh_obj->x_start)/h) + 1;
    r_cartesian_mesh_obj->n_x = round((r_cartesian_mesh_obj->y_end-r_cartesian_mesh_obj->y_start)/h) + 1;
}