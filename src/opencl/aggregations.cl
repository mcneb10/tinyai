/* OpenCL aggregation kernels */

// All aggregations take 2 parameters, a pointer to a double buffer

// the gid is the current index being processed in the array
kernel void agg_sum(global double* input, global double output) {
    uint gid = get_global_id(0);
    atomic_add(&output, input[gid]);
}

kernel void agg_product(global double* input, global double output) {
    uint gid = get_global_id(0);
    if (gid == 0) {
        output = input[0];
    } else {
        output *= input[gid];
    }
}

kernel void agg_max(global double* input, global double output) {
    uint gid = get_global_id(0);
    if(gid == 0) {
        output = input[0];
    } else {
        atomic_max(&output, input[gid]);
    }
}

kernel void agg_min(global double* input, global double output) {
    uint gid = get_global_id(0);
    if(gid == 0) {
        output = input[0];
    } else {
        atomic_min(&output, input[gid]);
    }
}

