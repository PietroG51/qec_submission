#ifndef KERNEL_SIMPLE_H
#define KERNEL_SIMPLE_H

#define NUM_NODES 100
#define NUM_REGIONS 10
#define NUM_NEIGHBORS 2
#define MAX 9223372036854775807

extern "C" void querk(ap_uint<32> detector_node, ap_uint<32> num_nodes, ap_uint<32> num_regions, ap_uint<32> * num_neighbors, ap_uint<64> * radius, ap_uint<32> * region_that_arrived_top, ap_uint<32> * wrapped_radius_cached, ap_uint<32> neighbors[][NUM_NEIGHBORS], ap_uint<32> neighbor_weights[][NUM_NEIGHBORS], ap_uint<64> neighbor_observables[][NUM_NEIGHBORS], ap_uint<32> * out_neighbor, ap_uint<64> * out_time);

#endif
