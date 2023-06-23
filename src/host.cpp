/******************************************
*MIT License
*
# *Copyright (c) [2020] [Beatrice Branchini, Luisa Cicolini, Giulia Gerometta, Marco Santambrogio]
*
*Permission is hereby granted, free of charge, to any person obtaining a copy
*of this software and associated documentation files (the "Software"), to deal
*in the Software without restriction, including without limitation the rights
*to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
*copies of the Software, and to permit persons to whom the Software is
*furnished to do so, subject to the following conditions:
*
*The above copyright notice and this permission notice shall be included in all
*copies or substantial portions of the Software.
*
*THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
*IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
*FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
*AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
*LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
*OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*SOFTWARE.
******************************************/

#include <fstream>
#include "xcl2.hpp"
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <climits>
#include <time.h>
#include <chrono>
#include <cstdint>

#define NUM_NODES 100
#define NUM_NEIGHBORS 2
#define NUM_REGIONS 10
#define MAX 9223372036854775807
#define PORT_WIDTH 32

#define NUM_KERNEL 1

#define NOW std::chrono::high_resolution_clock::now();

#define MAX_HBM_BANKCOUNT 32
#define BANK_NAME(n) n | XCL_MEM_TOPOLOGY
const int bank[MAX_HBM_BANKCOUNT] = {
    BANK_NAME(0),  BANK_NAME(1),  BANK_NAME(2),  BANK_NAME(3),  BANK_NAME(4),
    BANK_NAME(5),  BANK_NAME(6),  BANK_NAME(7),  BANK_NAME(8),  BANK_NAME(9),
    BANK_NAME(10), BANK_NAME(11), BANK_NAME(12), BANK_NAME(13), BANK_NAME(14),
    BANK_NAME(15), BANK_NAME(16), BANK_NAME(17), BANK_NAME(18), BANK_NAME(19),
    BANK_NAME(20), BANK_NAME(21), BANK_NAME(22), BANK_NAME(23), BANK_NAME(24),
    BANK_NAME(25), BANK_NAME(26), BANK_NAME(27), BANK_NAME(28), BANK_NAME(29),
    BANK_NAME(30), BANK_NAME(31)};


std::pair<size_t, uint64_t > find_next_event_at_node_occupied_by_growing_top_region(
	uint32_t detector_node,
	uint64_t rad1,
	uint32_t * num_neighbors,
	uint32_t neighbors[][NUM_NEIGHBORS],
	uint32_t neighbor_weights[][NUM_NEIGHBORS],
	uint32_t * region_that_arrived_top,
	uint32_t * wrapped_radius_cached,
	uint64_t * radius)

	{
	uint64_t best_time = MAX;
	uint32_t best_neighbor = MAX;
	uint32_t start = 0;
    if (!(num_neighbors[detector_node]==0) && neighbors[detector_node][0] == -1) {
        // Growing towards boundary
    	uint32_t weight = neighbor_weights[detector_node][0];
    	uint64_t collision_time = weight - ((rad1 >> 2)<<2);
        if (collision_time < best_time) {
            best_time = collision_time;
            best_neighbor = 0;
        }
        start++;
    }

    // Handle non-boundary neighbors.
    for (uint32_t i = start; i < num_neighbors[detector_node]; i++) {
        
    	uint32_t weight = neighbor_weights[detector_node][i];

    	uint32_t neighbor = neighbors[detector_node][i];

        if (region_that_arrived_top[detector_node] == region_that_arrived_top[neighbor]) {
            continue;
        }
        uint64_t rad2;

		if (region_that_arrived_top[neighbor] == -1) {
				rad2 = 0;
		} else {
			rad2 = radius[region_that_arrived_top[neighbor]] +(wrapped_radius_cached[neighbor] << 2);
		}

        if (rad2 & 2) {
            continue;
        }

        uint64_t collision_time = weight - ((rad1 >> 2) << 2) - ((rad2 >> 2) << 2);
        if (rad2 & 1) {
            collision_time >>= 1;
        }
        if (collision_time < best_time) {
            best_time = collision_time;
            best_neighbor = i;
        }
    }
    return {best_neighbor, best_time};
}

std::pair<size_t, uint64_t > find_next_event_at_node_not_occupied_by_growing_top_region(
	uint32_t detector_node,
	uint64_t rad1,
	uint32_t * num_neighbors,
	uint32_t neighbors[][NUM_NEIGHBORS],
	uint32_t neighbor_weights[][NUM_NEIGHBORS],
	uint32_t * region_that_arrived_top,
	uint32_t * wrapped_radius_cached,
	uint64_t * radius)

{
	uint64_t best_time = MAX;
	uint32_t best_neighbor = MAX;

	uint32_t start = 0;
    if (!(num_neighbors[detector_node]==0) && neighbors[detector_node][0] == -1)
        start++;

    // Handle non-boundary neighbors.
    for (uint32_t i = start; i < num_neighbors[detector_node]; i++) {
        
    	uint32_t weight = neighbor_weights[detector_node][i];

    	uint32_t neighbor = neighbors[detector_node][i];

    	uint64_t rad2;

		if (region_that_arrived_top[neighbor] == -1) {
				rad2 = 0;
		} else {
			rad2 = radius[region_that_arrived_top[neighbor]] +(wrapped_radius_cached[neighbor] << 2);
		}

        if (rad2 & 1) {
            auto collision_time = weight - ((rad1 >> 2) << 2) - ((rad2 >> 2) << 2);
            if (collision_time < best_time) {
                best_time = collision_time;
                best_neighbor = i;
            }
        }
    }
    return {best_neighbor, best_time};
}

std::pair<size_t, uint64_t > find_next_event_at_node_returning_neighbor_index_and_time(
    uint32_t detector_node,
	uint32_t * num_neighbors,
	uint32_t neighbors[][NUM_NEIGHBORS],
	uint32_t neighbor_weights[][NUM_NEIGHBORS],
	uint32_t * region_that_arrived_top,
	uint32_t * wrapped_radius_cached,
	uint64_t * radius)
{
	uint64_t rad1;

	if (region_that_arrived_top[detector_node] == -1) {
			rad1 = 0;
	} else {
		rad1 = radius[region_that_arrived_top[detector_node]] + (wrapped_radius_cached[detector_node] << 2);
	}


    if (rad1 & 1) {
        return find_next_event_at_node_occupied_by_growing_top_region(detector_node, rad1, num_neighbors, neighbors, neighbor_weights, region_that_arrived_top, wrapped_radius_cached, radius);
    } else {
        return find_next_event_at_node_not_occupied_by_growing_top_region(detector_node, rad1, num_neighbors, neighbors, neighbor_weights, region_that_arrived_top, wrapped_radius_cached, radius);
    }
}


int main(int argc, char *argv[]){
    
    std::string binaryFile = "querk.xclbin";
	std::string readsPath;
	
    cl::Context context;
 
    cl::CommandQueue commands;
    std::vector<uint32_t, aligned_allocator<uint32_t>> num_neighbors(NUM_NODES);
    num_neighbors = {1,1};
    std::vector<uint64_t, aligned_allocator<uint64_t>> radius(NUM_REGIONS);
    radius = {1,1};
    std::vector<uint32_t, aligned_allocator<uint32_t>> region_that_arrived_top(NUM_NODES);
    region_that_arrived_top = {0,1};
    std::vector<uint32_t, aligned_allocator<uint32_t>> wrapped_radius_cached(NUM_NODES);
    wrapped_radius_cached = {1,1};
    std::vector<uint32_t, aligned_allocator<uint32_t>> neighbors(NUM_NODES*NUM_NEIGHBORS);
    neighbors = {1,0};
    std::vector<uint32_t, aligned_allocator<uint32_t>> neighbor_weights(NUM_NODES*NUM_NEIGHBORS);
    neighbor_weights = {1,1};
    std::vector<uint64_t, aligned_allocator<uint64_t>> neighbor_observables(NUM_NODES*NUM_NEIGHBORS);
    neighbor_observables = {0,0};
    std::vector<uint32_t, aligned_allocator<uint32_t>> out_neighbor(1);
    out_neighbor = {-1};
    std::vector<uint64_t, aligned_allocator<uint64_t>> out_time(1);
    out_time = {-1};

    uint32_t detector_node = 0;
    uint32_t num_nodes = 2;
    uint32_t num_regions = 2;
	
    if (argc == 3) { //Input provided by file 

        binaryFile = argv[1];
        readsPath = argv[2];
    } else {
         binaryFile = argv[1];
    }   
    //     std::ifstream st(readsPath);
	//     std::string temp1;
	//     std::string temp2;
	    
    //     int idx = 0;
    //     while (st >> temp1 >> temp2){ //Reading input from file 
	// 	    std::string temp;
	// 	    temp = temp1.substr(1,std::string::npos);
	// 	    pattern_test[idx] = temp;
	// 	    temp = temp2.substr(1,std::string::npos);
    //         text_test[idx] = temp;
    //         idx++;
	//     }

    //     num = idx;

    //     for (int i=0; i<NUM; i++){ //Storing lengths
    //         patternLength[i] = pattern_test[i].size();
    //         textLength[i] = text_test[i].size();
    //         if (i>=num){
    //             patternLength[i] = 0;
    //             textLength[i] = 0;
    //         }
    //     }

    //     int iter_pattern = 0;
    //     int iter_text = 0;
    //     for (int i=0; i<num; i++){ //Reorganizing inputs for the core
            
    //         for (int j=0; j<pattern_test[i].size(); j++){
    //             pattern[iter_pattern] = pattern_test[i][j];
    //             iter_pattern++;
    //         }
            
    //         for (int k=0; k<text_test[i].size(); k++){
    //             text[iter_text] = text_test[i][k];
    //             iter_text++;
    //         }
    //     }
    
    // } else { //Random generation of inputs
 
    //     binaryFile = argv[1];

    //     std::cout<<"Seed: "<<seed<<std::endl;
	//     srand(seed);

	//     char alphabet[4] = {'A', 'C', 'G', 'T'};

	//     int iter = 0; 
	//     for (int i = 0; i < num; i++) {
	// 	    for (int j = 0; j < TEXT_SIZE; j++) {
	// 		    text_test[i].push_back(alphabet[rand()%4]);
	// 		    text[iter]=text_test[i][j];
	// 		    iter++;
	// 	    }
	//     }

	//     iter = 0;
	//     for (int i = 0; i < num; i++) {
	// 	    for (int j = 0; j < PATTERN_SIZE; j++) {
	// 		    pattern_test[i].push_back(alphabet[rand()%4]);
	// 		    pattern[iter] = pattern_test[i][j];
	// 		    iter++;
	// 	    }
	//     }

	//     for (int i = 0; i < num; i++) {
	// 	    patternLength[i] = PATTERN_SIZE;
	// 	    textLength[i] = TEXT_SIZE;
	//     }
    // }

    // // Set penalties
	// affine_penalties_t affine_penalties = {
	// 	.match = 0,
	// 	.mismatch = 3,
	// 	.gap_opening = 5,
	// 	.gap_extension = 1,
	// };

	// if ((affine_penalties.mismatch-affine_penalties.gap_extension) != (affine_penalties.gap_opening-affine_penalties.mismatch)) {
	// 	std::cout<<"Error: unsupported penalties!"<<std::endl;
	// 	std::cout<<"Must satisfy: (mismatch - gap extension) == (gap opening - mismatch)"<<std::endl;
	// 	std::cout<<"Test failed"<<std::endl;
	// 	exit(1);
	// }

	// printf("PENALTIES INITIALIZED: %d, %d, %d, %d \n", affine_penalties.match, affine_penalties.mismatch, affine_penalties.gap_opening, affine_penalties.gap_extension);

    std::string krnl_name = "querk";
    cl::Kernel krnl;

    // The get_xil_devices will return vector of Xilinx Devices
    auto devices = xcl::get_xil_devices();

    // read_binary_file() command will find the OpenCL binary file created using the
    // V++ compiler load into OpenCL Binary and return pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binaryFile);

    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    int valid_device = 0;

    cl_int err;

    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
            // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));
        OCL_CHECK(err, commands = cl::CommandQueue(context, device,
                            CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE, &err));

        std::cout << "Trying to program device[" << i 
                  << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
     
        cl::Program program(context, {device}, bins, NULL, &err);
        
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i
                        << "] with xclbin file!\n";                      
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            
            // Creating Kernel object using Compute unit names
            
            std::string cu_id = std::to_string(1);
            std::string krnl_name_full = krnl_name + ":{" + "querk_" + cu_id + "}";

            printf("Creating a kernel [%s] for CU(%d)\n", krnl_name_full.c_str(),  1);

            //Here Kernel object is created by specifying kernel name along with compute unit.
            //For such case, this kernel object can only access the specific Compute unit
            OCL_CHECK(err, krnl = cl::Kernel(program, krnl_name_full.c_str(), &err));
            

            valid_device++;
            break; // we break because we found a valid device
        }
        std::cout<<"dwvgae"<<std::endl;
    }

	std::cout<<"Kernel created"<<std::endl;
    
    if (valid_device == 0) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    // Create device buffers
    cl_mem_ext_ptr_t num_neighbors_ext;
    cl_mem_ext_ptr_t radius_ext;
    cl_mem_ext_ptr_t region_that_arrived_top_ext;
    cl_mem_ext_ptr_t wrapped_radius_cached_ext;
    cl_mem_ext_ptr_t neighbors_ext;
    cl_mem_ext_ptr_t neighbor_weights_ext;
    cl_mem_ext_ptr_t neighbor_observables_ext;
    cl_mem_ext_ptr_t out_neighbor_ext;
    cl_mem_ext_ptr_t out_time_ext;
    cl::Buffer num_neighbors_buffer;
    cl::Buffer radius_buffer;
    cl::Buffer region_that_arrived_top_buffer;
    cl::Buffer wrapped_radius_cached_buffer;
    cl::Buffer neighbors_buffer;
    cl::Buffer neighbor_weights_buffer;
    cl::Buffer neighbor_observables_buffer;
    cl::Buffer out_neighbor_buffer;
    cl::Buffer out_time_buffer;


    num_neighbors_ext.obj = num_neighbors.data();
    num_neighbors_ext.param = 0;
    num_neighbors_ext.flags = bank[0];

    radius_ext.obj = radius.data();
    radius_ext.param = 0;
    radius_ext.flags = bank[1];
    
    region_that_arrived_top_ext.obj = region_that_arrived_top.data();
    region_that_arrived_top_ext.param = 0;
    region_that_arrived_top_ext.flags = bank[2];

    wrapped_radius_cached_ext.obj = wrapped_radius_cached.data();
    wrapped_radius_cached_ext.param = 0;
    wrapped_radius_cached_ext.flags = bank[3];

    neighbors_ext.obj = neighbors.data();
    neighbors_ext.param = 0;
    neighbors_ext.flags = bank[4];

    neighbor_weights_ext.obj = neighbor_weights.data();
    neighbor_weights_ext.param = 0;
    neighbor_weights_ext.flags = bank[5];
    
    neighbor_observables_ext.obj = neighbor_observables.data();
    neighbor_observables_ext.param = 0;
    neighbor_observables_ext.flags = bank[6];

    out_neighbor_ext.obj = out_neighbor.data();
    out_neighbor_ext.param = 0;
    out_neighbor_ext.flags = bank[7];

    out_time_ext.obj = out_time.data();
    out_time_ext.param = 0;
    out_time_ext.flags = bank[8];

    

    OCL_CHECK(err, num_neighbors_buffer = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                CL_MEM_USE_HOST_PTR, sizeof(int)*NUM_NODES, &num_neighbors_ext, &err));
    OCL_CHECK(err, radius_buffer = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                CL_MEM_USE_HOST_PTR, sizeof(long int)*NUM_REGIONS, &radius_ext, &err));
    OCL_CHECK(err, region_that_arrived_top_buffer = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                CL_MEM_USE_HOST_PTR, sizeof(int)*NUM_NODES, &region_that_arrived_top_ext, &err));
    OCL_CHECK(err, wrapped_radius_cached_buffer = cl::Buffer(context, CL_MEM_READ_ONLY | CL_MEM_EXT_PTR_XILINX |
                                CL_MEM_USE_HOST_PTR, sizeof(int)*NUM_NODES, &wrapped_radius_cached_ext, &err));
    OCL_CHECK(err, neighbors_buffer = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
                                CL_MEM_USE_HOST_PTR, sizeof(int)*NUM_NODES*NUM_NEIGHBORS, &neighbors_ext, &err));
	
    OCL_CHECK(err, neighbor_weights_buffer = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
                                CL_MEM_USE_HOST_PTR, sizeof(int)*NUM_NODES*NUM_NEIGHBORS, &neighbor_weights_ext, &err));

    OCL_CHECK(err, neighbor_observables_buffer = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
                                CL_MEM_USE_HOST_PTR, sizeof(long int)*NUM_NODES*NUM_NEIGHBORS, &neighbor_observables_ext, &err));

    OCL_CHECK(err, out_neighbor_buffer = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
                                CL_MEM_USE_HOST_PTR, sizeof(int), &out_neighbor_ext, &err));
    
    OCL_CHECK(err, out_time_buffer = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
                              CL_MEM_USE_HOST_PTR, sizeof(long int), &out_time_ext, &err));





	commands.finish();

    // Write our data set into device buffers  
     
    err = commands.enqueueMigrateMemObjects({num_neighbors_buffer, radius_buffer, region_that_arrived_top_buffer, wrapped_radius_cached_buffer, neighbors_buffer, neighbor_weights_buffer, neighbor_observables_buffer, out_neighbor_buffer, out_time_buffer}, 0);

    if (err != CL_SUCCESS) {
            printf("Error: Failed to write to device memory!\n");
            printf("Test failed\n");
            exit(1);
    }

	commands.finish();
    
    // Set the arguments to our compute kernel
    OCL_CHECK(err, err = krnl.setArg(0, detector_node));
    OCL_CHECK(err, err = krnl.setArg(1, num_nodes));
    OCL_CHECK(err, err = krnl.setArg(2, num_regions));
    OCL_CHECK(err, err = krnl.setArg(3, num_neighbors_buffer));
    OCL_CHECK(err, err = krnl.setArg(4, radius_buffer));
    OCL_CHECK(err, err = krnl.setArg(5, region_that_arrived_top_buffer));
    OCL_CHECK(err, err = krnl.setArg(6, wrapped_radius_cached_buffer));
    OCL_CHECK(err, err = krnl.setArg(7, neighbors_buffer));
    OCL_CHECK(err, err = krnl.setArg(8, neighbor_weights_buffer));
    OCL_CHECK(err, err = krnl.setArg(9, neighbor_observables_buffer));
    OCL_CHECK(err, err = krnl.setArg(10, out_neighbor_buffer));
    OCL_CHECK(err, err = krnl.setArg(11, out_time_buffer));

    if (err != CL_SUCCESS) {
        printf("Error: Failed to set kernel arguments! %d\n", err);
        printf("Test failed\n");
        exit(1);
    }
    

	commands.finish();

    std::chrono::high_resolution_clock::time_point start = NOW;
    
    // Execute the kernel over the entire range of our 1d input data set
    // using the maximum number of work group items for this device
    
    err = commands.enqueueTask(krnl);


    if (err) {
        printf("Error: Failed to execute kernel! %d\n", err);
        printf("Test failed\n");
        exit(1);
    }

    commands.finish();
    std::chrono::high_resolution_clock::time_point end = NOW;
	std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);
    
    // Read back the results from the device to verify the output
  
    err = commands.enqueueMigrateMemObjects({out_neighbor_buffer, out_time_buffer}, CL_MIGRATE_MEM_OBJECT_HOST);  
    commands.finish();
    printf("Hardware results: %d %ld\n", (int) out_neighbor[0], (long int) out_time[0] );

    if (err != CL_SUCCESS) {
        printf("Error: Failed to read output array! %d\n", err);
        printf("Test failed\n");
        exit(1);
    }

	//printf("HW time: %lf\n", time);

	//Checking the results 


    start = NOW;

    auto golden = find_next_event_at_node_returning_neighbor_index_and_time(detector_node, num_neighbors.data(), (uint32_t (*) [NUM_NEIGHBORS]) neighbors.data(), (uint32_t (*) [NUM_NEIGHBORS]) neighbor_weights.data(), region_that_arrived_top.data(), wrapped_radius_cached.data(), radius.data());
	


    // for(unsigned i = 0 ; i < num; i++){
	// 	// Allocate MM
	// 	mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
	// 	// Init Affine-WFA
	// 	affine_wavefronts_t* affine_wavefronts = affine_wavefronts_new_complete(
	// 		pattern_test[i].size(), text_test[i].length(), &affine_penalties, NULL, mm_allocator);
	// 	// Align
	// 	affine_wavefronts_align(affine_wavefronts, pattern_test[i].c_str(), pattern_test[i].length(),
	// 		text_test[i].c_str(), text_test[i].length(), &sw_scores[i]);

	// 	affine_wavefronts_delete(affine_wavefronts);
	// 	mm_allocator_delete(mm_allocator);

	// }

	end = NOW;
	time = std::chrono::duration_cast<std::chrono::duration<double>>(end-start);

    std::cout << "Golden results: " <<  golden.first << " " << golden.second << std::endl;

	//printf("SW time: %lf\n", time);


	// for (int i=0; i<num; i++){
	// 	if (scores[i]!=sw_scores[i]){
    //         printf("HW: %d, SW: %d\n", scores[i], sw_scores[i]);
    //         test_score=false;
    //     }
	// }

	// if (test_score) 
	// 	std::cout<<"All scores correct"<<std::endl;
	// else 
	// 	std::cout<<"Test failed"<<std::endl;

}
