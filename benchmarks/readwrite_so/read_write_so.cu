#include <cuda.h>
#include <nvm_ctrl.h>
#include <nvm_types.h>
#include <nvm_queue.h>
#include <nvm_util.h>
#include <nvm_admin.h>
#include <nvm_error.h>
#include <nvm_cmd.h>
#include <string>
#include <stdexcept>
#include <vector>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <map>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <ctrl.h>
#include <buffer.h>
#include <event.h>
#include <queue.h>
#include <nvm_parallel_queue.h>
#include <nvm_io.h>
#include <page_cache.h>
#include <util.h>
#include <iostream>
#include <fstream>
#include <byteswap.h>
#ifdef __DIS_CLUSTER__
#include <sisci_api.h>
#endif
#include "read_write_so.h"

#define READ 0
#define WRITE 1
#define MIXED 2
#define VERIFY 3

using error = std::runtime_error;
using std::string;

 //const char* const ctrls_paths[] = {"/dev/libnvm0", "/dev/libnvm1", "/dev/libnvm2", "/dev/libnvm3", "/dev/libnvm4", "/dev/libnvm5", "/dev/libnvm6"};
const char* const ctrls_paths[] = {"/dev/libnvm0"};

/*
__device__ void read_data(page_cache_t* pc, QueuePair* qp, const uint64_t starting_lba, const uint64_t n_blocks, const unsigned long long pc_entry) {
    //uint64_t starting_lba = starting_byte >> qp->block_size_log;
    //uint64_t rem_bytes = starting_byte & qp->block_size_minus_1;
    //uint64_t end_lba = CEIL((starting_byte+num_bytes), qp->block_size);

    //uint16_t n_blocks = CEIL(num_bytes, qp->block_size, qp->block_size_log);
 
    nvm_cmd_t cmd;
    uint16_t cid = get_cid(&(qp->sq));
    //printf("cid: %u\n", (unsigned int) cid);
 
    nvm_cmd_header(&cmd, cid, NVM_IO_READ, qp->nvmNamespace);
    uint64_t prp1 = pc->prp1[pc_entry];
    uint64_t prp2 = 0;
    if (pc->prps)
        prp2 = pc->prp2[pc_entry];
    //printf("tid: %llu\tstart_lba: %llu\tn_blocks: %llu\tprp1: %p\n", (unsigned long long) threadIdx.x, (unsigned long long) starting_lba, (unsigned long long) n_blocks, (void*) prp1);
    nvm_cmd_data_ptr(&cmd, prp1, prp2);
    nvm_cmd_rw_blks(&cmd, starting_lba, n_blocks);
    uint16_t sq_pos = sq_enqueue(&qp->sq, &cmd);

    uint32_t cq_pos = cq_poll(&qp->cq, cid);
    sq_dequeue(&qp->sq, sq_pos);
    cq_dequeue(&qp->cq, cq_pos);
 
    put_cid(&qp->sq, cid);
 
}
*/

__global__
void sequential_access_kernel(Controller** ctrls, page_cache_d_t* pc,  uint32_t req_size, uint32_t n_reqs, //unsigned long long* req_count,
                                uint32_t num_ctrls, uint64_t reqs_per_thread, uint32_t access_type, uint64_t s_offset, uint64_t o_offset){
    //printf("in threads\n");
    uint64_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    // uint32_t bid = blockIdx.x;
    // uint32_t smid = get_smid();

    uint32_t ctrl = (tid/32) % (num_ctrls);
    uint32_t queue = (tid/32) % (ctrls[ctrl]->n_qps);
    uint64_t itr=0; 

//    printf("Num pages: %llu, s_offset: %llu n_reqs: %llu\t req_size: %llu\n", (unsigned long long int) pc->n_pages, (unsigned long long int) s_offset, (unsigned long long int) n_reqs, (unsigned long long) req_size); 
    for (;tid < pc->n_pages; tid = tid+n_reqs){
            uint64_t start_block = (o_offset+s_offset + tid*req_size) >> ctrls[ctrl]->d_qps[queue].block_size_log ;
            uint64_t pc_idx = (tid);
            //uint64_t start_block = (tid*req_size) >> ctrls[ctrl]->d_qps[queue].block_size_log;
            //start_block = tid;
            uint64_t n_blocks = req_size >> ctrls[ctrl]->d_qps[queue].block_size_log; /// ctrls[ctrl].ns.lba_data_size;;
//            printf("itr:%llu\ttid: %llu\tstart_block: %llu\tn_blocks: %llu\tpc_idx: %llu\n", (unsigned long long)itr, (unsigned long long) tid, (unsigned long long) start_block, (unsigned long long) n_blocks, (unsigned long long) pc_idx);
            itr = itr+1; 
            // uint8_t opcode;
            // for (size_t i = 0; i < reqs_per_thread; i++) {
                if (access_type == READ) {
                    read_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, pc_idx);
                    //if(tid ==( pc->n_pages - 1)){
                    //        printf("I am here\n");
                    //        hexdump(pc->base_addr+tid*req_size, 4096); 
                    //}
                }
                else {
                    write_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, pc_idx);
                }

            
            // }
            //read_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, tid);
            //read_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, tid);
            //read_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, tid);
            //__syncthreads();
            //read_data(pc, (ctrls[ctrl].d_qps)+(queue),start_block*2, n_blocks, tid);
            //printf("tid: %llu finished\n", (unsigned long long) tid);
    }
}

/*__global__
void random_access_kernel(Controller** ctrls, page_cache_t* pc,  uint32_t req_size, uint32_t n_reqs, unsigned long long* req_count, uint32_t num_ctrls, uint64_t* assignment, uint64_t reqs_per_thread, uint32_t access_type, uint8_t* access_type_assignment) {
    //printf("in threads\n");
    uint64_t tid = blockIdx.x * blockDim.x + threadIdx.x;
    uint32_t bid = blockIdx.x;
    uint32_t smid = get_smid();

    uint32_t ctrl = (tid/32) % (num_ctrls);
    uint32_t queue = (tid/32) % (ctrls[ctrl]->n_qps);


    if (tid < n_reqs) {
        uint64_t start_block = (assignment[tid]*req_size) >> ctrls[ctrl]->d_qps[queue].block_size_log;
        //uint64_t start_block = (tid*req_size) >> ctrls[ctrl]->d_qps[queue].block_size_log;
        //start_block = tid;
        uint64_t n_blocks = req_size >> ctrls[ctrl]->d_qps[queue].block_size_log; /// ctrls[ctrl].ns.lba_data_size;;
        //printf("tid: %llu\tstart_block: %llu\tn_blocks: %llu\n", (unsigned long long) tid, (unsigned long long) start_block, (unsigned long long) n_blocks);

        uint8_t opcode;
        for (size_t i = 0; i < reqs_per_thread; i++) {
            if (access_type == MIXED) {
                opcode = access_type_assignment[tid];
                access_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, tid, opcode);
            }
            else if (access_type == READ) {
                read_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, tid);

            }
            else {
                write_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, tid);
            }
        }
        //read_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, tid);
        //read_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, tid);
        //read_data(pc, (ctrls[ctrl]->d_qps)+(queue),start_block, n_blocks, tid);
        //__syncthreads();
        //read_data(pc, (ctrls[ctrl].d_qps)+(queue),start_block*2, n_blocks, tid);
        //printf("tid: %llu finished\n", (unsigned long long) tid);

    }

}
*/

__global__ 
void verify_kernel(uint64_t* orig_h, uint64_t* nvme_h, uint64_t n_elems,uint32_t n_reqs){
        uint64_t tid = blockIdx.x*blockDim.x + threadIdx.x; 

        for (;tid < n_elems; tid = tid+n_reqs){
           uint64_t orig_val = orig_h[tid]; 
           uint64_t nvme_val = nvme_h[tid]; 
           if(orig_val != nvme_val)
              printf("MISMATCH: at %llu\torig_val:%llu\tnvme_val:%llu\tn_reqs:%lu\tn_elms:%llu\n",tid, (unsigned long long)orig_val, (unsigned long long)nvme_h, n_reqs, n_elems);
        }
        __syncthreads();//really not needed. 
}


uint32_t nvme_dev_setting(uint32_t cudaDevice)
{
    cudaDeviceProp properties;
    if (cudaGetDeviceProperties(&properties, cudaDevice) != cudaSuccess)
    {
        fprintf(stderr, "Failed to get CUDA device properties\n");
        return 1;
    }

    cuda_err_chk(cudaSetDevice(cudaDevice));
    std::vector<Controller*> ctrls(1);
    for (size_t i = 0 ; i < 1; i++)
        ctrls[i] = new Controller(ctrls_paths[i], 1, cudaDevice, 1024, 1);
}

  uint32_t nvme_dev_read()
  {
    fprintf(stdout, "read\n");

    return 1;
  }

  uint32_t nvme_dev_write()
  {
    fprintf(stdout, "write\n");

    return 2;
  }