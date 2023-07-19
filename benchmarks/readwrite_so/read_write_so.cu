#include "read_write_so.h"

#define READ 0
#define WRITE 1
#define MIXED 2
#define VERIFY 3

using error = std::runtime_error;
using std::string;

 //const char* const ctrls_paths[] = {"/dev/libnvm0", "/dev/libnvm1", "/dev/libnvm2", "/dev/libnvm3", "/dev/libnvm4", "/dev/libnvm5", "/dev/libnvm6"};

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
void random_access_kernel(Controller** ctrls, page_cache_d_t* pc,  uint32_t req_size, uint32_t n_reqs, uint32_t num_ctrls, uint64_t* assignment, uint64_t reqs_per_thread, uint32_t access_type) {
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

        for (size_t i = 0; i < reqs_per_thread; i++) {
            if (access_type == READ) {
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


int dev_set(uint32_t cudaDevice, void* src_in)
{

    int fd_in;
    char* input_f = (char*) src_in;

    if((fd_in = open(input_f, O_RDONLY)) == -1)
    {
        fprintf(stderr, "Input file cannot be opened!\n");
        return 1;
    }
    fstat(fd_in, &sb_in);

    printf("sb_in.st_size: %ld\n", sb_in.st_size);

    map_in = mmap(NULL, sb_in.st_size, PROT_READ, MAP_SHARED, fd_in, 0);
    if(map_in == (void*)-1)
    {
        fprintf(stderr, "Input file map failed %d\n", map_in);
        return 1;
    }

    n_tsteps = ceil((float)(sb_in.st_size-0)/(float)total_cache_size);
    n_telem = ((sb_in.st_size-0)/sizeof(int64_t));

    // for(int id=0; id<30; id++)
    // {
    //     printf("id[%d] = %.3f\n", id, ((float*)map_in)[id]);
    // }

    cudaDeviceProp properties;
    if (cudaGetDeviceProperties(&properties, cudaDevice) != cudaSuccess)
    {
        fprintf(stderr, "Failed to get CUDA device properties\n");
        return 1;
    }

    cuda_err_chk(cudaSetDevice(cudaDevice));
    for (size_t i = 0 ; i < 1; i++)
        ctrls[i] = new Controller(ctrls_paths[i], 1, cudaDevice, queueDepth, queueNum);
    fprintf(stdout, "controller created\n");
    
    char st[15];
    cuda_err_chk(cudaDeviceGetPCIBusId(st, 15, cudaDevice));
    fprintf(stdout, "cudaDevice pcie: %s\n", st);

    h_pc = new page_cache_t(page_size, n_pages, cudaDevice, ctrls[0][0], (uint64_t) 64, ctrls);
    fprintf(stdout, "finished creating cache\n Total Cache size (MBs): %.2f\n", ((float)total_cache_size/(1024*1024)));

     //QueuePair* d_qp;
    d_pc = (page_cache_d_t*) (h_pc->d_pc_ptr);
    printf("n_tsteps: %lu, n_telem: %llu, n_pages:%llu\n", n_tsteps, n_telem, n_pages);

    return 0;
}

int nvme_dev_write()
{
    // strat write
    for (uint32_t cstep =0; cstep < n_tsteps; cstep++) {
        uint64_t cpysize = std::min(total_cache_size, sb_in.st_size-s_offset);
        printf("cstep: %lu  s_offset: %llu   cpysize: %llu pcaddr:%p, block size: %llu, grid size: %llu\n", cstep, s_offset, cpysize, h_pc->pdt.base_addr, b_size, g_size);

        cuda_err_chk(cudaMemcpy(h_pc->pdt.base_addr, map_in+s_offset+0, cpysize, cudaMemcpyHostToDevice));

        cudaEventCreate(&start_write); 
        cudaEventCreate(&stop_write);
        cudaEventRecord(start_write, 0);
        sequential_access_kernel<<<g_size, b_size>>>(h_pc->pdt.d_ctrls, d_pc, page_size, n_threads, //d_req_count,
        1, 1, WRITE, s_offset, 0);

        cuda_err_chk(cudaDeviceSynchronize());
        cudaEventRecord(stop_write, 0);
        cudaEventSynchronize(stop_write);

        float wcompleted = 100*(total_cache_size*(cstep+1))/(sb_in.st_size);
        cudaEventElapsedTime(&welapsed, start_write, stop_write);
        std::cout << "Write Completed:" << wcompleted << "%   Write Time:" <<welapsed << "ms" << std::endl;

         s_offset = s_offset + cpysize;
    }

    return 0;
}

uint64_t nvme_dev_read(void* read_offset, uint64_t idNUM, uint64_t read_size)
{
    // start read
    uint64_t* local_read_offset = (uint64_t*)read_offset;
    // n_tsteps = ceil((float)(read_size)/(float)total_cache_size);
    // for (uint32_t cstep =0; cstep < n_tsteps; cstep++) {
    // for (uint32_t cstep = 0; cstep < idNUM; cstep ++) {
        // uint64_t cpysize = std::min(total_cache_size, read_size);

        uint64_t* d_assignment;
        threadNum = idNUM;
        g_size = (threadNum + b_size - 1)/b_size;
        n_threads = b_size * g_size;
        cuda_err_chk(cudaMalloc(&d_assignment, n_threads*sizeof(uint64_t)));
        cuda_err_chk(cudaMemcpy(d_assignment, local_read_offset,  n_threads*sizeof(uint64_t), cudaMemcpyHostToDevice));

        // uint64_t data_read_offset = local_read_offset[cstep] * read_size;
        // printf("cstep: %lu  s_offset: %llu   cpysize: %llu pcaddr:%p, block size: %llu, grid size: %llu\n", cstep, s_offset, cpysize, h_pc->pdt.base_addr, b_size, g_size);

        cuda_err_chk(cudaMemset(h_pc->pdt.base_addr, 0, total_cache_size));
        
        cudaEventCreate(&start_read); 
        cudaEventCreate(&stop_read);
        cudaEventRecord(start_read, 0);
        // sequential_access_kernel<<<g_size, b_size>>>(h_pc->pdt.d_ctrls, d_pc, page_size, n_threads, //d_req_count,
        // 1, 1, READ, data_read_offset, 0);

        random_access_kernel<<<g_size, b_size>>>(h_pc->pdt.d_ctrls, d_pc, page_size,
        n_threads, 1, d_assignment, 1, READ);
        
        cuda_err_chk(cudaDeviceSynchronize());

        cudaEventRecord(stop_read, 0);
        cudaEventSynchronize(stop_read);
        
        // float rcompleted = 100*(cpysize*(cstep+1))/(read_size);
        cudaEventElapsedTime(&relapsed, start_read, stop_read);
        // std::cout << "Read Completed:" << rcompleted << "%   Read Time:" <<relapsed << "ms" << std::endl;
        std::cout << "Read Time:" <<relapsed << "ms" << std::endl;

        // local_read_offset = local_read_offset + cpysize;

        // printf("cuda addr: %p\n", h_pc->pdt.base_addr);
    // }

    return reinterpret_cast<uint64_t>(h_pc->pdt.base_addr);
}

int free_dev()
{
    for (size_t i = 0; i < 1; i++)
            delete ctrls[i];
    delete h_pc;

    return 0;
}
