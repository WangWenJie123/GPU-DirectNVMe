#ifndef __READ_WRITE_SO__
#define __READ_WRITE_SO__

#include <stdio.h>
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

const char* const ctrls_paths[] = {"/dev/libnvm0"};
std::vector<Controller*> ctrls(1);
uint64_t queueNum = 32;
uint64_t queueDepth = 64;
uint64_t threadNum = 4096;

cudaEvent_t start_read, stop_read, start_write, stop_write;
float relapsed = 0, welapsed = 0;

uint64_t b_size = 64;//64;
uint64_t g_size = (threadNum + b_size - 1)/b_size;//80*16;
uint64_t n_threads = b_size * g_size;

// data cache memory size in cuda is page_size * n_pages=32MB
uint64_t page_size = 4096;
uint64_t n_pages = 8192;
uint64_t total_cache_size = (page_size * n_pages);

uint32_t n_tsteps = 0;  
uint64_t n_telem = 0; 
uint64_t s_offset = 0; 

page_cache_t* h_pc = NULL;
page_cache_d_t* d_pc = NULL;

struct stat sb_in;
void* map_in;

extern "C" {
 
  int dev_set(uint32_t cudaDevice, void* src_in);

  uint64_t nvme_dev_read(uint64_t read_offset, uint64_t read_size);

  int nvme_dev_write();

  int free_dev();
 
}

#endif
