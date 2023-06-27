#ifndef __READ_WRITE_SO__
#define __READ_WRITE_SO__

#include <stdio.h>

extern "C" {
 
  int main(uint32_t cudaDevice);

  uint32_t nvme_dev_read();

  uint32_t nvme_dev_write();
 
}

#endif