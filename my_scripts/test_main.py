from ctypes import *

nvme_rwlib = CDLL("/home/wwj/GNN_Training_Acceleration/Ref_projs/bam/build/lib/libread_write.so")

set_val = nvme_rwlib.dev_set(0)
write_val = nvme_rwlib.nvme_dev_write()
read_val = nvme_rwlib.nvme_dev_read()
free_val = nvme_rwlib.free_dev()