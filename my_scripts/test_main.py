from ctypes import *

nvme_rwlib = CDLL("/home/wwj/GNN_Training_Acceleration/Ref_projs/bam/build/lib/libread_write.so")

set_val = nvme_rwlib.main(0)
# print("set_val: ", set_val)

# read_val = nvme_rwlib.nvme_dev_read()
# print("read_val: ", read_val)

# write_val = nvme_rwlib.nvme_dev_write()
# print("write_val: ", write_val)