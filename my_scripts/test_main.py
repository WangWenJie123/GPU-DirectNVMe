from ctypes import *
import numpy as np
import cupy as cp
import torch
from cupy.cuda import MemoryPointer, UnownedMemory

nvme_rwlib = CDLL("/home/wwj/GNN_Training_Acceleration/Ref_projs/bam/build/lib/libread_write.so")
nvme_rwlib.nvme_dev_read.restype = c_uint64

feature_data_path = b"/home/wwj/GNN_Training_Acceleration/Ginex/dataset/ogbn-papers100M-ginex/features.dat"

set_val = nvme_rwlib.dev_set(0, feature_data_path)

# write_val = nvme_rwlib.nvme_dev_write()

feature_size = 128*4

# create id array and pass to cuda api
id_len = 1000
id_list = []
for i in range(id_len):
    id_list.append(i)
id_list_ptr = (c_uint64*id_len)(*id_list)
read_val = nvme_rwlib.nvme_dev_read(id_list_ptr, id_len, feature_size)

mem = UnownedMemory(read_val, feature_size*id_len, None)
memptr = MemoryPointer(mem, 0)

cuda_data = cp.ndarray(shape=int(feature_size*id_len/4), dtype=np.float32, memptr=memptr)
    
torch_data = torch.utils.dlpack.from_dlpack(cuda_data.toDlpack()).float()

for i in range(30):
    print("torch_data: ", torch_data[i])
free_val = nvme_rwlib.free_dev()