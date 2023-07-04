cmake_minimum_required (VERSION 3.10)
project (libnvm-benchmarks)

set (THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)

find_package (CUDA 11.7)

# make_benchmark (readwrite-benchmark-module readwrite-bench "read_write_so.cu;")
make_benchmark (read_write "read_write_so.cu;")

make_benchmark_choice (readwrite-benchmark readwrite-benchmark-smartio readwrite-benchmark-module)
