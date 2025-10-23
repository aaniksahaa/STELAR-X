NVCC = nvcc
NVCC_FLAGS = -O3 -Xcompiler -fPIC -shared

# Memory-efficient CUDA library
TARGET = libmemory_efficient_weight_calc.so
SRC = memory_efficient_weight_calc.cu

$(TARGET): $(SRC)
	$(NVCC) $(NVCC_FLAGS) -o $@ $<

clean:
	rm -f $(TARGET) 