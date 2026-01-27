NVCC = nvcc
NVCC_FLAGS = -O3 -Xcompiler -fPIC -shared

TARGET = libweight_calc.so
SRC = weight_calc.cu

$(TARGET): $(SRC)
	$(NVCC) $(NVCC_FLAGS) -o $@ $<

clean:
	rm -f $(TARGET) 