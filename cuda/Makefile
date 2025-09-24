NVCC = nvcc
NVCC_FLAGS = -O3 -Xcompiler -fPIC -shared

TARGET = libweight_calc.so
TARGET_RANGE = librange_weight_calc.so
SRC = weight_calc.cu
RANGE_SRC = range_weight_calc.cu

all: $(TARGET) $(TARGET_RANGE)

$(TARGET): $(SRC)
	$(NVCC) $(NVCC_FLAGS) -o $@ $<

$(TARGET_RANGE): $(RANGE_SRC)
	$(NVCC) $(NVCC_FLAGS) -o $@ $<

clean:
	rm -f $(TARGET) $(TARGET_RANGE)

.PHONY: all clean 