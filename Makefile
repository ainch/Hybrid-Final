EXEC = 2D_Hybrid_PIC

CC = gcc
NVCC = nvcc
OPENMP=-fopenmp
OPENMP_CUDA=-Xcompiler -fopenmp -lgomp
# CFLAGS = -O3
CFLAGS =
CUFLAGS= -arch=sm_75 -use_fast_math --ptxas-options=-v
# CUFLAGS= -gencode arch=compute_52,code=sm_52 -use_fast_math --ptxas-options=-v
LFLAGS = -lm
CUDA_DIR=/usr/local/cuda-10.2
IDIR = -I$(CUDA_DIR)/include -Iinc
LIBS  = -L$(CUDA_DIR)/lib64 -lm -ldl -lcurand -lcublas -lcusparse -lcusolver -lcufft #-lefence
SRC_DIR := ./src
OBJ_DIR := ./obj


C_OBJS = $(OBJ_DIR)/main.o \
          $(OBJ_DIR)/parson.o \
          $(OBJ_DIR)/start.o \
          $(OBJ_DIR)/PhysicsData.o \
          $(OBJ_DIR)/load.o \
          $(OBJ_DIR)/Tecplot.o \
          $(OBJ_DIR)/Particle.o
#          $(OBJ_DIR)/savedata.o \
#          $(OBJ_DIR)/mccTools.o \
#          $(OBJ_DIR)/fluid.o \
#          $(OBJ_DIR)/fdtd.o

CU_OBJS = $(OBJ_DIR)/cuda_main.o \
          $(OBJ_DIR)/cuda_start.o 
#          $(OBJ_DIR)/cuda_run.o \
#          $(OBJ_DIR)/cuda_field.o \
#          $(OBJ_DIR)/cuda_move.o \
#          $(OBJ_DIR)/cuda_rand.o \
#          $(OBJ_DIR)/cuda_mccAr_cell.o \
#          $(OBJ_DIR)/cuda_diag.o \
#          $(OBJ_DIR)/cuda_sortboundary.o \
#          $(OBJ_DIR)/cuda_deposit.o \
#          $(OBJ_DIR)/cuda_beam.o

C_H_FILES = inc/xypic.h \
            inc/main.h \
            inc/start.h \
            inc/def.h \
            inc/PhysicsData.h \
            inc/load.h \
            inc/Tecplot.h \
            inc/Particle.h

CU_H_FILES = inc/xypic.cuh \
             inc/cuda_main.cuh \
             inc/cuda_start.cuh \
             inc/helper_cuda.h \
             inc/main.h \
             inc/interop.cuh \
             inc/def.h

all: $(OBJ_DIR) $(EXEC)

debug: CFLAGS += -g -G
debug: CUFLAGS += -g -G
debug: $(OBJ_DIR) $(EXEC)

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(C_H_FILES)
	$(NVCC) $(CFLAGS) $(CUFLAGS) $(IDIR) $(OPENMP_CUDA) $(LIBS) -c -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu $(CU_H_FILES)
	$(NVCC) $(CFLAGS) $(CUFLAGS) $(IDIR) $(OPENMP_CUDA) $(LIBS) -c -o $@ $<

$(EXEC): $(C_OBJS) $(CU_OBJS)
	$(NVCC) $(CFLAGS) $(CUFLAGS) $(IDIR) $(OPENMP_CUDA) $(LIBS) -o $@ $^

clean:
	@rm -rf $(OBJ_DIR) Potential.txt $(EXEC)

tar:
	@tar cvzf ./2D_PIC_cuda_com.tgz ./Viewer ./Trace *.c *.h makefile *.inp *.cu