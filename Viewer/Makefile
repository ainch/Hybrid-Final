EXEC = cuda_viewer

CC = gcc
NVCC = gcc
CFLAGS = -O3 
LFLAGS = -lm
Xgrafix_DIR=/usr/local/xgrafix
Xgrafix2_DIR=/usr/X11R6
IDIR = -I$(Xgrafix_DIR)/include -I$(Xgrafix2_DIR)/include/X11 -Iinc
LIBS  = -L/usr/X11R6/lib -L/usr/local/xgrafix/lib  -L/usr/local/lib -L/usr/lib64\
 -lXGC250 -ltk8.5 -ltcl8.5 -lXpm -lX11 -lm -ldl 
SRC_DIR := ./src
OBJ_DIR := ./obj
OBJ_DIR2 := ../obj
HEADER_DIR := ../inc

C_OBJS = $(OBJ_DIR)/Main_Viewer.o \
          $(OBJ_DIR)/Viewer.o \
          $(OBJ_DIR2)/parson.o \
          $(OBJ_DIR2)/start.o \
          $(OBJ_DIR2)/PhysicsData.o \
          $(OBJ_DIR2)/load.o \
          $(OBJ_DIR2)/Diagnostics.o \
          $(OBJ_DIR2)/Efield.o \
          $(OBJ_DIR2)/Fluid.o \
          $(OBJ_DIR2)/Tecplot.o \
          $(OBJ_DIR2)/Particle.o

C_H_FILES = inc/xypic.h \
            inc/Main_Viewer.h \
            inc/Viewer.h \
            $(HEADER_DIR)/start.h \
            $(HEADER_DIR)/def.h \
            $(HEADER_DIR)/PhysicsData.h \
            $(HEADER_DIR)/load.h \
            $(HEADER_DIR)/Efield.h \
            $(HEADER_DIR)/Fluid.h \
            $(HEADER_DIR)/Tecplot.h \
            $(HEADER_DIR)/Particle.h

all: $(OBJ_DIR) $(EXEC)

debug: CFLAGS += -g -G
debug: $(OBJ_DIR) $(EXEC)

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(C_H_FILES)
	$(CC) -c -o $@ $< $(CFLAGS) $(IDIR) $(LIBS)

$(EXEC): $(C_OBJS) 
	$(CC) -o $@ $^ $(CFLAGS) $(IDIR) $(LIBS) 

clean:
	@rm -rf $(OBJ_DIR) Potential.txt $(EXEC)

tar:
	@tar cvzf ./2D_PIC_cuda_com.tgz ./Viewer ./Trace *.c *.h makefile *.inp *.cu