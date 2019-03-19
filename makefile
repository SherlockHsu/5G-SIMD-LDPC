SRCS = $(wildcard *.c)

OBJS = $(SRCS:.c = .o)

CC = icc

INCLUDES = -I./inc

LIBS = -lm -lpthread -lmkl_rt -fopenmp ./lib/libsimd_ldpc.a

CCFLAGS = -Wall -march=core-avx512 -std=c99 -g

OUTPUT = main

all:$(OUTPUT)

$(OUTPUT) : $(OBJS)
	$(CC) $^ -o $@ $(INCLUDES) $(LIBS)

%.o : %.c
	$(CC) -c $< $(CCFLAGS)

clean:
	rm -rf main *.o *.txt    #清除中间文件及生成文件

.PHONY:clean
