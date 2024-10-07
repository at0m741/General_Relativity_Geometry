NAME = Schwarzschild_geometry

CC = clang

SRCS = Matrix.c \
	   Metric.c \
	   Tensor.c \
	   Utils.c \
	   main.c

CFLAGS = -O3 -lm -ffast-math -mavx2 -mfma -march=native 

all: $(NAME)

$(NAME): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $(NAME)

clean:
	rm -f $(NAME)

re: clean all

.PHONY: all clean re

