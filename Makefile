CC = gcc
LFLAGS = -lm
CFLAGS = -g
OBJ = sph_collapse.o
EXEC = sph_collapse.out

Main: $(OBJ)
	$(CC) -Wall -o $(EXEC) $(OBJ) $(CFLAGS) $(LFLAGS)

runmain: Main
	./$(EXEC)

