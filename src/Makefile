CC = gcc
CFLAGS = -Wall -g
LDFLAGS = -lm

#source files
SRC = raytracer.c
OBJ = $(SRC:.c=.o)
EXEC = raytracer

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $(EXEC)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

#run program with input.txt

run: $(EXEC)
	./$(EXEC) input.txt

#Clean up
clean:
	rm -f $(OBJ) $(EXEC)
