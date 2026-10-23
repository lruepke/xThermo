CC = g++
PROM = waterEOS
binDIR=.
# CFLAGS = -I../src
DEPS = $(shell find . -name "*.H")
SRC = $(shell find . -name "*.cpp")
OBJ = $(SRC:%.cpp=%.o)

$(PROM): $(OBJ)
	$(CC) -o $(binDIR)/$(PROM) $(OBJ) $(CFLAGS)  -L${FOAM_USER_LIBBIN} -lfreesteam

%.o:%.cpp $(DEPS)
	$(CC) -c $< -o $@ $(CFLAGS) 

clean:
	rm *.o