CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11
LDFLAGS = -larmadillo

EXEC = main.exe
OBJ = main.o random.o popolazione.o individuo.o

$(EXEC): $(OBJ)
	$(CC) $(OBJ) -o $(EXEC) $(LDFLAGS)

main.o: main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)

random.o: random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

popolazione.o: popolazione.cpp popolazione.h individuo.h random.h
	$(CC) -c popolazione.cpp -o popolazione.o $(CFLAGS)

individuo.o: individuo.cpp individuo.h random.h
	$(CC) -c individuo.cpp -o individuo.o $(CFLAGS)

clean:
	rm -f *.o $(EXEC) seed.out
