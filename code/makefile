OBJS = Measurements.cpp GenHam_Heis.o  Lanczos_07.o lapack.o graphs.o 
CC = g++
CFLAGS = -O2 
#CFLAGS = -O2 -arch x86_64
#LIBS = -lm -framework veclib
LIBS = -framework Accelerate

2dEE.out: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o 2dEE.out $(LIBS)

Measurements.o : Measurements.cpp GenHam.h Lanczos_07.h lapack.h simparam.h 
	$(CC) $(CFLAGS) -c Measurements.cpp

GenHam_Heis.o: GenHam_Heis.cpp GenHam.h Lanczos_07.h entropy.h
	$(CC) $(CFLAGS) -c GenHam_Heis.cpp

Lanczos_07.o: Lanczos_07.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c Lanczos_07.cpp

lapack.o: lapack.cpp lapack.h 
	$(CC) $(CFLAGS) -c lapack.cpp

graphs.o: graphs.cpp graphs.h
	$(CC) $(CFLAGS) -c graphs.cpp

clean :
	rm *.o
