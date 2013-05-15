OBJS = NLC_2D_EE.cpp GenHam.o  Lanczos_07.o lapack.o graphs.o 
CC = g++
#CFLAGS = -O2 
CFLAGS = -O2 -arch x86_64
#LIBS = -lm -framework veclib
LIBS = -framework Accelerate

2dEE.out: $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o 2dEE.out $(LIBS)

NLC_2D_EE.o : NLC_2D_EE.cpp GenHam.h Lanczos_07.h lapack.h simparam.h 
	$(CC) $(CFLAGS) -c NLC_1D_TFIM.cpp

GenHam.o: GenHam.cpp GenHam.h Lanczos_07.h entropy.h
	$(CC) $(CFLAGS) -c GenHam.cpp

Lanczos_07.o: Lanczos_07.cpp GenHam.h Lanczos_07.h
	$(CC) $(CFLAGS) -c Lanczos_07.cpp

lapack.o: lapack.cpp lapack.h 
	$(CC) $(CFLAGS) -c lapack.cpp

graphs.o: graphs.cpp graphs.h
	$(CC) $(CFLAGS) -c graphs.cpp

clean :
	rm *.o
