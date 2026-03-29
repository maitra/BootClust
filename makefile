CC = gcc

CFLAGS = -Wall -O4 -std=c99 -pedantic 
OBJ = src/BootCluster.o src/EigValDec.o src/hc.o src/hclass.o src/hclassify.o src/inverse.o src/kmeans.o src/mat_vec-VM.o src/order.o src/quantile.o src/sorted.o src/svd.o src/kmeans-initials_hclust.o src/eigens.o src/statutils.o src/readopt.o

BootClust: src/BC.c $(OBJ)
	$(CC) $(CFLAGS) src/BC.c -o BootClust $(OBJ) -lm -llapack -lRmath
clean:  
	rm -rf *~ src/*~ src/*.o BootClust
