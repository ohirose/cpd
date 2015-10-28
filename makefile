all:
	gcc -Wall -o cpd -O3 util.c rot.c affine.c cpd.c main.c -framework accelerate
