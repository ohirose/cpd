all:
	gcc -Wall -o cpd -O3 cpd.c util.c rigid.c affine.c info.c main.c -llapack -lm
