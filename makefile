all:
	gcc -Wall -o cpd    -O3 cpd.c    util.c -framework accelerate
	gcc -Wall -o affine -O3 affine.c util.c -framework accelerate
	gcc -Wall -o rot    -O3 rot.c    util.c -framework accelerate
