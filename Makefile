nbody: nbody.o nbody.h
	gcc -o nbody nbody.o -lm -lX11 -lplplotd

nbody.o: nbody.c nbody.h
	gcc -c nbody.c
