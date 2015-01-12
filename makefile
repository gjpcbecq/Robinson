ARCHI=-m64
NORM=-pedantic


ER.o: ER.c ER.h
	gcc -c ER.c

testER.o: testER.c ER.h	
	gcc -c testER.c 
	
ER.dylib: ER.o ER.h
	gcc -shared $(ARCHI) $(NORM) -o ER.dylib ER.o 
	
test: ER.h ER.o testER.o
	gcc -o testER ER.o testER.o


