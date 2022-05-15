all: main.o
	gcc main.o -o main -lm
main.o: main.c stb_image.h stb_image_write.h
	gcc -c  main.c
