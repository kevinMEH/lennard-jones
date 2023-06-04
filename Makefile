all: build/$(file)
	./build/$(file)

build/$(file): $(file).c src/*.c
	gcc $(file).c src/*.c $(args) -O3 -o build/$(file)

unop:
	gcc $(file).c src/*.c $(args) -o build/$(file)
	./build/$(file)

clean:
	rm build/*