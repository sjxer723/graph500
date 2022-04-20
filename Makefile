all:
	gn gen -C out
	ninja -C out

format:
	gn format *.gn src/*.gn test/*.gn
	
clean:
	rm -rf out/

