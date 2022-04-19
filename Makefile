all:
	tools/gn gen -C out
	ninja -C out

format:
	tools/gn format *.gn src/*.gn test/*.gn
	
clean:
	rm -rf out/

