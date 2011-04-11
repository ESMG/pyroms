all: grid.0 grid.1 grid.2 grid.3 grid.4 grid.5 grid.6

grid.0:
	../gridgen prm.0 -v
grid.1:
	../gridgen prm.1 -v
grid.2:
	../gridgen prm.2 -v
grid.3:
	../gridgen prm.3 -v
grid.4:
	../gridgen prm.4 -v
grid.5:
	../gridgen prm.5 -v
grid.6:
	echo "1" | awk -f mkprotogrid.6.awk > proto.6
	../gridgen prm.6 -v
	./addheader.sh proto.6 grid.6
clean:
	rm -f *~ grid.* rect.* sigmas.*
