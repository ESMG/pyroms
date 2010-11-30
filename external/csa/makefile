#
# makefile.in - Input makefile for the configure script.
#

SHELL = /bin/sh

#### Start of system configuration section. ####

prefix = /usr/local
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
libdir = /usr/local/lib
includedir = ${prefix}/include

INSTALLDIRS =\
$(bindir)\
$(libdir)\
$(includedir)

DISTDIR = ~/pub_web

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL}
INSTALL_DATA = ${INSTALL} -m 644

#### End of system configuration section. ####

CC = gcc
CFLAGS = -g -O2 -Wall -pedantic
LDFLAGS = 

AR = ar
ARFLAGS = cru

MLIB = -lm

SRC =\
csa.c\
minell.c\
svd.c

# I do not put minell.c into the library to avoid conflicts with libnn
#
LIBSRC =\
csa.c\
svd.c

HEADERS =\
csa.h\
minell.h\
nan.h\
svd.h\
version.h

PROGRAMS =\
csabathy

TESTS=\
minell\
lsq_test\
svd_test

LIBOBJS = ${LIBSRC:.c=.o}

all: libcsa.a $(PROGRAMS)

csabathy: $(SRC)
	$(CC) -o csabathy $(CFLAGS) -DCSA_STANDALONE -I. $(SRC) $(LDFLAGS) $(MLIB)

standalone: override LDFLAGS+=-static
standalone: $(PROGRAMS)
	strip $(PROGRAMS)

libcsa.a: $(LIBOBJS)
	$(AR) $(ARFLAGS) libcsa.a $(LIBOBJS)
	chmod go+r libcsa.a

tests: $(TESTS)

minell: minell.c
	$(CC) -o minell minell.c -DME_STANDALONE $(CFLAGS) -I. $(LDFLAGS) $(MLIB)

lsq_test: svd.c svd.h
	$(CC) -o lsq_test $(CFLAGS) -DLSQ_TEST svd.c $(LDFLAGS) $(MLIB)

svd_test: svd.c svd.h
	$(CC) -o svd_test $(CFLAGS) -DSVD_TEST svd.c $(LDFLAGS) $(MLIB)

installdirs:
	$(SHELL) mkinstalldirs $(INSTALLDIRS)

install: all installdirs
	for i in libcsa.a; do \
	   $(INSTALL_DATA) $$i $(libdir)/$$i; \
	done

	for i in csa.h; do \
	   $(INSTALL_DATA) $$i $(includedir); \
	done

	for i in csabathy; do \
	   fname=`basename $$i`; \
	   $(INSTALL_PROGRAM) $$i $(bindir); \
	done
clean:	
	rm -f *.o libcsa.a $(PROGRAMS) $(TESTS) *~

configclean:
	rm -f config.h makefile config.cache config.status config.log

distclean: clean configclean

indent:
	indent -T FILE -T point -T csa -T square -T triangle -T minell $(SRC) $(HEADERS)
	rm -f *~

dist:
	WDIR=`pwd`;\
	DIR=`basename $$WDIR`;\
	cp -r $$WDIR /tmp;\
	cd /tmp/$$DIR/examples/1;\
	make clean;\
	rm -f example.pdf;\
	cd /tmp/$$DIR/examples/2;\
	make clean;\
	rm -f example.pdf;\
	cd /tmp/$$DIR/examples/3;\
	make clean;\
	rm -f example.pdf;\
	cd /tmp/$$DIR;\
	rm -rf `find /tmp/$$DIR -name CVS`;\
	if [ -e makefile ] ; then\
		make distclean;\
	fi;\
	/tmp/"$$DIR"/configure;\
	make standalone;\
	chmod a+r+x $(PROGRAMS);\
	mv -f $(PROGRAMS) $(DISTDIR);\
	make distclean;\
	cd $(DISTDIR);\
	mv -f $$DIR".tar.gz" $$DIR".tar.gz.prev";\
	gzip -f $(PROGRAMS);\
	cd /tmp;\
	tar czvf $$DIR".tar.gz" $$DIR;\
	chmod a+r $$DIR".tar.gz";\
	mv -f $$DIR".tar.gz" $(DISTDIR)/$$DIR".tar.gz";\
	rm -rf $$DIR
