#
CERN_ROOT = /apps/cernlib/x86_64_rhel6/2005

all: ana_compton

ana_compton: analyze_compton.o Makefile
	g++ analyze_compton.o -o ana_compton -L/usr/lib64/root \
	`root-config --cflags --glibs`
#
analyze_compton.o: analyze_compton.cc Makefile $(wildcard ./*.h)
	@echo Compiling analyze_compton.cc; g++ -c analyze_compton.cc \
	`root-config --cflags --glibs`
#
