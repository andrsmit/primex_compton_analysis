#
CERN_ROOT = /apps/cernlib/x86_64_rhel6/2005

all: ana_compton ana_compton_sys

ana_compton: analyze_compton.o Makefile ComptonAna.o
	g++ analyze_compton.o ComptonAna.o -o ana_compton -L/usr/lib64/root \
	`root-config --cflags --glibs`
#
ana_compton_sys: analyze_compton_sys.o Makefile ComptonAna.o ComptonAna_sys.o
	g++ analyze_compton_sys.o ComptonAna.o ComptonAna_sys.o -o ana_compton_sys -L/usr/lib64/root \
	`root-config --cflags --glibs`
#
ComptonAna.o: ComptonAna.cc Makefile $(wildcard ./*.h)
	@echo Compiling ComptonAna.cc; g++ -c ComptonAna.cc \
	`root-config --cflags --glibs`
#
ComptonAna_sys.o: ComptonAna_sys.cc Makefile $(wildcard ./*.h)
	@echo Compiling ComptonAna_sys.cc; g++ -c ComptonAna_sys.cc \
	`root-config --cflags --glibs`
#
analyze_compton.o: analyze_compton.cc Makefile $(wildcard ./*.h)
	@echo Compiling analyze_compton.cc; g++ -c analyze_compton.cc \
	`root-config --cflags --glibs`
#
analyze_compton_sys.o: analyze_compton_sys.cc Makefile $(wildcard ./*.h)
	@echo Compiling analyze_compton_sys.cc; g++ -c analyze_compton_sys.cc \
	`root-config --cflags --glibs`
#
