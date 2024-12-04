# This assumes GNU make!
# You'll need to have installed the fftw3 and hdf5 libraries
# and might need to change the following two lines to point to them

# On Ubuntu, the following should work:
#CPPFLAGS=-O3 -Dreal=float -I/usr/include -DH5_USE_16_API
#LIBS=-L/usr/lib64/ -lhdf5 -lhdf5_hl -lfftw3 -lm

# On MacOS with fftw3 and hdf5 installed using Homebrew, the following should work:
CPPFLAGS=-O3 -Dreal=float -I/opt/homebrew/include -DH5_USE_16_API
LIBS=-L/opt/homebrew/lib -lhdf5 -lhdf5_hl -lfftw3 -lm

mestel2d : mestel2d.o 2d_bodies.o 2d_box.o h5_files.o
	$(CXX) -o mestel2d mestel2d.o 2d_bodies.o 2d_box.o h5_files.o $(LIBS)

clean :
	@rm mestel2d.o 2d_bodies.o 2d_box.o h5_files.o mestel2d
