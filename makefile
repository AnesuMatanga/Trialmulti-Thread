cxx = g++
cflags = -O3 -fopenmp
inc = -I./ -I"C:/SFML-2.6.1/include"
sfml-libs = -L"C:/SFML-2.6.1/lib" -lsfml-window -lsfml-graphics -lsfml-system

all: nbody-nographics nbody-graphics

nbody-nographics : nbody.cpp
	${cxx} nbody.cpp -o nbody-nographics ${cflags} ${inc}

nbody-graphics : nbody.cpp
	${cxx} nbody.cpp -o nbody-graphics ${cflags} ${inc} ${sfml-libs} -DGRAPHICS

clean:
	del nbody-graphics nbody-nographics