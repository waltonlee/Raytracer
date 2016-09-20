all: *.cpp
	g++ *.cpp -L/usr/local/lib -I/usr/local/Cellar/glui/2.3.6/include/ -framework OpenGL -framework GLUT -lglui -w -o run
