OBJS = Clock.o VolumeScene.o Main.o
CC = clang++
CFLAGS = -Wall -c $(DEBUG) -ferror-limit=1
LDFLAGS = -Wall $(DEBUG) -framework GLUT -framework OpenGL -framework Cocoa

volumeray: $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o volumeray

Clock.o: Clock.cpp
	$(CC) $(CFLAGS) Clock.cpp

VolumeScene.o: VolumeScene.cpp
	$(CC) $(CFLAGS) VolumeScene.cpp

Main.o: Main.cpp
	$(CC) $(CFLAGS) Main.cpp

.PHONY: clean
clean:
	rm -f *.o volumeray
	rm -f *~
