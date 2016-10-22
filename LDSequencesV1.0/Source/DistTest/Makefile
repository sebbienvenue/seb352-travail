SRC=BaseDefinitions.cpp matrix.cpp SimSeqs.cpp Distributions.cpp main.cpp
OBJS=$(SRC:.cpp=.o)
CC=c++  
CCFLAGS=-g3 -Wall

LNKFLAGS=-g3 -Wall
TARGET=disttest


all: $(TARGET)

.cpp.o:
	$(CC) $(CCFLAGS) -c $<
.c.o:
	$(CC) $(CCFLAGS) -c $<

$(TARGET): $(OBJS)
	$(CC) $(LNKFLAGS) -o $(TARGET) $(OBJS)

clean:
	rm -rf $(OBJS) $(TARGET)

