SRC=BaseDefinitions.cpp matrix.cpp SimSeqs.cpp main.cpp 
OBJS=$(SRC:.cpp=.o)
CC=c++  
CCFLAGS=-O2 -Wall

LNKFLAGS=-O2 -Wall
TARGET=SeqTest.x


all: $(TARGET)

.cpp.o:
	$(CC) $(CCFLAGS) -c $<
.c.o:
	$(CC) $(CCFLAGS) -c $<

$(TARGET): $(OBJS)
	$(CC) $(LNKFLAGS) -o $(TARGET) $(OBJS)

clean:
	rm -rf $(OBJS) $(TARGET)

