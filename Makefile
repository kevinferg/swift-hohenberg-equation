CC = gcc
CFLAGS = -O2
LDFLAGS = -lm
TARGET = main.exe
SRCS = $(wildcard *.c)

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) $(LDFLAGS) -o $(TARGET)

clean:
	del /Q $(TARGET) 2>nul || rm -f $(TARGET)