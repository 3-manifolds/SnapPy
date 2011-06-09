@OS = $(shell echo %OS%)

ifeq ($(OS),Windows_NT)
	OUTPUT = Twister.exe
	OS_ENVIRONMENT = 0
else
	OUTPUT = Twister.out
	OS_ENVIRONMENT = 1
endif

CC = g++
CFLAGS = -D OS_ENVIRONMENT=$(OS_ENVIRONMENT) -W -Wall
MAIN = twister_main.cpp twister.cpp parsing.cpp global.cpp
HEAD = twister.h parsing.h global.h

all : $(MAIN) $(HEAD)
	$(CC) -o $(OUTPUT) $(MAIN) $(CFLAGS)

.PHONY : all
