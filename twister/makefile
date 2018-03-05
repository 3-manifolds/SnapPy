@OS = $(shell echo %OS%)

ifeq ($(OS),Windows_NT)
	OUTPUT = Twister.exe
	OS_ENVIRONMENT = 0
else
	OUTPUT = Twister.out
	OS_ENVIRONMENT = 1
endif

CC = g++
MAIN = ./lib/kernel/*.cpp
HEAD = ./lib/kernel/*.h
CFLAGS = -D OS_ENVIRONMENT=$(OS_ENVIRONMENT) -W -Wall

all : $(MAIN)
	$(CC) -o $(OUTPUT) $(MAIN) $(CFLAGS)

.PHONY : all
