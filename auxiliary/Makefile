#
#
# It is very general Makefile for program, 
# which is done only from one file.cc in dir
#
#

OB = $(patsubst %.cc, %.o, $(wildcard *.cc))
PR_NAME= gmv2mpeg 
CFLAGS= -g
$(PR_NAME): $(OB)
	@g++ -o $@ $(OB)
$(OB): $(wildcard *.cc) 
	@g++ $(CFLAGS) -c $*.cc
	@echo "Exe file is " $(PR_NAME)

.PHONY: clean
clean : 
	rm $(PR_NAME) $(OB)
