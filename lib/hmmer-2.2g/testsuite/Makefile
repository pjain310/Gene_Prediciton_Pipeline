# Generated automatically from Makefile.in by configure.
################################################################
# Makefile for HMMER testsuite
# CVS $Id: Makefile.in,v 1.13 2001/08/05 23:44:43 eddy Exp $
##########
# HMMER - Biological sequence analysis with profile HMMs
# Copyright (C) 1992-1999 Washington University School of Medicine
# All Rights Reserved
# 
#     This source code is distributed under the terms of the
#     GNU General Public License. See the files COPYING and LICENSE
#     for details.
###########

CC        = gcc
CFLAGS    = -g -O2
MDEFS     =  -DHMMER_THREADS  -DHAVE_PTHREAD_ATTR_SETSCOPE=1 -DHAVE_NTOHS=1 -DHAVE_NTOHL=1 -DHAVE_HTONS=1 -DHAVE_HTONL=1 

# Configuration for optional pthreads multiprocessor support
#
PTHREAD_LIBS   = 
PTHREAD_CFLAGS = -pthread

SHELL     = /bin/sh
MYLIBS    = -lhmmer -lsquid
LIBS      =  -lm

SHIVA = alignalign_test\
	evd_test\
	masks_test\
	parsingviterbi_test\
	tophits_test\
	trace_test\
	viterbi_exercise\
	weeviterbi_test

#######
## Targets defining how to make Shiva executables.
#######

.c.o:
	$(CC) $(CFLAGS) $(PTHREAD_CFLAGS) $(MDEFS) -I../squid -I../src -c $<		

all: 	$(SHIVA)

$(SHIVA): %: %.o
	$(CC) $(CFLAGS) $(PTHREAD_CFLAGS) $(MDEFS) -o $@ -L../squid -L../src $@.o $(MYLIBS) $(PTHREAD_LIBS) $(LIBS)

#######
## `make check` actually runs the tests.
#######

check: $(SHIVA)
	@echo 
	@echo Running compiled Shiva exercises:
	@echo Warning: some tests may take several minutes to complete. 
	@for shivatest in $(SHIVA); do\
	   if ./$$shivatest; then\
              echo $$shivatest: ok;\
           else\
              echo $$shivatest: FAILED;\
           fi;\
	done
	@echo
	@echo Running scripted Shiva exercises:
	@echo Warning: This also may take several minutes.
	perl ./Optiontests.pl
	sh ./Exercises.sh

#######
## Miscellaneous
#######

clean:
	-rm -f *.o *~ Makefile.bak core $(SHIVA) TAGS gmon.out

distclean:
	make clean
	-rm -f Makefile

TAGS:
	etags -t *.c *.h Makefile.in


