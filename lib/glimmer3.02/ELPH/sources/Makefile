CLASSDIR := .

# Directories to search for header files
SEARCHDIRS := -I- -I${CLASSDIR}


SYSTYPE :=     $(shell uname)

# C compiler

CC      := g++
CFLAGS  = -Wall ${SEARCHDIRS} -fno-exceptions -fno-rtti -D_REENTRANT -g

%.o : %.c
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cc
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.C
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cpp
	${CC} ${CFLAGS} -c $< -o $@

%.o : %.cxx
	${CC} ${CFLAGS} -c $< -o $@

# C/C++ linker

LINKER    := g++
LDFLAGS    =
LOADLIBES := 

.PHONY : all
all:    elph

elph:  ./elph.o ${CLASSDIR}/motif.o ${CLASSDIR}/GBase.o ${CLASSDIR}/GString.o ${CLASSDIR}/GArgs.o
	${LINKER} ${LDFLAGS} -o $@ ${filter-out %.a %.so, $^} ${LOADLIBES}

# target for removing all object files

.PHONY : tidy
tidy::
	@${RM} core elph *.o ${CLASSDIR}/GBase.o ${CLASSDIR}/GString.o ${CLASSDIR}/GArgs.o

# target for removing all object files

.PHONY : clean
clean:: tidy
	@${RM} core elph *.o ${CLASSDIR}/GBase.o ${CLASSDIR}/GString.o ${CLASSDIR}/GArgs.o


