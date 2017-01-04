#!/bin/csh
#------------------------------------------------------------------------------
# Version 11-February-2002
#------------------------------------------------------------------------------
# To install octopus properly, you must follow these steps:
# make include
# make octopus
# make clean
#------------------------------------------------------------------------------
# verify that the required libraries are properly set
PGPDIR  = /usr/local/pgplot
X11DIR  = /usr/X11R6/lib
FIODIR  = /usr/local/cfitsio
FCOMPIL = g77 -O3 -g -Wall
CCOMPIL = gcc -O3 -g -Wall
#------------------------------------------------------------------------------
# Nothing SHOULD be modified below this comment line
#------------------------------------------------------------------------------
# macro definitions
FSOURCE = button.f buttqbr.f buttqcf.f buttqch.f buttqex.f \
          buttqit.f buttqpr.f buttqxb.f buttqyb.f buttqytext.f buttsbr.f \
          buttscf.f buttsch.f buttsex.f buttsit.f buttspr.f buttsxb.f \
          buttsyb.f buttsytext.f ifbutton.f rpgband.f rpgbegin.f rpgbegok.f \
          rpgenv.f rpgeras.f rpgerasb.f rpgerasw.f \
          iofunctions.f \
          octopus.f \
          subprece.f
FOBJECT = button.o buttqbr.o buttqcf.o buttqch.o buttqex.o \
          buttqit.o buttqpr.o buttqxb.o buttqyb.o buttqytext.o buttsbr.o \
          buttscf.o buttsch.o buttsex.o buttsit.o buttspr.o buttsxb.o \
          buttsyb.o buttsytext.o ifbutton.o rpgband.o rpgbegin.o rpgbegok.o \
          rpgenv.o rpgeras.o rpgerasb.o rpgerasw.o \
          iofunctions.o \
          octopus.o \
          subprece.o
# Default rule to create program
octopus:  $(FOBJECT)
#	$(FCOMPIL) -o $@ $(FOBJECT) -L$(PGPDIR) -L$(FIODIR) -L$(X11DIR) -lpgplot -lcfitsio -lX11 -lnsl -lsocket
	$(FCOMPIL) -o $@ $(FOBJECT) -L$(PGPDIR) -L$(FIODIR) -L$(X11DIR) -lpgplot -lcfitsio -lX11
# Target to clean object modules
clean:    $(FOBJECT)
	rm -f $(FOBJECT)
# Target to touch source modules
touch:
	touch $(FSOURCE)
# Target to create the file octopus_dir.inc
include:
	rm -f octopus_dir.inc
	echo "        CHARACTER*255 OCTOPUS_DIR" > octopus_dir.inc
	echo "        PARAMETER(OCTOPUS_DIR=" >> octopus_dir.inc
	echo "     +   '`pwd`')" >> octopus_dir.inc
# second level dependencies
.f.o: $(FSOURCE)
	$(FCOMPIL) -c $?
# definitions
.PRECIOUS: octopus
