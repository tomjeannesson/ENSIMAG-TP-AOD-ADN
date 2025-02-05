# Makefile qui genere l'executable distanceEdition et fait des tests de verification
#
#
CC=gcc
LATEXC=pdflatex
DOCC=doxygen
CFLAGS=-g -Wall 

REFDIR=.
SRCDIR=$(REFDIR)/src
BINDIR=$(REFDIR)/bin
DOCDIR=$(REFDIR)/doc
TESTDIR=$(REFDIR)/tests
REPORTDIR=$(REFDIR)/report

LATEXSOURCE=$(wildcard $(REPORTDIR)/*.tex)
CSOURCE=$(wildcard $(SRCDIR)/*.c)
PDF=$(LATEXSOURCE:.tex=.pdf)

all: binary report doc binary_perf

binary: $(BINDIR)/distanceEdition

binary_perf: $(BINDIR)/distanceEdition-perf

$(BINDIR)/distanceEdition-perf: $(SRCDIR)/distanceEdition.c $(BINDIR)/Needleman-Wunsch-recmemo.o
	$(CC) $(OPT) -D__PERF_MESURE__ -I$(SRCDIR) -o $(BINDIR)/distanceEdition-perf $(BINDIR)/Needleman-Wunsch-recmemo.o $(SRCDIR)/distanceEdition.c 

report: $(PDF) 

doc: $(DOCDIR)/index.html


$(BINDIR)/distanceEdition: $(SRCDIR)/distanceEdition.c $(BINDIR)/Needleman-Wunsch-recmemo.o
	$(CC) $(OPT) -I$(SRCDIR) -o $(BINDIR)/distanceEdition $(BINDIR)/Needleman-Wunsch-recmemo.o $(SRCDIR)/distanceEdition.c -lm

$(BINDIR)/Needleman-Wunsch-recmemo.o: $(SRCDIR)/Needleman-Wunsch-recmemo.h $(SRCDIR)/Needleman-Wunsch-recmemo.c $(SRCDIR)/characters_to_base.h
	$(CC) $(OPT) -I$(SRCDIR) -c  -o $(BINDIR)/Needleman-Wunsch-recmemo.o $(SRCDIR)/Needleman-Wunsch-recmemo.c
	
$(BINDIR)/extract-fasta-sequences-size: $(SRCDIR)/extract-fasta-sequences-size.c
	$(CC) $(OPT) -I$(SRCDIR) -o $(BINDIR)/extract-fasta-sequences-size $(SRCDIR)/extract-fasta-sequences-size.c

clean:
	rm -rf $(DOCDIR) $(BINDIR)/* $(REPORTDIR)/*.aux $(REPORTDIR)/*.log  $(REPORTDIR)/rapport.pdf $(TESTDIR)/*.output $(TESTDIR)/cachegrind.out.*

#$(BINDIR)/distanceEdition: $(CSOURCE)
#	$(CC) $(CFLAGS)  $^ -o $@ 

$(BINDIR)/distanceEditiondebug: $(CSOURCE)
	$(CC) $(CFLAGS)  $^ -o $@ -DDEBUG

%.pdf: $(LATEXSOURCE)
	$(LATEXC) -output-directory $(REPORTDIR) $^ 

$(DOCDIR)/index.html: $(SRCDIR)/Doxyfile $(CSOURCE)
	$(DOCC) $(SRCDIR)/Doxyfile


test: clean $(BINDIR)/distanceEdition $(TESTDIR)/Makefile-test
	cd $(TESTDIR) ; make -f Makefile-test all 
	
test-valgrind: clean $(BINDIR)/distanceEdition $(TESTDIR)/Makefile-test
	cd $(TESTDIR) ; make -f Makefile-test all-valgrind
	
.PHONY: all doc bin report 

