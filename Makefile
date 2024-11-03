CC=	g++
OPT=	-std=c++17 -O3 -I include
SRCS=	GERMLINE_0001.cpp GERMLINE.cpp Share.cpp Chromosome.cpp ChromosomePair.cpp HMIndividualsExtractor.cpp MarkerSet.cpp Individual.cpp Individuals.cpp InputManager.cpp MatchFactory.cpp MatchesBuilder.cpp NucleotideMap.cpp PEDIndividualsExtractor.cpp Match.cpp PolymorphicIndividualsExtractor.cpp SNP.cpp SNPPositionMap.cpp SNPs.cpp
OBJS=	GERMLINE_0001.o GERMLINE.o Chromosome.o Share.o ChromosomePair.o HMIndividualsExtractor.o MarkerSet.o Individual.o Individuals.o InputManager.o MatchFactory.o MatchesBuilder.o NucleotideMap.o PEDIndividualsExtractor.o Match.o PolymorphicIndividualsExtractor.o SNP.o SNPPositionMap.o SNPs.o
MAIN=	germline
BMATCH=	parse_bmatch

all: clean setup germline bmatch test

bmatch: setup
	$(CC) $(BMATCH).cpp -o bin/$(BMATCH)
	md5sum bin/$(BMATCH) > $(BMATCH).md5

$(OBJS): $(SRCS)
	$(CC) $(OPT) -c $*.cpp

germline: $(OBJS) setup
	$(CC) $(OPT) -o bin/$(MAIN) $(OBJS) -lstdc++fs
	md5sum bin/$(MAIN) > $(MAIN).md5

setup:
	-mkdir -p bin

clean:
	-rm -f *.o bin/$(MAIN) $(BMATCH) test/generated.match test/generated.log test/generated.err test/generated.out
test: test_plink

test_plink:
	-mkdir -p test/output
	-@rm -f test/output/*
	-@./bin/$(MAIN) -silent -bits 50 -min_m 1 -err_hom 2 -samples_to_compare_to test/old_humans -new_samples test/new_humans -err_het 0 < test/restricted.run > test/output/restricted.out 2> test/output/restricted.err | echo -e "---\nRunning Test Case\n---"
	-@./bin/$(MAIN) -silent -bits 50 -min_m 1 -err_hom 2 -err_het 0 < test/test.run > test/output/generated.out 2> test/output/generated.err | echo -e "---\nRunning Test Case\n---"
	diff -q -s test/expected.match test/output/generated.match
	diff -q -s test/restricted.match test/output/restricted.match
