CC=	g++
OPT=	-O3 -I include
SRCS=	GERMLINE_0001.cpp GERMLINE.cpp Share.cpp Chromosome.cpp ChromosomePair.cpp HMIndividualsExtractor.cpp MarkerSet.cpp Individual.cpp Individuals.cpp InputManager.cpp MatchFactory.cpp MatchesBuilder.cpp NucleotideMap.cpp PEDIndividualsExtractor.cpp Match.cpp PolymorphicIndividualsExtractor.cpp SNP.cpp SNPPositionMap.cpp SNPs.cpp
OBJS=	GERMLINE_0001.o GERMLINE.o Chromosome.o Share.o ChromosomePair.o HMIndividualsExtractor.o MarkerSet.o Individual.o Individuals.o InputManager.o MatchFactory.o MatchesBuilder.o NucleotideMap.o PEDIndividualsExtractor.o Match.o PolymorphicIndividualsExtractor.o SNP.o SNPPositionMap.o SNPs.o
MAIN=	germline
BMATCH=	parse_bmatch

all: clean germline bmatch test

bmatch:
	$(CC) $(BMATCH).cpp -o $(BMATCH)

$(OBJS): $(SRCS)
	$(CC) $(OPT) -c $*.cpp

germline: $(OBJS)
	$(CC) $(OPT) -o $(MAIN) $(OBJS)

clean:
	-rm -f *.o $(MAIN) $(BMATCH) test/generated.match test/generated.log test/generated.err test/generated.out
test: test_plink

test_plink:
	-mkdir -p test/output
	-@rm -f test/test_outputs/generated.match test/test_outputs/generated.log test/output/generated.err test/output/generated.out
	-@./$(MAIN) -silent -bits 50 -min_m 1 -err_hom 2 -old_samples test/old_humans -new_samples test/new_humans -err_het 0 < test/restricted.run > test/output/restricted.out 2> test/output/restricted.err | echo -e "---\nRunning Test Case\n---"
	-@./$(MAIN) -silent -bits 50 -min_m 1 -err_hom 2 -err_het 0 < test/test.run > test/output/generated.out 2> test/output/generated.err | echo -e "---\nRunning Test Case\n---"
	diff -q -s test/expected.match test/output/generated.match
	diff -q -s test/restricted.match test/output/restricted.match
