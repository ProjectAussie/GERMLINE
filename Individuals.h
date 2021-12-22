// Individuals.h: A collection of individuals

#ifndef INDIVIDUALS_H
#define INDIVIDUALS_H

#include "BasicDefinitions.h"
#include "Chromosome.h"
#include "Individual.h"
#include <set>
#include <ostream>
using namespace std;

class Individual;
class Individuals
{
public:

	// Individuals(): default constructor
	// Precondition: None.
	// Postcondition: individuals is empty.
	Individuals();
	~Individuals();

	// addIndividual(): adds an Individual object
	// Precondition: None.
    // Postcondition: ind has been added to individuals
	void addIndividual( Individual * ind );
	Individual * getIndividual ( size_t id ) { return pedigree[ id ]; }

	bool more();
	Individual* next();
	void begin();
	size_t size() { return pedigree.size(); }
	void initialize();
	void initializeOutputFileHandles(string chromosome);
	void print( ostream& );
	
	void freeMatches();
	void freeMarkers();
	void loadOldIndividuals(string f);
	void loadNewIndividuals(string f);
	bool isOld(string);
	bool isNew(string);
	bool hasRestrictions();

	bool useEmbarkRFGermlineOutput;
	string chromosome; // used for individual file handle outputs
	string individualOutputFolder; // top level folder for individual outputs

private:

	void permuteMarkerSet(Chromosome *, int, MarkerSet);
	// stores the individuals
	vector< Individual * > pedigree;
	size_t iter;

	long sets;
	set<string> samples_to_compare_to; 
	set<string> new_samples;
};

#endif

// end Individuals.h
