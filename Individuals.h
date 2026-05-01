// Individuals.h: A collection of individuals

#ifndef INDIVIDUALS_H
#define INDIVIDUALS_H

#include "BasicDefinitions.h"
#include "Chromosome.h"
#include "Individual.h"
#include <fstream>
#include <set>
#include <ostream>
#include <unordered_map>
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
	// Flush + close all per-dog ofstreams and run the in-place sort once per
	// dog. Safe to call multiple times. Called explicitly from ~Individuals
	// before the pedigree is destroyed.
	void closeOutputFileHandles();
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

	// In haploid mode each dog is loaded as two Individual objects (haplotype
	// .0 and .1) sharing one single_id. Per-dog match/homoz output files MUST
	// be written through a single ofstream so the userspace buffer never
	// interleaves with another buffer pointed at the same file (SCICO-1241).
	// Owned here so each file is closed and sorted exactly once.
	unordered_map<string, ofstream*> match_file_by_single_id;
	unordered_map<string, ofstream*> homoz_file_by_single_id;
	unordered_map<string, string> match_path_by_single_id;
	unordered_map<string, string> homoz_path_by_single_id;
};

#endif

// end Individuals.h
