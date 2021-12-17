// Individuals.cpp: A collection of individuals

#include "Individuals.h"
#include <iostream>
using namespace std;


// Individuals(): default constructor
Individuals::Individuals()
{}

Individuals::~Individuals()
{
	for(begin();more();next())
		delete pedigree[ iter ];
}

void Individuals::initialize()
{
	for ( iter = 0 ; iter < pedigree.size() ; iter++ ) pedigree[ iter ]->reserveMemory();
}

void Individuals::initializeOutputFileHandles()
{
	for ( iter = 0 ; iter < pedigree.size() ; iter++ ) {
		if ( pedigree[ iter ]->is_new ) {
			pedigree[ iter ]->setIndividualMatchFile();
			pedigree[ iter ]->setIndividualHomozFile();
		}
	}
}

void Individuals::freeMatches()
{
	for(begin();more();next()) pedigree[ iter ]->freeMatches();
}

void Individuals::freeMarkers()
{
	for(begin();more();next()) { pedigree[ iter ]->clearMarkers(); }
}

void Individuals::print( ostream& out )
{
	for(begin();more();next())
		out << pedigree[ iter ]->getID() << endl;
}

void Individuals::begin()
{
	iter = 0;
}

bool Individuals::more()
{
	return iter < pedigree.size();
}

Individual * Individuals::next()
{
	return pedigree[ iter++ ];
}

// addIndividual(): adds an Individual object
void Individuals::addIndividual(Individual * ind)
{
	pedigree.push_back(ind);
	ind->setNumericID( (unsigned int) num_samples++ );
}


void _loadIndividuals(string f, set<string> *samples) {
	ifstream s(f.c_str());
	if (!s) {
		cerr << "WARNING: List of individuals \"" << f << "\" could not be loaded" << endl;
		return;
	}

	string fam_id, ind_id;

	while(!s.eof()) {
		s >> fam_id >> ind_id;
		samples->insert(fam_id + " " + ind_id);
	}
}

void Individuals::loadOldIndividuals(string f) {
	_loadIndividuals(f, &samples_to_compare_to);
}


void Individuals::loadNewIndividuals(string f) {
	_loadIndividuals(f, &new_samples);
}


bool Individuals::hasRestrictions() {
  return !samples_to_compare_to.empty();
}

bool Individuals::isOld(string indBaseID) {
  return samples_to_compare_to.find(indBaseID) != samples_to_compare_to.end();
}


bool Individuals::isNew(string indBaseID) {
  return new_samples.find(indBaseID) != new_samples.end();
}

// end Individuals.cpp
