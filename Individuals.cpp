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

void Individuals::loadOldIndividuals(string f) {
	ifstream s_old(f.c_str());
	if (!s_old) {
		cerr << "WARNING: List of old individuals \"" << f << "\" could not be loaded" << endl;
		return;
	}

	string fam_id, ind_id;

	while(!s_old.eof()) {
		s_old >> fam_id >> ind_id;
                //cout << fam_id << ind_id << endl;
                old_samples.insert(fam_id + " " + ind_id);
	}
}


void Individuals::loadNewIndividuals(string f) {
	ifstream s_new(f.c_str());
	if (!s_new) {
		cerr << "WARNING: List of new individuals \"" << f << "\" could not be loaded" << endl;
		return;
	}

	string fam_id, ind_id;

	while(!s_new.eof()) {
		s_new >> fam_id >> ind_id;
		new_samples.insert(fam_id + " " + ind_id);
	}
}


bool Individuals::hasRestrictions() {
  return !old_samples.empty();
}

bool Individuals::isOld(string indBaseID) {
  //cout << ind.getBaseID() << endl;
  return old_samples.find(indBaseID) != old_samples.end();
}


bool Individuals::isNew(string indBaseID) {
  return new_samples.find(indBaseID) != new_samples.end();
}

// end Individuals.cpp
