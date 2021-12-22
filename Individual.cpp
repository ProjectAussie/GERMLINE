// Individual.cpp: An individual with genetic data

#include "Individual.h"
using namespace std;

// Individual(): default constructor
Individual::Individual()
{
	if ( HAPLOID ) {
		if (DEBUG) cout << "Individual() constructor; haploid mode" << endl;
		chromosome = new Chromosome[1];
	} else {
		chromosome = new Chromosome[2];
	}
	numeric_id = 0;
}

Individual::~Individual()
{
	delete[] chromosome;
	delete[] all_matches;
}

void Individual::freeMatches()
{
	for ( size_t iter = 0 ; iter < num_samples ; iter++ )
		if ( all_matches[ iter ] != NULL ) deleteMatch( iter );
}

Match * Individual::getMatch( size_t id )
{
	return all_matches[ id ];
}

void Individual::assertHomozygous()
{
	size_t iter = this->getNumericID();
	Match * m;
	if ( all_matches[ iter ] != NULL )
	{
		// increment this match
		all_matches[ iter ]->end_ms = position_ms;

	} else
	{	
		// this is a new match
		m = new Match();
		if (DEBUG) cout << "new Match() in Individual.cpp::assertHomozygous, assigning start_ms and end_ms to " << position_ms << endl;
		m->end_ms = m->start_ms = position_ms;
		m->node[0] = m->node[1] = this;
		if (DEBUG) cout << "extendBack called in Individual.cpp" << endl;
		m->extendBack();
		all_matches[ iter ] = m;
	}
}

void Individual::assertShares()
{
	Match * m;
	set<Individual*>::iterator cip;

	// try to extend previous matches that did not match currently
	for( size_t iter = 0 ; iter < num_samples ; iter++ )
	{
		if ( all_matches[ iter ] == NULL ) continue;

		m = all_matches[ iter ];
		// Can we increment?
		if ( m->approxEqual() ) m->end_ms = position_ms;
		else deleteMatch( iter );
	}
}

void Individual::clearMatch( size_t id )
{
	all_matches[ id ] = NULL;
}
void Individual::deleteMatch( size_t id )
{
	// try to print it
	// cout << "Writing results for: " << single_id << endl;
	all_matches[ id ]->print( MATCH_FILE );
	delete all_matches[ id ];

	// erase from the list
	clearMatch( id );
}

void Individual::addMatch( size_t id , Match * m)
{
	all_matches[ id ] = m;
}

void Individual::reserveMemory()
{
	all_matches = new Match * [ num_samples ];
	for ( size_t i = 0 ; i < num_samples ; i++ ) all_matches[ i ] = NULL;
}

void Individual::print(ostream& out,long start,long end)
{
	short tot;
	if ( HAPLOID ) tot=1; else tot=2;
	for(int i=0;i<tot;i++)
	{
		out << getID() << '\t';
		chromosome[i].print(out,start,end);
		out << endl;
	}
}

int Individual::numHet()
{
	if ( HAPLOID ) return 0;
	else return int(( chromosome[0].getMarkerSet()->getMarkerBits() ^ chromosome[1].getMarkerSet()->getMarkerBits() ).count());
}

bool Individual::isHeterozygous()
{
	if ( HAPLOID ) return false;
	else return !( chromosome[0].getMarkerSet()->equal( chromosome[1].getMarkerSet() ) );
}

bool Individual::isHeterozygous(int i)
{
	if ( HAPLOID ) return false;
	else return chromosome[0].getMarkerSet()->getMarker(i) != chromosome[1].getMarkerSet()->getMarker(i);
}

void Individual::setOffset(streamoff o)
{
	offset = o;
}

streamoff Individual::getOffset()
{
	return offset;
}

// getID(): accessor for ID
string Individual::getID() const
{
	return ID;
}

string Individual::getBaseID() const
{
	return BaseID;
}

Chromosome * Individual::getAlternateChromosome( Chromosome * c)
{
	if ( HAPLOID ) return &(chromosome[0]);
	else
	{
		if( &(chromosome[0]) == c ) return &(chromosome[1]); else return &(chromosome[0]);
	}
}

// returns a pointer to a chromosome?
Chromosome * Individual::getChromosome(int ct)
{
	if ( HAPLOID ) {
		// AG: if haploid, hardcode "ct" to zero. Not sure why.
		// return memory address of h[0] ?
		return &(chromosome[0]);
	} else {
		return &(chromosome[ct]);
	}
}

unsigned int Individual::getNumericID()
{
	return numeric_id;
}

void Individual::setNumericID( unsigned int id )
{
	numeric_id = id;
}

// setID(): mutator for ID.
void Individual::setID(string id)
{
	ID = id;
}

void Individual::setBaseID(string id)
{
	BaseID = id;
}

void Individual::clearMarkers()
{
	chromosome[0].clear();
	if ( !HAPLOID ) {
		chromosome[1].clear();
	}
}

// addMarkerSet(): adds MarkerSet to a chromosome
void Individual::addMarkerSet(int ct, MarkerSet * marker_set)
{
	if (DEBUG) cout << "Individual::addMarkerSet" << endl;
	if ( HAPLOID ) {
		ct = 0;
	}
	chromosome[ct].addMarkerSet(marker_set);
}

// operator<<(): overloaded stream insertion operator
ostream& operator<<(ostream &fout, Individual& ind)
{
	fout << ind.getID() << endl;
	fout << ind.getChromosome(0) << endl;
	fout << ind.getChromosome(1) << endl;
	return fout;
}

void Individual::setIndividualMatchFile(string chromosome)
{
	string ext = ".tsv";
	string dir = ALL_SAMPLES.individualOutputFolder + "/dog_level_match_files/" + single_id;
	experimental::filesystem::path _dir(dir);
	if ( !experimental::filesystem::exists(_dir) ) {
		experimental::filesystem::create_directories(_dir);
	}
	string fileHandleName = dir + "/chr" + chromosome + ext;
	// cout << fileHandleName << endl;
	individualMatchFile = new ofstream(fileHandleName, ofstream::app);
}

void Individual::setIndividualHomozFile(string chromosome)
{
	string ext = ".tsv";
	string dir = ALL_SAMPLES.individualOutputFolder + "/dog_level_homoz_files/" + single_id;
	experimental::filesystem::path _dir(dir);
	if ( !experimental::filesystem::exists(_dir) ) {
		experimental::filesystem::create_directories(_dir);
	}
	string fileHandleName = dir + "/chr" + chromosome + ext;
	// cout << fileHandleName << endl;
	individualHomozFile = new ofstream(fileHandleName, ofstream::app);
}

ofstream* Individual::getIndividualMatchFile()
{
	return individualMatchFile;
}

ofstream* Individual::getIndividualHomozFile()
{
	return individualHomozFile;
}


// end Individual.cpp
