// Individual.cpp: An individual with genetic data

#include "Individual.h"
#include <cstdlib>
#include <stdexcept>
#include <string>
using namespace std;

static void sortFileInPlace(const string& path)
{
	if (path.empty()) return;
	if (path.find('\'') != string::npos)
		throw runtime_error("Cannot sort output file (path contains single quote): " + path);
	string cmd = "LC_ALL=C sort -S 128M -o '" + path + "' '" + path + "'";
	if (std::system(cmd.c_str()) != 0)
		throw runtime_error("sort(1) failed on " + path);
}

Individual::Individual()
{
	if ( HAPLOID ) {
		if (DEBUG) cout << "Individual() constructor; haploid mode" << endl;
		chromosome = new Chromosome[1];
	} else {
		chromosome = new Chromosome[2];
	}
	numeric_id = 0;
	is_new = false;
	is_old = false;
	individualMatchFile = nullptr;
	individualHomozFile = nullptr;
}

Individual::~Individual()
{
	delete[] chromosome;
	for ( auto& [id, m] : all_matches ) delete m;
	// Per-individual TSVs are written in hash-map order (same as MATCH_FILE).
	// Sort them in place after the stream flushes so each individual's match
	// / homoz file is deterministic for downstream consumers. Skipped when
	// -unsorted_output is set, matching the global .match behavior.
	delete individualMatchFile;
	if ( !UNSORTED_OUTPUT ) sortFileInPlace(individualMatchFilePath);
	delete individualHomozFile;
	if ( !UNSORTED_OUTPUT ) sortFileInPlace(individualHomozFilePath);
}

void Individual::freeMatches()
{
	for ( auto& [id, m] : all_matches ) { m->print( MATCH_FILE ); delete m; }
	all_matches.clear();
}

Match * Individual::getMatch( size_t id )
{
	auto it = all_matches.find( (unsigned int)id );
	if ( it != all_matches.end() ) return it->second;
	return nullptr;
}

void Individual::assertHomozygous()
{
	unsigned int iter = this->getNumericID();
	Match * m;
	auto it = all_matches.find( iter );
	if ( it != all_matches.end() )
	{
		it->second->end_ms = position_ms;

	} else
	{
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
	auto it = all_matches.begin();
	while ( it != all_matches.end() )
	{
		if ( it->second->approxEqual() ) { it->second->end_ms = position_ms; ++it; }
		else { it->second->print( MATCH_FILE ); delete it->second; it = all_matches.erase( it ); }
	}
}

void Individual::deleteMatch( size_t id )
{
	auto it = all_matches.find( (unsigned int)id );
	if ( it == all_matches.end() ) { if (DEBUG) cerr << "deleteMatch: id " << id << " not found" << endl; return; }
	it->second->print( MATCH_FILE );
	delete it->second;
	all_matches.erase( it );
}

void Individual::addMatch( size_t id , Match * m)
{
	auto [ it, inserted ] = all_matches.emplace( (unsigned int)id, m );
	if ( !inserted ) { delete it->second; it->second = m; }
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

// Userspace buffer size for per-individual ofstreams. The default
// libstdc++ filebuf is ~8 KB; at 14 GB of per-chromosome output
// distributed across potentially thousands of files, that's hundreds
// of thousands of write() syscalls. 32 KB keeps syscall count and RSS
// growth (per_stream * n_individuals * 2 files) in balance.
static constexpr size_t INDIV_FILE_BUF_SIZE = 32 * 1024;

static void attachLargeBuffer(ofstream* ofs, vector<char>& buffer)
{
	buffer.resize(INDIV_FILE_BUF_SIZE);
	ofs->rdbuf()->pubsetbuf(buffer.data(), buffer.size());
}

void Individual::setIndividualMatchFile(string chromosome)
{
	string ext = ".tsv";
	string dir = ALL_SAMPLES.individualOutputFolder + "/dog_level_match_files/" + single_id;
	try { filesystem::create_directories(dir); }
	catch (const filesystem::filesystem_error& e) { throw runtime_error("Cannot create output directory '" + dir + "': " + e.what()); }
	individualMatchFilePath = dir + "/chr" + chromosome + ext;
	individualMatchFile = new ofstream();
	attachLargeBuffer(individualMatchFile, individualMatchFileBuffer);
	individualMatchFile->open(individualMatchFilePath, ofstream::app);
}

void Individual::setIndividualHomozFile(string chromosome)
{
	string ext = ".tsv";
	string dir = ALL_SAMPLES.individualOutputFolder + "/dog_level_homoz_files/" + single_id;
	try { filesystem::create_directories(dir); }
	catch (const filesystem::filesystem_error& e) { throw runtime_error("Cannot create output directory '" + dir + "': " + e.what()); }
	individualHomozFilePath = dir + "/chr" + chromosome + ext;
	individualHomozFile = new ofstream();
	attachLargeBuffer(individualHomozFile, individualHomozFileBuffer);
	individualHomozFile->open(individualHomozFilePath, ofstream::app);
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
