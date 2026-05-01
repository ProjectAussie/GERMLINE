// Individuals.cpp: A collection of individuals

#include "Individuals.h"
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
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

static ofstream* openPerDogFile(const string& root,
                                const string& subdir,
                                const string& single_id,
                                const string& chromosome,
                                string& out_path)
{
	string dir = root + "/" + subdir + "/" + single_id;
	try { filesystem::create_directories(dir); }
	catch (const filesystem::filesystem_error& e) {
		throw runtime_error("Cannot create output directory '" + dir + "': " + e.what());
	}
	out_path = dir + "/chr" + chromosome + ".tsv";
	auto* ofs = new ofstream(out_path, ofstream::app);
	if (!ofs->is_open()) {
		delete ofs;
		throw runtime_error("Cannot open per-dog output file: " + out_path);
	}
	return ofs;
}

// Individuals(): default constructor
Individuals::Individuals()
{}

Individuals::~Individuals()
{
	// Safety net only — GERMLINE::mine() should call closeOutputFileHandles()
	// on the clean path so a sort(1) failure surfaces as a runtime error
	// instead of being swallowed during global teardown.
	try { closeOutputFileHandles(); } catch (...) {}
	for(begin();more();next())
		delete pedigree[ iter ];
}

void Individuals::initialize()
{
}

void Individuals::initializeOutputFileHandles(string chromosome)
{
	// Dedup ofstream allocation by single_id: in haploid mode the pedigree
	// holds two Individual objects per dog with the same single_id. Both
	// haplotypes share the underlying file, so they must share one ofstream
	// (and one userspace buffer). Independent buffers writing through O_APPEND
	// fds produced the NF=13 / NF=0 corruption in v1.7 (SCICO-1236).
	for ( iter = 0 ; iter < pedigree.size() ; iter++ ) {
		if ( !pedigree[ iter ]->is_new ) continue;
		const string& sid = pedigree[ iter ]->single_id;

		auto match_it = match_file_by_single_id.find(sid);
		if (match_it == match_file_by_single_id.end()) {
			string path;
			ofstream* ofs = openPerDogFile(individualOutputFolder, "dog_level_match_files", sid, chromosome, path);
			match_file_by_single_id[sid] = ofs;
			match_path_by_single_id[sid] = path;
			match_it = match_file_by_single_id.find(sid);
		}
		pedigree[ iter ]->setIndividualMatchFile(match_it->second);

		auto homoz_it = homoz_file_by_single_id.find(sid);
		if (homoz_it == homoz_file_by_single_id.end()) {
			string path;
			ofstream* ofs = openPerDogFile(individualOutputFolder, "dog_level_homoz_files", sid, chromosome, path);
			homoz_file_by_single_id[sid] = ofs;
			homoz_path_by_single_id[sid] = path;
			homoz_it = homoz_file_by_single_id.find(sid);
		}
		pedigree[ iter ]->setIndividualHomozFile(homoz_it->second);
	}
}

void Individuals::closeOutputFileHandles()
{
	// Idempotent best-effort close. Two callers: an explicit drive-by from
	// GERMLINE::mine (the clean path that surfaces sort failures) and the
	// destructor (the safety net). The destructor cannot throw, so any
	// sortFileInPlace failure is captured here, the rest of the cleanup
	// proceeds, and the first stored exception is rethrown only after every
	// stream is closed and every map is cleared.
	std::exception_ptr stored;

	auto run_sort = [&](const string& path) {
		if (UNSORTED_OUTPUT) return;
		try { sortFileInPlace(path); }
		catch (...) { if (!stored) stored = std::current_exception(); }
	};

	for (auto& [sid, ofs] : match_file_by_single_id) {
		delete ofs;
		run_sort(match_path_by_single_id[sid]);
	}
	match_file_by_single_id.clear();
	match_path_by_single_id.clear();

	for (auto& [sid, ofs] : homoz_file_by_single_id) {
		delete ofs;
		run_sort(homoz_path_by_single_id[sid]);
	}
	homoz_file_by_single_id.clear();
	homoz_path_by_single_id.clear();

	for (auto* ind : pedigree) {
		if (!ind) continue;
		ind->setIndividualMatchFile(nullptr);
		ind->setIndividualHomozFile(nullptr);
	}

	if (stored) std::rethrow_exception(stored);
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
