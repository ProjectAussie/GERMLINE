// Chromosome.cpp.haplotyped markers for a chromosome

#include "Chromosome.h"
#include <iostream>
using namespace std;


// Chromosome(): default constructor
Chromosome::Chromosome()
{}

MarkerSet * Chromosome::getMarkerSet()
{
	return &chromosome[position_ms];
}

MarkerSet * Chromosome::getMarkerSet(unsigned int pos)
{
	return &chromosome[pos];
}

void Chromosome::clear()
{
	chromosome.clear();
}

void Chromosome::addMarkerSet(MarkerSet * marker_set)
{
	if (DEBUG) cout << "Chromosome.addMarkerSet called" << endl;
	// Pre-reserve capacity on first insert of a load so the vector doesn't
	// reallocate (and move every prior MarkerSet) as marker sets stream in.
	// num_sets is the marker-set count for the current chromosome, set in
	// GERMLINE::mine before buildMatches runs. Vector::clear() keeps
	// capacity, so subsequent chromosomes only reserve again if they are
	// larger than any previously seen one.
	if (chromosome.empty() && num_sets > chromosome.capacity())
		chromosome.reserve(num_sets);
	chromosome.push_back(std::move(*marker_set));
	delete marker_set;
}

void Chromosome::print_snps(ostream& out, unsigned int start, unsigned int end)
{
	unsigned int p_ms = position_ms;

	unsigned int ms_start = start / MARKER_SET_SIZE;
	unsigned int ms_end = end / MARKER_SET_SIZE;
	if( start % MARKER_SET_SIZE != 0 ) { position_ms = ms_start; chromosome[ms_start++].print(out,start % MARKER_SET_SIZE,MARKER_SET_SIZE); out << ' '; }
	print(out,ms_start,ms_end);
	if( end % MARKER_SET_SIZE != 0 ) { out << ' '; chromosome[ms_end].print(out,0,end % MARKER_SET_SIZE); }

	position_ms = p_ms;
}

void Chromosome::print(ostream& out,unsigned int start,unsigned int end)
{
	for(position_ms=start;position_ms<end;position_ms++)
	{
		if( position_ms > start ) out << ' ';
		chromosome[position_ms].print(out);
	}
}

ostream& operator<<(ostream &fout, Chromosome& c)
{
	fout << c.getMarkerSet();
	return fout;
}

// end Chromosome.cpp
