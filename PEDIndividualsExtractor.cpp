// KOSIndividualsExtractor.cpp: Manages input from Plink files

#include "PEDIndividualsExtractor.h"
#include <cctype>
#include <iostream>
using namespace std;

// PEDIndividualsExtractor(): default constructor
PEDIndividualsExtractor::PEDIndividualsExtractor()
{}

void PEDIndividualsExtractor::stripWhitespace()
{
	if (stream.is_open())
	{
		char c;
		while((c=stream.peek())!=EOF && isspace(c))
		stream.get();
	}
	else
	{
		cerr << "WARNING:PolymorphicIndividualsExtractor::stripWhiteSpace():stream is not open" << endl;
		valid_flag = false;
	}
}

void PEDIndividualsExtractor::loadInput()
{
	individualsP = &ALL_SAMPLES;
	ALL_SNPS.processMAPFile();
	ALL_SNPS.beginChromosome();
	numberOfMarkers = ALL_SNPS.size();

	while (stream)  // Adam G. changed this to get germline 1.5.1 from columbia's website to run on ubuntu
	{
		getIndividuals();
		stream.seekg(numberOfMarkers*4 + 1,ios::cur);
	}
	
	individualsP->initialize();
	stream.clear();
}

// getInput(): gets individuals from .ped file
void PEDIndividualsExtractor::getInput()
{
	cout << "Please enter the MAP file name" << endl;
	cin >> map_file;
	cout << "Please enter the PED file name" << endl;
	cin >> ped_file;
	
	if ( !ALL_SNPS.setFile( map_file ) )
	{
		cerr << "WARNING:PEDIndividualsExtractor::getInput():cannot open map file" << endl;
		valid_flag = false;
		return;
	}

	stream.open( ped_file.c_str() );
	if ( !stream )
	{
		cerr << "WARNING:PEDIndividualsExtractor::getInput():cannot open ped file" << endl;
		valid_flag = false;
		return;
	}
}

// getIndividuals(): gets the next nuclear family from stream
void PEDIndividualsExtractor::getIndividuals()
{
	string discard, ID, famID;
	stream >> famID >> ID >> discard >> discard >> discard >> discard;
	if(!stream.good()) return;
	if ( HAPLOID )
	{
		Individual * new_ind[2];
		new_ind[0] = new Individual();
		new_ind[1] = new Individual();
		new_ind[0]->setOffset( stream.tellg() );
		new_ind[1]->setOffset( stream.tellg() );
		new_ind[0]->setID(famID + " " + ID + ".0" );
		new_ind[1]->setID(famID + " " + ID + ".1" );
		
		individualsP->addIndividual( new_ind[0] );
		individualsP->addIndividual( new_ind[1] );
	} else
	{
		Individual * new_ind = new Individual;
		new_ind->setOffset(stream.tellg());
		new_ind->setID(famID + " " + ID);
		individualsP->addIndividual(new_ind);
	}
}

void PEDIndividualsExtractor::getCompleteMarkerSet(Individual * p)
{
	if (DEBUG) cout << "PEDIndividualsExtractor::getCompleteMarkerSet" << endl;
	stream.seekg(p->getOffset() + 4*ALL_SNPS.getROIStart().getMarkerNumber() + 4*position_ms*MARKER_SET_SIZE + 1);
	MarkerSet * marker_sets[2];
	marker_sets[0] = new MarkerSet();
	marker_sets[1] = new MarkerSet();

	readMarkerSet( marker_sets );

	p->addMarkerSet(UNTRANS, marker_sets[0]);
	p->addMarkerSet(TRANS, marker_sets[1]);
}

void PEDIndividualsExtractor::readMarkerSet( MarkerSet ** marker_set )
{
	unsigned int maxsize = ALL_SNPS.currentSize();
	if (DEBUG) cout << "PEDIndividualsExtractor::readMarkerSet; max_size: " << maxsize << " MARKER_SET_SIZE: " << MARKER_SET_SIZE << endl;

	for (int position = 0; position < MARKER_SET_SIZE; position++)
	{
		if (position_ms * MARKER_SET_SIZE + position >= maxsize) {
			break;
		}
		int overall_position = position_ms * MARKER_SET_SIZE + position;
		if (DEBUG) cout << "position " << position << "; overall position: " << overall_position << endl;
		for (int allele = 0; allele < 2; allele++) { // AG: what is "al"?
			stripWhitespace();
			char marker = stream.peek();
			if (DEBUG) cout << "allele " << allele << ": " << marker << endl;
			int allele_as_binary = ALL_SNPS.mapNucleotideToBinary(marker, overall_position);
			if ( allele_as_binary == 1 ) {
				marker_set[allele]->set(position, true); // AG: are we collapsing the sets here?
			}
			stream.get();
		}
	}
}

void PEDIndividualsExtractor::getCompleteMarkerSet(Individual * p0 , Individual * p1 )
{
	if (DEBUG) cout << "PEDIndividualsExtractor::getCompleteMarkerSet with two params" << endl;
	stream.seekg(p0->getOffset() + 4*ALL_SNPS.getROIStart().getMarkerNumber() + 4*position_ms*MARKER_SET_SIZE + 1);
	MarkerSet * marker_sets[2];
	marker_sets[0] = new MarkerSet();
	marker_sets[1] = new MarkerSet();

	if (DEBUG) cout << "readMarkerSet(marker_sets)" << endl;
	readMarkerSet(marker_sets);

	if (DEBUG) cout << "addMarkerSet marker_sets[0]" << endl;
	p0->addMarkerSet(TRANS, marker_sets[0]);

	if (DEBUG) cout << "addMarkerSet marker_sets[1]" << endl;
	p1->addMarkerSet(TRANS, marker_sets[1]);
}

// end PEDIndividualsExtractor.cpp

