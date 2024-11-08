#include "GERMLINE.h"
#include <string.h>

const std::string VERSION = "1.6.0-embark";
double MIN_MATCH_LEN = 3;
int MARKER_SET_SIZE = 128;
bool PRINT_MATCH_HAPS = false;
bool ROI = false;
bool HAP_EXT = false;
bool WIN_EXT = false;
bool ALLOW_HOM = false;
bool HOM_ONLY = false;
bool HAPLOID = false;
bool SILENT = false;
bool DEBUG = false;
bool BINARY_OUT = false;
int MAX_ERR_HOM = 4;
int MAX_ERR_HET = 1;

// main(): runs GERMLINE
int main(int argc, char* argv[])
{
	// parse arguments
	string rs_range[2], map, samples_to_compare_to, new_samples, chromosome;
	samples_to_compare_to = "";
	new_samples = "";
	map = "";
	rs_range[0] = "";
	rs_range[1] = "";
	string params = argv[0];
	ALL_SAMPLES.useEmbarkRFGermlineOutput = false;
	ALL_SAMPLES.chromosome = "";
	ALL_SAMPLES.individualOutputFolder = "";

	bool bad_param = false;
	for(int i=1; i<argc; i++) {
		params += " " + string(argv[i]);
		if (strncmp(argv[i], "-min_m", strlen("-min_m")) == 0 && i < argc-1) {
			MIN_MATCH_LEN = atof(argv[++i]);
		}
		else if (strncmp(argv[i], "-err_hom", strlen("-max_err")) == 0 && i < argc-1) {
			MAX_ERR_HOM = atoi(argv[++i]);
		}
		else if (strncmp(argv[i], "-err_het", strlen("-max_err")) == 0 && i < argc-1) {
			MAX_ERR_HET = atoi(argv[++i]);
		}
		else if (strncmp(argv[i], "-from_snp", strlen("-from_snp")) == 0 && i < argc-1 ) {
			rs_range[0] = argv[++i];
		}
		else if (strncmp(argv[i], "-to_snp", strlen("-to_snp")) == 0 && i < argc-1 ) {
			rs_range[1] = argv[++i];
		}
		else if (strncmp(argv[i], "-print", strlen("-print")) == 0 ) {
			PRINT_MATCH_HAPS = true;
		}
		else if (strncmp(argv[i], "-silent", strlen("-silent")) == 0 ) {
			SILENT = true;
		}
		else if (strncmp(argv[i], "-debug", strlen("-debug")) == 0 ) {
			DEBUG = true;
		}
		else if (strncmp(argv[i], "-map", strlen("-map")) == 0 && i < argc-1) {
			map = argv[++i];
		}
		else if (strncmp(argv[i], "-bits", strlen("-bits")) == 0 && i < argc-1) {
			MARKER_SET_SIZE = atoi(argv[++i]);
		}
		else if (strncmp(argv[i], "-homoz-only", strlen("-homoz-only")) == 0 ) {
			ALLOW_HOM = true;
			HOM_ONLY = true;
		}
		else if (strncmp(argv[i], "-homoz", strlen("-homoz")) == 0 ) {
			ALLOW_HOM = true;
		}
		else if (strncmp(argv[i], "-bin_out", strlen("-bin_out")) == 0 ) {
			BINARY_OUT = true;
		}
		else if (strncmp(argv[i], "-haploid", strlen("-haploid")) == 0 ) {
			HAPLOID = true;
			HAP_EXT = true;
		}
		else if( strncmp(argv[i], "-h_extend", strlen("-h_extend")) == 0 ) {
			HAP_EXT = true;
		}
		else if (strncmp(argv[i], "-w_extend", strlen("-w_extend")) == 0 ) {
			WIN_EXT = true;
		}
		else if( strncmp(argv[i], "-version", strlen("-version")) == 0 ) {
			cout << "version" << " " << VERSION << endl;
			return 0;
		}
		else if (strncmp(argv[i], "-samples_to_compare_to", strlen("-samples_to_compare_to")) == 0 && i < argc-1) {
			samples_to_compare_to = argv[++i];
		}
		else if (strncmp(argv[i], "-new_samples", strlen("-new_samples")) == 0 && i < argc-1) {
			new_samples = argv[++i];
		}
		else if (strncmp(argv[i], "-chromosome", strlen("-chromosome")) == 0 && i < argc-1) {
			chromosome = argv[++i];
			chromosome.insert(chromosome.begin(), 2 - chromosome.length(), '0');
			ALL_SAMPLES.chromosome = chromosome;
			cout << "Generating outputs for chromosome: " << chromosome << endl;
		}
		else if (strncmp(argv[i], "-individual_outputs", strlen("-individual_outputs")) == 0) {
			ALL_SAMPLES.individualOutputFolder = argv[++i];
			ALL_SAMPLES.useEmbarkRFGermlineOutput = true;
			cout << "Generating individual outputs" << endl;
		}
		else {
			cout << argv[i] << endl;
			bad_param = true;
		}
	}

	if(MIN_MATCH_LEN < 0)
	{
		cerr << "-min_m must be non-negative" << endl << endl;
		bad_param = true;
	} else if(MAX_ERR_HOM < 0 || MAX_ERR_HET < 0 )
	{
		cerr << "-err_hom,-err_het must be non-negative" << endl << endl;
		bad_param = true;
	}

	if(bad_param)
	{
		cerr << "usage: " << argv[0] << "<flags (optional)>" << endl
		<< "flags:" << endl
		<< '\t' << "-version" << '\t' << "Print version information." << endl
		<< '\t' << "-silent" << '\t' << "Suppress all output except for warnings and prompts." << endl
		<< '\t' << "-bin_out" << '\t' << "Output in binary format to save space." << endl
		<< '\t' << "-min_m" << '\t' << "Minimum length for match to be used for imputation (in cM or MB)." << endl
		<< '\t' << "-err_hom" << '\t' << "Maximum number of mismatching homozygous markers (per slice)." << endl
		<< '\t' << "-err_het" << '\t' << "Maximum number of mismatching heterozygous markers (per slice)." << endl
		<< '\t' << "-from_snp" << '\t' << "Start SNP (rsID)." << endl
		<< '\t' << "-to_snp" << '\t' << "End SNP (rsID)." << endl
		<< '\t' << "-haps" << '\t' << "Print the resolved haplotypes in a seperate HAPS file." << endl
		<< '\t' << "-map" << '\t' << "Genetic distance map." << endl
		<< '\t' << "-bits" << '\t' << "Slice size." << endl
		<< '\t' << "-homoz" << '\t' << "Allow self matches (homozygosity)" << endl
		<< '\t' << "-homoz-only" << '\t' << "Look for autozygous/homozygous segments only, does not detect IBD" << endl
		<< '\t' << "-haploid" << '\t' << "Treat input individual as two fully phased chromosomes with no recombination\n\t\toutput IDs with 0/1 suffix for chromosome destinction" << endl
		<< '\t' << "-h_extend" << '\t' << "Extend from seeds if *haplotypes* match" << endl
		<< '\t' << "-w_extend" << '\t' << "Extend, one marker at a time, beyong the boundaries of a found match" << endl
		<< '\t' << "-samples_to_compare_to" << '\t' << "List of individuals to cross compare." << endl
		<< '\t' << "-new_samples" << '\t' << "List of new individuals." << endl
		<< '\t' << "-chromosome" << '\t' << "Chromosome number for individual match/homoz file outputs." << endl
		<< '\t' << "-individual_outputs" << '\t' << "Top level directory to generate individual match/homoz file outputs to" << endl;
		return 1;
	}

	if( rs_range[0] != "" && rs_range[1] != "" )
	{
		ROI = true;
		ALL_SNPS.setROI(rs_range);
	}

	if(map != "") {
		ALL_SNPS.loadGeneticDistanceMap(map);
	}
	if (samples_to_compare_to != "") {
		ALL_SAMPLES.loadOldIndividuals(samples_to_compare_to);
	}
	if (new_samples != "") {
		ALL_SAMPLES.loadNewIndividuals(new_samples);
	}

	GERMLINE germline;
    germline.mine( params );
    return 0;
}


// end GERMLINE_0001.cpp
