#include "Match.h"

void Match::extendBack()
{
	// save old position values
	unsigned int SAVE_pms = position_ms;
	cout << "extendBack called. position_ms: " << position_ms << endl;
	
	// iterate backwards through genome
	while(position_ms > 0)
	{
		position_ms--;
		if( !approxEqual() )
		{
			position_ms++;
			break;
		}
	}
	cout << "assigning start_ms value: " << position_ms << endl;
	start_ms = position_ms;
	// restore saved values
	position_ms = SAVE_pms;
}

bool Match::approxEqual()
{
	// homozygosity check
	if ( node[0] == node[1] )
	{
		if ( ALLOW_HOM )
		{
			if ( (int) ( node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).count() 
				 <= ( MAX_ERR_HOM + MAX_ERR_HET ) ) return true; else return false;
		}
		else
		{
			return false;
		}
	}
	else
	{
		// 1. Haplotype extension
		if ( HAPLOID )
		{
				if ( (int)(node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()).count() <= MAX_ERR_HOM ) return true;
		} else
		{
			for ( int a = 0 ; a < 2 ; a++ ) {
				for ( int b = 0 ; b < 2 ; b++ ) { 
					if ( (int)(node[0]->getChromosome( a )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( b )->getMarkerSet()->getMarkerBits()).count() <= MAX_ERR_HOM )
					{
						return true;
					}
				}
			}
		}

		if ( HAPLOID || HAP_EXT ) return false;

		// 2. Genotype extension
		// identify common homozygous SNPs
		boost::dynamic_bitset<> mask
			= ( node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).flip()
			& ( node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet()->getMarkerBits() ).flip();

		// assert that homozygous SNPs are identical
		if ( (int) ((node[0]->getChromosome( 0 )->getMarkerSet()->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet()->getMarkerBits()) & mask).count() <= MAX_ERR_HET )
		{
			return true;
		}
		else return false;
	}
}

int Match::scanLeft( unsigned int marker_set_number )
{
	bool found_mismatch = false;
	int marker = MARKER_SET_SIZE - 1; // e.g. 40
	cout << "initial marker in scanLeft: " << marker << endl;
	cout << "marker_set_number in scanLeft: " << marker_set_number << endl;

	if (HAPLOID) {
		cout << "SCANNING LEFT IN HAPLOID MODE BUT NOT HAP_EXT" << endl;
		cout << "for marker from 41 to 0 (or error)" << endl;

		// when HAPLOID is true, getChromosome resets its parameter to zero, so the input argument doesn't matter
		boost::dynamic_bitset<>& node_zero_bits = node[0]->getChromosome(999)->getMarkerSet(marker_set_number)->getMarkerBits();
		boost::dynamic_bitset<>& node_one_bits = node[1]->getChromosome(999)->getMarkerSet(marker_set_number)->getMarkerBits();

		cout << "node 0 markerBits: " << node_zero_bits << endl;
		cout << "node 1 markerBits: " << node_one_bits << endl;

		for ( marker = MARKER_SET_SIZE - 1 ; marker >= 0 && !found_mismatch; marker-- ) {
			cout << "for-loop; marker (updated with brackets): " << marker << endl;

			// AG: scanLeft goes all the way to the end here, rather than stopping where the match tract goes het
			int node_zero_marker = node_zero_bits[marker];
			int node_one_marker = node_one_bits[marker];

			cout << "node 0 markerBits[" << marker << "]: " << node_zero_marker << endl;
			cout << "node 1 markerBits[" << marker << "]: " <<  node_one_marker << endl;

			if ( node_zero_marker != node_one_marker ) {
				cout << "setting error to true at marker " << marker << endl;
				found_mismatch = true;
			}
		}
	} else if (HAP_EXT)
	{
		cout << "SCANNING LEFT IN HAP_EXT" << endl;
		int cur_marker;
		// AG: looks like we run a function here over 0,0, 0,1, 1,0, 1,1
		for ( int a = 0 ; a < 2 ; a++ ) {
			for ( int b = 0 ; b < 2 ; b++ ) { 
				found_mismatch = false;
				for ( cur_marker = MARKER_SET_SIZE - 1; cur_marker >= 0 && !found_mismatch; cur_marker-- )
				{
					if ( node[0]->getChromosome( a )->getMarkerSet()->getMarkerBits()[cur_marker] != node[1]->getChromosome( b )->getMarkerSet()->getMarkerBits()[cur_marker] )
						found_mismatch = true; // err may actually mean end-of-homozygous-tract?
				}
				// So looks like we start at 40 and advance cur_marker (reducing it?), reassigning marker along the way, and returning it
				// in bits=41 world, we get down to marker=24 and then return it
				if ( cur_marker < marker ) marker = cur_marker;
			}
		}
	} else
	{
		cout << "SCANNING LEFT IN ELSE CLAUSE" << endl;
		cout << "marker_set: " << marker_set_number << endl;

		// AG: big question here: what are nodes? what are chromosomes?

		boost::dynamic_bitset<>& node_zero_chromosome_zero_bits = node[0]->getChromosome(0)->getMarkerSet(marker_set_number)->getMarkerBits();
		boost::dynamic_bitset<>& node_zero_chromosome_one_bits = node[0]->getChromosome(1)->getMarkerSet(marker_set_number)->getMarkerBits();
		boost::dynamic_bitset<>& node_one_chromosome_zero_bits = node[1]->getChromosome(0)->getMarkerSet(marker_set_number)->getMarkerBits();
		boost::dynamic_bitset<>& node_one_chromosome_one_bits = node[1]->getChromosome(1)->getMarkerSet(marker_set_number)->getMarkerBits();

		cout << "node_zero_chromosome_zero_bits: " << node_zero_chromosome_zero_bits << endl;
		cout << "node_zero_chromosome_one_bits: " << node_zero_chromosome_one_bits << endl;
		cout << "node_one_chromosome_zero_bits: " << node_one_chromosome_zero_bits << endl;
		cout << "node_one_chromosome_one_bits: " << node_one_chromosome_one_bits << endl;

		// AG: ^ is XOR. .flip() flips ones and zeros in a bit-set
		boost::dynamic_bitset<> node_zero_is_heterozygous = (node_zero_chromosome_zero_bits ^ node_zero_chromosome_one_bits);
		boost::dynamic_bitset<> node_zero_is_homozygous = (node_zero_chromosome_zero_bits ^ node_zero_chromosome_one_bits).flip();
		boost::dynamic_bitset<> node_one_is_homozygous = (node_one_chromosome_zero_bits ^ node_one_chromosome_one_bits).flip();
		boost::dynamic_bitset<> chromosome_zero_is_heterozygous = (node_zero_chromosome_zero_bits ^ node_one_chromosome_zero_bits);
		boost::dynamic_bitset<> both_nodes_are_homozygous_and_chromosome_zero_is_heterozygous = node_zero_is_homozygous & node_one_is_homozygous & chromosome_zero_is_heterozygous;

		cout << "node_zero_is_heterozygous: " << node_zero_is_heterozygous << endl;
		cout << "node_zero_is_homozygous: " << node_zero_is_homozygous << endl;
		cout << "node_one_is_homozygous: " << node_one_is_homozygous << endl;
		cout << "chromosome_zero_is_heterozygous: " << chromosome_zero_is_heterozygous << endl;
		cout << "both_nodes_are_homozygous_and_chromosome_zero_is_heterozygous: " << both_nodes_are_homozygous_and_chromosome_zero_is_heterozygous << endl;

		boost::dynamic_bitset<> final_mask = both_nodes_are_homozygous_and_chromosome_zero_is_heterozygous;

		if (HOM_ONLY) {
			cout << "HOM_ONLY; CUSTOMIZING MASK TO DETECT HET MARKER WITHIN NODE ZERO" << endl;
			final_mask = node_zero_is_heterozygous;
		}

		cout << "final mask: " << final_mask << endl;

		for (marker = MARKER_SET_SIZE - 1; marker >= 0 && !found_mismatch; marker--) {
			cout << "testing marker " << marker << endl;
			// AG: I'm not sure why this works. Seems like we should be iterating backwards
			// through the mask array, which has the heterozygous SNP (snp 5) 5 places from the end
			if (final_mask[marker]) {
				cout << "setting err to true at marker " << marker << endl;
				found_mismatch = true;
			}
		}
	}

	// AG hypothesis: marker gets set to -1 by the for-loop marker--, then the loop quits
	// so we return -1 as marker
	// cout << "returning marker: " << marker << endl;	
	return marker;
}

int Match::scanRight( unsigned int marker_set_number )
{
	bool found_mismatch = false;
	int marker = 0;

	if (HAPLOID) {
		cout << "SCANNING RIGHT IN HAPLOID MODE" << endl;
		for (marker = 0; marker < MARKER_SET_SIZE && !found_mismatch; marker++) {
			if (node[0]->getChromosome(0)->getMarkerSet()->getMarkerBits()[marker] != node[1]->getChromosome(0)->getMarkerSet()->getMarkerBits()[marker] ) {
				found_mismatch = true;
			}
		}
	} else if (HAP_EXT)
	{
		cout << "SCANNING RIGHT IN HAP_EXT MODE" << endl;
		int cur_marker;
		for (int a = 0; a < 2; a++) {
			for (int b = 0; b < 2; b++) { 
				found_mismatch = false;
				for (cur_marker = 0; cur_marker < MARKER_SET_SIZE && !found_mismatch; cur_marker++) {
					if (node[0]->getChromosome(a)->getMarkerSet()->getMarkerBits()[cur_marker] != node[1]->getChromosome(b)->getMarkerSet()->getMarkerBits()[cur_marker]) {
						found_mismatch = true;
					}
				}
				if (cur_marker > marker) marker = cur_marker;
			}
		}
	} else
	{
		cout << "ELSE-CLAUSE OF SCAN-RIGHT" << endl;
		boost::dynamic_bitset<> node_zero_chromosome_zero_bits = node[0]->getChromosome(0)->getMarkerSet(marker_set_number)->getMarkerBits();
		boost::dynamic_bitset<> node_zero_chromosome_one_bits = node[0]->getChromosome(1)->getMarkerSet(marker_set_number)->getMarkerBits();
		boost::dynamic_bitset<> node_one_chromosome_zero_bits = node[1]->getChromosome(0)->getMarkerSet(marker_set_number)->getMarkerBits();
		boost::dynamic_bitset<> node_one_chromosome_one_bits = node[1]->getChromosome(1)->getMarkerSet(marker_set_number)->getMarkerBits();

		boost::dynamic_bitset<> node_zero_is_heterozygous = (node_zero_chromosome_zero_bits ^ node_zero_chromosome_one_bits);
		boost::dynamic_bitset<> node_zero_is_homozygous = (node_zero_chromosome_zero_bits ^ node_zero_chromosome_one_bits).flip();
		boost::dynamic_bitset<> node_one_is_heterozygous = (node_one_chromosome_zero_bits ^ node_one_chromosome_one_bits);
		boost::dynamic_bitset<> node_one_is_homozygous = (node_one_chromosome_zero_bits ^ node_one_chromosome_one_bits).flip();
		boost::dynamic_bitset<> chromosome_zero_is_heterozygous = (node_zero_chromosome_zero_bits ^ node_one_chromosome_zero_bits);

		boost::dynamic_bitset<> mask = node_zero_is_homozygous & node_one_is_homozygous & chromosome_zero_is_heterozygous;

		if (HOM_ONLY) {
			cout << "HOM_ONLY; CUSTOMIZING MASK TO DETECT HET MARKER WITHIN NODE ZERO" << endl;
			mask = node_zero_is_heterozygous;
		}

		cout << "mask: " << mask << endl;
		for (marker = 0; marker < MARKER_SET_SIZE && !found_mismatch; marker++) {
			if (mask[marker]) {
				found_mismatch = true;
			}
		}
	}
	return marker;
}

int Match::diff( unsigned int ms )
{
	boost::dynamic_bitset<> mask
		= ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[0]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip()
		& ( node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).flip();
	mask = ( node[0]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[1]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits()) & mask;

	return int(mask.count());
}

bool Match::isHom( int n , unsigned int ms )
{
	return (int) ( node[n]->getChromosome( 0 )->getMarkerSet(ms)->getMarkerBits() ^ node[n]->getChromosome( 1 )->getMarkerSet(ms)->getMarkerBits() ).count() <= ( MAX_ERR_HOM + MAX_ERR_HET );
}

void Match::print( ostream& fout )
{
	// extend this match from both ends
	cout << "Match::print" << endl;
	cout << "WIN_EXT: " << WIN_EXT << endl;
	cout << "start_ms:" << endl;
	cout << start_ms << endl;
	cout << "MARKER_SET_SIZE:" << endl;
	cout << MARKER_SET_SIZE << endl;
	unsigned int snp_start = ALL_SNPS.getROIStart().getMarkerNumber() + start_ms * MARKER_SET_SIZE;
	unsigned int snp_end = ALL_SNPS.getROIStart().getMarkerNumber() + ( end_ms + 1 ) * MARKER_SET_SIZE - 1;

	cout << "snp_start:" << endl;
	cout << snp_start << endl;
	cout << "snp_end:" << endl;
	cout << snp_end << endl;

	int marker;
	
	if (WIN_EXT)
	{
		cout << "WIN_EXT is true; trying to extend match" << endl;
		// backwards
		if (start_ms > 0)
		{
			// parameter to scanLeft is called "ms", presumably "marker_set"
			marker = scanLeft(start_ms - 1);
			cout << "marker returned from scanLeft:" << marker << endl;
			snp_start -= (MARKER_SET_SIZE - marker - 2);
			// marker e.g. 24 (how do we get 24? start_ms = 1, so we may be doing scanLeft(0))
			// scanLeft(0) returns 24..?
			// then doing snp_start = 41 - (41 - 24 - 2) -> 41 - 15 -> 26
			// start_ms:
			// 1
			// MARKER_SET_SIZE:
			// 41
			// snp_start:
			// 41
			// snp_end:
			// 163
			// scanLeft marker:
			// 24
			// snp_start after reassigning:
			// 26
			// The homozygous interval stops there -- germline does not seek again to the left of SNP 26
			// This is why the homoz outputs vary by bits-parameter
			cout << "snp_start after reassigning:" << endl;
			cout << snp_start << endl;
		}
	}
	cout << "end_ms: " << end_ms << endl;
	cout << "num_sets: " << num_sets << endl;
	if ( WIN_EXT || end_ms == num_sets - 2 )
	{
		cout << "WIN_EXT or end_ms == num_sets - 2, so we're attempting to extend to the right / forwards" << endl;
		// forwards
		if( end_ms < num_sets - 1 )
		{
			cout << "end_ms < num_sets - 1, so we call scanRight(end_ms + 1)" << endl;
			marker = scanRight( end_ms + 1 );
			snp_end += marker - 1;
		} else {
			cout << "end_ms was not less than num_sets - 1, so we didn't call scanRight" << endl;
		}
	} else {
		cout << "end_ms / num_sets: " << end_ms << " / " << num_sets << endl;
	}
	

	bool genetic;
	float distance;
	if ( ( distance = ALL_SNPS.getDistance(snp_start,snp_end,genetic)) < MIN_MATCH_LEN ) return;
	// print

	// get hamming distance & ignored bit count
	int dif = 0;
	for( unsigned int i = start_ms; i <= end_ms ; i++) { dif += diff( i ); }
	
	// calculate if homozygous
	bool hom[2];
	if ( node[0] == node[1] ) { hom[0] = hom[1] = 1; }
	else
	{
		for ( int n = 0 ; n < 2 ; n++ )
		{
			hom[n] = true;
			for ( unsigned int i = start_ms ; i<= end_ms && hom ; i++ )
			{
				hom[n] = isHom( n , i );
			}
		}
	}

	if ( BINARY_OUT )
	{
		unsigned int pid[2];
		pid[0] = node[0]->getNumericID();
		pid[1] = node[1]->getNumericID();
		unsigned int sid[2];
		sid[0] = ALL_SNPS.getSNP(snp_start).getMarkerNumber();
		sid[1] = ALL_SNPS.getSNP(snp_end).getMarkerNumber();
		fout.write( (char*) &pid[0] , sizeof( unsigned int ) );
		fout.write( (char*) &pid[1] , sizeof( unsigned int ) );
		fout.write( (char*) &sid[0] , sizeof( unsigned int ) );
		fout.write( (char*) &sid[1] , sizeof( unsigned int ) );
		fout.write( (char*) &dif , sizeof( int ) );
		fout.write( (char*) &hom[0] , sizeof( bool ) );
		fout.write( (char*) &hom[1] , sizeof( bool ) );
	} else
	{
		fout << node[0]->getID() << '\t';
		fout << node[1]->getID() << '\t';
		fout << ALL_SNPS.getSNP(snp_start).getChr() << '\t';
		fout << ALL_SNPS.getSNP(snp_start).getPhysPos() << ' ';
		fout << ALL_SNPS.getSNP(snp_end).getPhysPos() << '\t';
		fout << ALL_SNPS.getSNP(snp_start).getSNPID() << ' ';
		fout << ALL_SNPS.getSNP(snp_end).getSNPID() << '\t';
		fout << ( snp_end - snp_start + 1) << '\t';
		fout << setiosflags(ios::fixed) << setprecision(2) << distance << '\t';
		if ( genetic ) fout << "cM" << '\t'; else fout << "MB" << '\t';
		fout << dif;
		for ( int n = 0 ; n < 2 ; n++ )
			if ( hom[n] ) fout << '\t' << 1; else fout << '\t' << 0;
		fout << endl;
	}
	num_matches++;
}

