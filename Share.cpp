
#include "Individual.h"

Share::Share( Individual * cip )
{
	add( cip );
}

void Share::assertMatches()
{
	Match * m;
	for ( size_t i = 0 ; i < matches.size() ; i++ )
	{
		if (DEBUG) cout << "i in share::assertMatches outer loop: " << matches[i] << endl;
		for ( size_t ii = i + 1 ; ii < matches.size() ; ii++ )
		{
			if ((!ALL_SAMPLES.hasRestrictions()) || (matches[i]->is_old && matches[ii]->is_new) || (matches[i]->is_new && matches[ii]->is_old)) {
				if (DEBUG) cout << "ii in share::assertMatches inner loop: " << matches[ii] << endl;
				// Check if this pair matched in previous word (symmetrically)
				m = matches[i]->getMatch( matches[ii]->getNumericID() );
				if (DEBUG) cout << "numericID in share::assertMatches inner loop: " << matches[ii]->getNumericID() << endl;
				if ( m == NULL ) m = matches[ii]->getMatch( matches[i]->getNumericID() );
				if ( m != NULL )
				{
					// This match can be incremented
					if (DEBUG) cout << "incrementing match end_ms to " << position_ms << endl;
					m->end_ms = position_ms;
				}
				else
				{
					// This match must be created
					m = createMatch( matches[i] , matches[ii] );
					// Extend the match backwards
					if (DEBUG) cout << "extendBack() called from Share.cpp" << endl;
					m->extendBack();
					// Mark asserted
					matches[i]->addMatch( matches[ii]->getNumericID() , m );
				}
			}
		}
	}
}


Match * Share::createMatch(Individual * c1 , Individual * c2)
{
	Match * new_match = new Match();
	new_match->end_ms = new_match->start_ms = position_ms;

	new_match->node[0] = c1;
	new_match->node[1] = c2;

	return new_match;
}

void Share::add(Individual * cip)
{
	matches.push_back( cip );
}
