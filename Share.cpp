
#include "Individual.h"

Share::Share( Individual * cip )
{
	add( cip );
}

void Share::assertMatches()
{
	list<Individual*>::iterator i , ii;
	Match * m;

	for ( i = matches.begin() ; i != matches.end() ; i++ )
	{
		if (DEBUG) cout << "i in share::assertMatches outer loop: " << *i << endl;
		ii = i;
		for ( ++ii ; ii != matches.end() ; ii++ )
		{
                  if ((ALL_SAMPLES.isOld(**i) && ALL_SAMPLES.isNew(**ii)) || (ALL_SAMPLES.isNew(**i) && ALL_SAMPLES.isOld(**ii))) {
				continue;
			}
			if (DEBUG) cout << "ii in share::assertMatches inner loop: " << *ii << endl;
			// Check if this pair matched in previous word (symmetrically)
			m = (*i)->getMatch( (*ii)->getNumericID() );
			if (DEBUG) cout << "numericID in share::assertMatches inner loop: " << (*ii)->getNumericID() << endl;
			if ( m == NULL ) m = (*ii)->getMatch( (*i)->getNumericID() );

			if ( m != NULL )
			{
				// This match can be incremented
				if (DEBUG) cout << "incrementing match end_ms to " << position_ms << endl;
				m->end_ms = position_ms;
			}
			else
			{
				// This match must be created
				m = createMatch( *i , *ii );
				// Extend the match backwards
				if (DEBUG) cout << "extendBack() called from Share.cpp" << endl;
				m->extendBack();
				// Mark asserted
				(*i)->addMatch( (*ii)->getNumericID() , m );
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
