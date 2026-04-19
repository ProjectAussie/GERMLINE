// MatchFactory.h: Generates matches from individuals

#ifndef MATCHFACTORY_H
#define MATCHFACTORY_H

#include "MarkerSet.h"
#include "Individual.h"
#include <unordered_map>
#include <vector>
#include <functional>

using namespace std;

struct DynamicBitsetHash
{
	size_t operator()(const boost::dynamic_bitset<>& bs) const
	{
		using block_type = boost::dynamic_bitset<>::block_type;
		size_t seed = bs.size();
		vector<block_type> blocks(bs.num_blocks());
		boost::to_block_range(bs, blocks.begin());
		for ( auto block : blocks )
			seed ^= hash<block_type>()(block) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		return seed;
	}
};

class MatchFactory
{

public:

	// MatchFactory(): default constructor
	// Precondition: None.
	// Postcondition: segments and matches are empty and position is -1.
	MatchFactory();

	int size();

    // initialize(): initializes object
	// Precondition:  None.
	// Postcondition: If 0=<pos, then position is set to pos and map is empty;
	//  otherwise an error message is printed.
	void initialize();

	void hash(Individual *);
	void assertShares();

private:

	// stores data to check for matches
	unordered_map < boost::dynamic_bitset<> , Share, DynamicBitsetHash > segments;
	unordered_map < boost::dynamic_bitset<> , Share, DynamicBitsetHash >::iterator iter;
};

#endif

// end MatchFactory.h
