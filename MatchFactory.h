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
		size_t n = bs.num_blocks();
		// Stack buffer covers -bits up to 1024 (16 × 64-bit blocks).
		// Default -bits 128 uses 2 blocks; -bits 512 uses 8 blocks.
		// Falls back to heap for unusually large word sizes.
		constexpr size_t MAX_STACK_BLOCKS = 16;
		if ( n <= MAX_STACK_BLOCKS )
		{
			block_type buf[MAX_STACK_BLOCKS];
			boost::to_block_range(bs, buf);
			for ( size_t i = 0; i < n; ++i )
				seed ^= hash<block_type>()(buf[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
		else
		{
			vector<block_type> blocks(n);
			boost::to_block_range(bs, blocks.begin());
			for ( auto block : blocks )
				seed ^= hash<block_type>()(block) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
		}
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
