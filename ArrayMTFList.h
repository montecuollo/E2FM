/*
 * ArrayMTFList.h
 *
 *  Created on: 24/dic/2014
 *      Author: fernando
 */

#ifndef ARRAYMTFLIST_H_
#define ARRAYMTFLIST_H_

#include <stdint.h>

namespace std {

class ArrayMTFList {
public:
	ArrayMTFList(uint32_t alphabetSize);
	virtual ~ArrayMTFList();

	uint32_t moveToFront(uint32_t position);
	uint32_t get(uint32_t position);
	int64_t indexOf(uint32_t item);
	void dumpList();

private:
	uint32_t *elementData;
	uint32_t alphabetSize;
};

} /* namespace std */

#endif /* ARRAYMTFLIST_H_ */
