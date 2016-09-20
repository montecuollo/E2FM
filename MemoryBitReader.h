/*
 * MemoryBitReader.h
 *
 *  Created on: 19/nov/2015
 *      Author: fernando
 */

#ifndef MEMORYBITREADER_H_
#define MEMORYBITREADER_H_
#include "BitReader.h"
namespace std {

class MemoryBitReader : public BitReader{
public:
	MemoryBitReader(char *data,uint64_t size);
	MemoryBitReader();
	virtual ~MemoryBitReader();

	void open();

	void close();

	void gotoBookmark(Bookmark *bookmark);

	virtual void gotoNextByteStart();

	Bookmark *getBookmark();

	virtual uint8_t getByteFromBlockBuffer(bool *errorFlag);

	virtual int64_t read(uint8_t n);

	virtual int32_t getInt();

	virtual int64_t getInt64();

	static char *loadDataFromFile(string filePath,uint64_t &size);



protected:
	int64_t position;
};

} /* namespace std */

#endif /* MEMORYBITREADER_H_ */
