/*
 * BitReader.h
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#ifndef BITREADER_H_
#define BITREADER_H_

#include "Bookmark.h"
namespace std {

class BitReader {
public:
	BitReader();
	virtual ~BitReader();

	/**
	 *
	 * @return a bookmark to the actual position in the file (0-based)
	 * @throws Exception
	 */
	virtual Bookmark *getBookmark() = 0;

	/**
	 * Goes to the bookmark passed as parameter
	 * The position contained in the bookmark must be 0-based
	 * @param bookmark
	 * @throws Exception
	 */
	virtual void gotoBookmark(Bookmark *bookmark)=0;


	//Gets a byte from the block buffer
	virtual uint8_t getByteFromBlockBuffer(bool *errorFlag)=0;

	/**
	 * Retrieves a block of size blockSize from the data and stores it
	 * in the block buffer
	 * @return the real size of the blockBuffer, less than blockSize if
	 *         there aren't enough data
	 * @throws Exception
	 */
	//virtual uint32_t readBlockBuffer()=0;

	virtual void open();


	virtual void close();

	void init();

	virtual int64_t read(uint8_t n);

	//bool getBit();

	void clearBlockBuffer();

	uint8_t getUByte();

	virtual int32_t getInt()=0;

	virtual int64_t getInt64()=0;

	int64_t integerDecode (uint8_t headBits);

	virtual void gotoNextByteStart()=0;

	uint32_t getBlockSize();

	void setBlockSize(uint32_t blockSize);

protected:
	int64_t bitsBuffer;
	int8_t liveBits;
	uint32_t blockSize;
	char *blockBuffer;
	int64_t blockBufferEffectiveLength;
	int64_t liveBytes;
	int64_t readBlocks;

};

} /* namespace std */

#endif /* BITREADER_H_ */
