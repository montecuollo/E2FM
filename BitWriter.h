/*
 * BitWriter.h
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#ifndef BITWRITER_H_
#define BITWRITER_H_
#include "Bookmark.h"

namespace std {

class BitWriter {
public:
	BitWriter();
	virtual ~BitWriter();


	/** Get the bookmark to the actual position **/
	virtual Bookmark *getBookmark()=0;

	/**
	 * Goto bookmark, restoring file position and bit buffer state
	 * @param bookmark
	 * @throws IOException
	 */
	virtual void gotoBookmark(Bookmark *bookmark)=0;

	/** Write the block buffer on resilient storage */
	virtual void writeBlockBuffer() = 0;

	virtual void open() =0;

	virtual void close() =0;

	virtual void putByteIntoBlockBuffer(uint8_t ch);

	void finishedWithStream();

	void write(uint8_t n, int64_t v);

	void writeUByte(uint8_t c);

	void writeInt(int32_t u);

	void writeInt64(int64_t u);


	/* -----------------------------------------------------------------------------
	   Codifica di Interi con una parte fissa ed una variabile.
	   Ha prestazioni migliori rispetto a write7x8
	   Mutuata da FM-Index
	   ------------------------------------------------------------------------------
	 */

	 void integerEncode(int64_t u,uint8_t log2log2Maxvalue);
	/**
	 * flush the bit buffer: after calling this function we are sure next data written is byte-aligned
	 * @throws Exception
	 */
	void flush();

	uint32_t getBlockSize() const;

	void setBlockSize(uint32_t blockSize);


protected:
	void init();
	void clearBlockBuffer();

	int64_t bitsBuffer;
	int8_t liveBits;
	uint32_t blockSize;
	char *blockBuffer;
	int32_t blockBufferEffectiveLength;
	int32_t liveBytes;
	int32_t writtenBlocks;
};

} /* namespace std */

#endif /* BITWRITER_H_ */
