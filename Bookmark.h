/*
 * Bookmark.h
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#ifndef BOOKMARK_H_
#define BOOKMARK_H_
#include <string>

namespace std {

class Bookmark {
public:
	Bookmark();

	Bookmark(uint64_t bitAbsolutePosition);

	Bookmark(uint64_t bytePosition,uint8_t liveBits);

	Bookmark(uint64_t bytePosition,int64_t bitsBuffer,uint8_t liveBits);

	virtual ~Bookmark();

	uint64_t getPosition();

	void setPosition(uint64_t position);

	int64_t getBitsBuffer();

	void setBitsBuffer(int64_t bitsBuffer);

	uint8_t getLiveBits();

	void setLiveBits(uint8_t liveBits);

	/*
	 * Alla posizione 8*position (allineata al byte) tolgo i bit ancora da leggere
	 */
	uint64_t getAbsoluteBitPosition() const;

	string toString() const;

	friend ostream& operator<<(ostream& os, const Bookmark& b);

	friend class MemoryBitReader;
	friend class FileBitReader;

private:
	uint64_t position;    //file pointer position
	int64_t bitsBuffer;  //bits buffer
	uint8_t liveBits;	 //live bits in bit buffer
};

} /* namespace std */

#endif /* BOOKMARK_H_ */
