/*
 * Bookmark.cpp
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#include "Bookmark.h"
#include <math.h>
#include <string>
#include <sstream>

namespace std {

    /*
    Bookmark::Bookmark(uint64_t bitAbsolutePosition){
		position=(uint64_t) ceil((double)(bitAbsolutePosition+1)/(double)8);
		liveBits=8 - bitAbsolutePosition%8;
		bitsBuffer=0;
	}*/

    Bookmark::Bookmark(uint64_t bitAbsolutePosition){
		position=bitAbsolutePosition/8;
		liveBits=8 - bitAbsolutePosition%8;
		bitsBuffer=0;
	}

	Bookmark::Bookmark(uint64_t bytePosition,uint8_t liveBits){
		this->position=bytePosition;
		this->liveBits=liveBits;
		this->bitsBuffer=0;
	}

	Bookmark::Bookmark(uint64_t bytePosition,int64_t bitsBuffer,uint8_t liveBits){
		this->position=bytePosition;
		this->bitsBuffer=bitsBuffer;
		this->liveBits=liveBits;
	}

	Bookmark::~Bookmark() {
		// TODO Auto-generated destructor stub
	}


	uint64_t Bookmark::getPosition() {
		return position;
	}
	void Bookmark::setPosition(uint64_t position) {
		this->position = position;
	}

	int64_t Bookmark::getBitsBuffer() {
		return bitsBuffer;
	}

	void Bookmark::setBitsBuffer(int64_t bitsBuffer) {
		this->bitsBuffer = bitsBuffer;
	}

	uint8_t Bookmark::getLiveBits() {
		return liveBits;
	}

	void Bookmark::setLiveBits(uint8_t liveBits) {
		this->liveBits = liveBits;
	}

	/*
	 * Alla posizione 8*position (allineata al byte) tolgo i bit ancora da leggere
	 */
	uint64_t Bookmark::getAbsoluteBitPosition() const{
		return 8*position - liveBits;
	}

	string Bookmark::toString() const{
		stringstream stream;
		stream << "[" << position  << "," << liveBits << "] <=> " << getAbsoluteBitPosition();
		string  s = stream.str();
		return s;
	}


	ostream& operator<<(ostream& os,  const Bookmark& b) {
			os << b.toString();
			return os;
	}


} /* namespace std */
