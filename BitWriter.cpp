/*
 * BitWriter.cpp
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#include "BitWriter.h"

namespace std {

BitWriter::BitWriter() {
	blockSize=0;

}

BitWriter::~BitWriter() {
	// TODO Auto-generated destructor stub
}


void BitWriter::init(){
	bitsBuffer=0;
	liveBits=0;
	liveBytes=0;
	if (blockSize==0)
		blockSize=512;
	blockBuffer=new char[blockSize];
	writtenBlocks=0;
	blockBufferEffectiveLength=0;
}


void BitWriter::clearBlockBuffer(){
	for (uint32_t i=0;i<blockSize;i++)
		blockBuffer[i]=0;
}

void BitWriter::putByteIntoBlockBuffer(uint8_t ch){
	if (blockSize==liveBytes){
		writeBlockBuffer();
		clearBlockBuffer();
		liveBytes=0;
		writtenBlocks++;
	}
	blockBuffer[liveBytes]=(char)((short)ch&0xFF);
	liveBytes++;
}

void BitWriter::finishedWithStream() {
	while (liveBits > 0) {
		uint8_t ch = bitsBuffer >> 56;
		putByteIntoBlockBuffer(ch);
		bitsBuffer <<= 8;
		liveBits -= 8;
	}
	if (liveBytes>0){
		writeBlockBuffer();
		clearBlockBuffer();
		liveBytes=0;
		writtenBlocks++;
	}
}

/** Modified to act as fm_bit_write24 function of FM-index V2 implementation
 * Il massimo valore di n deve essere 56, in quanto, per come è implementata la funzione write,
 * nel bit buffer ci sono al più 7 bit live*/

void BitWriter::write(uint8_t n,int64_t v){
	bitsBuffer = bitsBuffer | (v << (64 - liveBits - n));
	liveBits = liveBits + n;

	while (liveBits >= 8) {
		uint8_t ch=bitsBuffer >> 56;   //preleva dal bit buffer gli 8 bit più significativi
		putByteIntoBlockBuffer(ch);
		bitsBuffer <<= 8;		   //elimina dal bit buffer gli 8 bit più significativi, che ormai sono stati scritti su disco
		liveBits -= 8;
	}
}

void BitWriter::writeUByte(uint8_t c){
     write(8, c);
}

void BitWriter::writeInt(int32_t u){
        write(8, (u >> 24) & 0xff);
        write(8, (u >> 16) & 0xff);
        write(8, (u >> 8) & 0xff);
        write(8, u & 0xff);
}

void BitWriter::writeInt64(int64_t u){
		write(8, (u >> 56) & 0xff);
		write(8, (u >> 48) & 0xff);
		write(8, (u >> 40) & 0xff);
		write(8, (u >> 32) & 0xff);
        write(8, (u >> 24) & 0xff);
        write(8, (u >> 16) & 0xff);
        write(8, (u >> 8) & 0xff);
        write(8, u & 0xff);
}


/* -----------------------------------------------------------------------------
Codifica di Interi con una parte fissa ed una variabile.
Questa codifica non so se esisteva ma ha prestazioni buone rispetto a write7x8
Mutuata da FM-Index
------------------------------------------------------------------------------
*/
void BitWriter::integerEncode(int64_t u, uint8_t log2log2Maxvalue){
int64_t  k;
int32_t i = 0;
switch (u){
	/* primi 4 casi speciali */
	case 0: write(log2log2Maxvalue,0);
			break;
	case 1: write(log2log2Maxvalue,1);
			break;
	case 2: write((log2log2Maxvalue+1),4);
			break;
	case 3: write((log2log2Maxvalue+1),5);
			break;
/*	Altri casi. Si calcola i = log2(occ) __numero
	di bit per rappresentare occ. si calcola  k=
	occ-2^i cioe' la distanza dalla potenza di due
	inferiore. Scrive i su log2log2BuckSize bits
	e k su i-1 bits indica __numero in piu' alla
	potenza di due precedente 						*/
	default: {
		int64_t pow = 1;
		/* 	calcola distanza dalla potenza di due
			precedente i indica l'esponente di
			tale potenza */
		i = 3;
		pow = 4;
		while(pow <= u){
			pow= pow<<1;
			i = i+1;
		}
		i--; pow = pow>>1;
		k = u - pow;
		write(log2log2Maxvalue,i);
		write(i-1,k);
	}
  }
}


/**
 * flush the bit buffer: after calling this function we are sure next data written is byte-aligned
 * @throws Exception
 */
void BitWriter::flush(){
	if (liveBits != 0)
		write((8 - (liveBits % 8)), 0);	// pad with zero !
}

uint32_t BitWriter::getBlockSize() const{
	return blockSize;
}

void BitWriter::setBlockSize(uint32_t blockSize) {
	this->blockSize = blockSize;
}








} /* namespace std */
