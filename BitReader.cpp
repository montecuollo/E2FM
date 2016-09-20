/*
 * BitReader.cpp
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#include "BitReader.h"
#include <stdexcept>

namespace std {

BitReader::BitReader() {
	blockSize=0;
}

BitReader::~BitReader() {
	// TODO Auto-generated destructor stub
}


	void BitReader::open(){
		bitsBuffer=0;
		liveBits=0;
    }


    void BitReader::close(){

    }


    void BitReader::init(){
	}

	int64_t BitReader::read(uint8_t n){
        int8_t bsLiveShadow = liveBits;
        int64_t bsBuffShadow = bitsBuffer;
        bool errorFlag;

        if (bsLiveShadow < n) {
            do {
            	uint8_t thech = getByteFromBlockBuffer(&errorFlag);
                if (errorFlag)
                    throw runtime_error("unexpected end of stream");
                //accoda i bit del carattere appena letto all'attuale contenuto del bit buffer
                bsBuffShadow = (bsBuffShadow << 8) | thech;
                bsLiveShadow += 8;
            } while (bsLiveShadow < n);

            bitsBuffer = bsBuffShadow;
        }
        liveBits = bsLiveShadow - n;
        //effettua uno shift right in modo che i primi n bit si trovino a destra
        //ed effettua l'AND bit a bit con una maschera costituita da n cifre 1 a destra e con le rimanenti cifre a 0
        int64_t shifted=(bsBuffShadow >> (bsLiveShadow - n));
        int64_t mask=1;
        mask=(mask << n)-1;

        int64_t retVal=shifted & mask;
        return retVal;
    }





    /*
	bool BitReader::getBit(){
    	int8_t bsLiveShadow = liveBits;
    	int64_t bsBuffShadow = bitsBuffer;
    	bool errorFlag;

        if (bsLiveShadow < 1) {
        	uint8_t thech=getByteFromBlockBuffer(&errorFlag);

            if (errorFlag) {
                throw runtime_error("unexpected end of stream");
            }

            bsBuffShadow = (bsBuffShadow << 8) | thech;
            bsLiveShadow += 8;
            bitsBuffer = bsBuffShadow;
        }

        liveBits = bsLiveShadow - 1;
        return ((bsBuffShadow >> (bsLiveShadow - 1)) & 1) != 0;
    }*/


    void BitReader::clearBlockBuffer(){
    	for (uint32_t i=0;i<blockSize;i++)
    		blockBuffer[i]=0;
    }




    uint8_t BitReader::getUByte() {
        return  read(8);
    }


    int64_t BitReader::integerDecode (uint8_t headBits){
    	int64_t k;
    	int32_t i;
    	i=read(headBits);
    	/*
    	 * casi base speciali
    	 */

    	if (i == 0)
    		return 0;
    	if (i == 1)
    		return 1;
    	if (i == 2)
    	{
    		i=read (1);
    		k = 2 + i;
    		return k;
    	}
    	/*
    	 * i indica il numero di bits che seguono num = 2^i + k
    	 */
    	i--;
    	/*
    	 * bit_read24 e' corretta fino a superbuckets di dimensione inferiore
    	 * a 16384*1024. E' comunque meglio utilizzare questa perche molto
    	 * piu' veloce
    	 */
    	k=read (i);
    	return ((1 << i) + k);
    }


	uint32_t BitReader::getBlockSize() {
		return blockSize;
	}

	void BitReader::setBlockSize(uint32_t blockSize) {
		if (blockSize!=0)
			this->blockSize = blockSize;
		else
			this->blockSize=512;
	}


} /* namespace std */


