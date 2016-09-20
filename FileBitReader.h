/*
 * FileBitReader.h
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#ifndef FILEBITREADER_H_
#define FILEBITREADER_H_
#include <fstream>
#include "Bookmark.h"
#include "BitReader.h"
namespace std {

class FileBitReader: public BitReader {
public:

	FileBitReader(string fileName);
	FileBitReader();

	virtual ~FileBitReader();

	void open();

	void open(string fileName);

	void close();

	uint32_t readBlockBuffer();

	virtual void gotoNextByteStart();

	void gotoBookmark(Bookmark *bookmark);

	Bookmark *getBookmark();

	virtual uint8_t getByteFromBlockBuffer(bool *errorFlag);



	virtual int32_t getInt();

	virtual int64_t getInt64();

protected:
	string fileName;
	ifstream in;
};

} /* namespace std */

#endif /* FILEBITREADER_H_ */
