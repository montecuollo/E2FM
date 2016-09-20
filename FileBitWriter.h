/*
 * FileBitWriter.h
 *
 *  Created on: 08/ott/2014
 *      Author: fernando
 */

#ifndef FILEBITWRITER_H_
#define FILEBITWRITER_H_

#include "Bookmark.h"
#include "BitWriter.h"
#include <fstream>

namespace std {

class FileBitWriter: public BitWriter {
public:

	/* Needed to have a ofstream class member */
	FileBitWriter(const FileBitWriter&) = delete;
	FileBitWriter& operator=(const FileBitWriter&) = delete;


	virtual ~FileBitWriter();

	FileBitWriter(string fileName);

	virtual Bookmark *getBookmark();

	/**
	 * Goto bookmark, restoring file position and bit buffer state
	 * @param bookmark
	 * @throws IOException
	 */

	void gotoBookmark(Bookmark *bookmark);

	void open();

	void close();

	void writeBlockBuffer();


protected:
	 string fileName;
	 ofstream out;

};

} /* namespace std */

#endif /* FILEBITWRITER_H_ */
