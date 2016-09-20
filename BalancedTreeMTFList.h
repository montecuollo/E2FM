/*
 * BalancedTreeMTFList.h
 *
 *  Created on: 14/ago/2015
 *      Author: fernando
 */

#ifndef BALANCEDTREEMTFLIST_H_
#define BALANCEDTREEMTFLIST_H_
#include <stddef.h>

namespace std {


class Pair{
public:
	Pair(int64_t timestamp,uint32_t character){
		this->timestamp=timestamp;
		this->character=character;
	}

	int64_t timestamp;
	uint32_t character;  //character associated to the timestamp
};

enum BalancedTreeNodeType
{
    empty = 0,
    twoNode,
    threeNode
};


class BalancedTreeNode{
public:
	BalancedTreeNode(){
		smallKey=NULL;
		bigKey=NULL;
		type=empty;
		parent=NULL;
		leftChild=NULL;
		middleChild=NULL;
		rightChild=NULL;
		size=0;
	}

	BalancedTreeNode(BalancedTreeNode *leftChild, BalancedTreeNode * rightChild, Pair *smallKey){
			this->smallKey=smallKey;
			bigKey=NULL;
			type=twoNode;
			parent=NULL;
			this->leftChild=leftChild;
			middleChild=NULL;
			this->rightChild=rightChild;

			size=1;
			if (leftChild!=NULL){
				leftChild->parent=this;
				size+=leftChild->size;
			}
			if (rightChild!=NULL){
				rightChild->parent=this;
				size+=rightChild->size;
			}


	}

	void addPair(Pair *pair);
	BalancedTreeNode* searchInsertionLeaf(int64_t timestamp);
	BalancedTreeNode* search(int64_t timestamp,int64_t *index);

	inline void updateNodeSize();
	void dumpSubTree(int64_t level);

	//Return the element having a given position in the sorted list rooted in this node
	uint32_t select(int64_t position,BalancedTreeNode **node,Pair **pair);


	Pair *smallKey;
	Pair *bigKey;
	BalancedTreeNodeType type;
	BalancedTreeNode *parent;
	BalancedTreeNode *leftChild;
	BalancedTreeNode *middleChild;
	BalancedTreeNode *rightChild;

	uint32_t size;



};


class BalancedTreeMTFList {
public:
	BalancedTreeMTFList(uint32_t alphabetSize);
	virtual ~BalancedTreeMTFList();
	int64_t indexOf(uint32_t character,BalancedTreeNode **node,Pair **pair);
	int64_t moveToFront(uint32_t character,int64_t newTimestamp,BalancedTreeNode *node,Pair *pair); //used in compression
	uint32_t getFirstCharacter();  //get(0) returns the element in position 0
	uint32_t getCharacter(int64_t position,BalancedTreeNode **node,Pair **pair); //returns the character in a given position of the list

	void dumpList();
	void dumpTree();
private:
	uint32_t alphabetSize;
	int *timestamps; //array of timestamps for each character
	BalancedTreeNode *root;
	BalancedTreeNode *leftmostLeaf;
	inline void insertIntoTree(int64_t timestamp,uint32_t character);
	inline void deleteFromTree(int64_t timestamp);
	inline void deleteFromTree(int64_t timestamp,BalancedTreeNode* node,Pair* pair);
	inline BalancedTreeNode* swapWithSuccessor(BalancedTreeNode* node,Pair* pair);
	inline BalancedTreeNode* getSuccessorLeaf(BalancedTreeNode* node, Pair *pKey);
	inline bool tryLocalRotation(BalancedTreeNode *node);
	inline bool tryLocalMerge(BalancedTreeNode *node);
	inline bool tryGlobalMerge(BalancedTreeNode *pNode,BalancedTreeNode **ppNextNodeToFix);

	inline BalancedTreeNode* searchInTree(int64_t timestamp, Pair **pair, int64_t *index);

	inline Pair * splitNode(Pair * pair,BalancedTreeNode * nodeToSplit,BalancedTreeNode ** addedNode);



};

} /* namespace std */

#endif /* BALANCEDTREEMTFLIST_H_ */
