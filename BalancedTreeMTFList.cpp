/*
 * BalancedTreeMTFList.cpp
 *
 *  Created on: 14/ago/2015
 *      Author: fernando
 */
#include <chrono>
#include "BalancedTreeMTFList.h"
#include <string.h>
#include <iostream>
namespace std {

BalancedTreeMTFList::BalancedTreeMTFList(uint32_t alphabetSize) {
	this->alphabetSize=alphabetSize;
	//initialize timestamps for each alphabet's character and initialize tree
	timestamps=new int[alphabetSize];
	root=NULL;
	leftmostLeaf=NULL;
	for (int64_t c=alphabetSize-1;c>=0;c--){
		timestamps[c]=c;
		insertIntoTree(c,c);
		//dumpTree();
	}

	 //esempio di utilizzo per la Move to Front
/*
	Pair *pair;
	BalancedTreeNode *node;
	int64_t index=indexOf(3,&node,&pair);
	moveToFront(3,-1,node,pair);
	dumpTree();
	uint32_t c=getFirstCharacter();

	for (uint32_t i=0;i<alphabetSize;i++){
		c=getCharacter(i,&node,&pair);
		cout << c << " " << pair->timestamp <<endl;
	}

	index=indexOf(2,&node,&pair);
	moveToFront(2,-2,node,pair);
	dumpTree();
	c=getFirstCharacter();

	index=indexOf(2,&node,&pair);
	moveToFront(2,-2,node,pair);
	dumpTree();
	c=getFirstCharacter();

	index=indexOf(3,&node,&pair);
	moveToFront(3,-3,node,pair);
	dumpTree();
	c=getFirstCharacter();

	index=indexOf(4,&node,&pair);
	moveToFront(4,-4,node,pair);
	dumpTree();
	c=getFirstCharacter();

*/
}



BalancedTreeMTFList::~BalancedTreeMTFList() {
	delete[] timestamps;
}



inline void BalancedTreeNode::updateNodeSize(){
	if (type==twoNode)
		size=1+(leftChild!=NULL?leftChild->size:0) + (rightChild!=NULL?rightChild->size:0);
	else if (type=threeNode)
		size=2+(leftChild!=NULL?leftChild->size:0) + (middleChild!=NULL?middleChild->size:0) +
		(rightChild!=NULL?rightChild->size:0);
	else {
		size=0;
	}
}

void BalancedTreeMTFList::insertIntoTree(int64_t timestamp, uint32_t character) {
	Pair *pair=new Pair(timestamp,character);

	/* If no root present, create one. */
	if(root == NULL)
	{
		root = new BalancedTreeNode();
		root->addPair(pair);
		leftmostLeaf=root;
	}
	else{
		/* Find place in tree to insert key.*/
		//BalancedTreeNode* insertionLeaf=root->searchInsertionLeaf(timestamp);
		BalancedTreeNode* insertionLeaf=leftmostLeaf;

		if(insertionLeaf->type == twoNode){
			insertionLeaf->addPair(pair);
			BalancedTreeNode *p=insertionLeaf;
			while (p->parent!=NULL){
				p->parent->size++;
				p=p->parent;
			}
		}
			/* In case of three node we need to split node. */
		else if(insertionLeaf->type == threeNode)
		{
			bool splittingNotNeeded=false;
			BalancedTreeNode* addedNode=NULL;

			Pair *p= pair;
			BalancedTreeNode *nodeToSplit=insertionLeaf;
			bool leafLevel=true;
			do
			{
				/* Split node. */
				p=splitNode(p,nodeToSplit,&addedNode);

				if (leafLevel){
					leftmostLeaf=addedNode;
					leafLevel=false;
				}

				/* Node splitting has propagated up to root node, create new root. */
				if(nodeToSplit->parent == NULL)
				{
					BalancedTreeNode* newRoot = new BalancedTreeNode(addedNode, nodeToSplit, p);
					root = newRoot;
					splittingNotNeeded = true;
				}
				else
				{
					nodeToSplit = nodeToSplit->parent;
					if(nodeToSplit->type == twoNode)
					{
						if( p->timestamp < nodeToSplit->smallKey->timestamp)
						{
							nodeToSplit->addPair(p);
							nodeToSplit->middleChild=nodeToSplit->leftChild;
							nodeToSplit->leftChild=addedNode;
							//nodeToSplit->size+=addedNode->size;
							splittingNotNeeded = true;
							nodeToSplit->updateNodeSize();

						}
						else if( p->timestamp > nodeToSplit->smallKey->timestamp)
						{
							nodeToSplit->addPair(p);
							nodeToSplit->middleChild=addedNode;
							splittingNotNeeded = true;
							//nodeToSplit->size+=addedNode->size;
							nodeToSplit->updateNodeSize();
						}

					}

				}
			} while (!splittingNotNeeded);

			//update all the ancestor node sizes
			while (nodeToSplit->parent!=NULL){
				nodeToSplit->parent->updateNodeSize();
				nodeToSplit=nodeToSplit->parent;
			}
		}
	}
}

Pair * BalancedTreeMTFList::splitNode(Pair* pair, BalancedTreeNode* nodeToSplit,BalancedTreeNode ** addedNode) {

		BalancedTreeNode *node = NULL;
		Pair *toUpperLevel;
		/* Three ways to split a node, depending on key value. */
		if (pair->timestamp > nodeToSplit->bigKey->timestamp){
			node=new BalancedTreeNode(nodeToSplit->leftChild,
								 nodeToSplit->middleChild,
								 nodeToSplit->smallKey);
			nodeToSplit->smallKey=pair;
			toUpperLevel=nodeToSplit->bigKey;
			nodeToSplit->leftChild=*addedNode;
		}
		else if(pair->timestamp < nodeToSplit->smallKey->timestamp){
			node = new BalancedTreeNode(*addedNode,
										nodeToSplit->leftChild,
										pair);
			toUpperLevel=nodeToSplit->smallKey;
			nodeToSplit->smallKey=nodeToSplit->bigKey;
			nodeToSplit->leftChild=nodeToSplit->middleChild;
		}
		else{ /* Split in the middle. */
			node = new BalancedTreeNode(nodeToSplit->leftChild,
					*addedNode,
					nodeToSplit->smallKey);
			nodeToSplit->smallKey=nodeToSplit->bigKey;
			nodeToSplit->leftChild= nodeToSplit->middleChild;
			toUpperLevel=pair;

		}
		nodeToSplit->bigKey=NULL;
		nodeToSplit->middleChild=NULL;
		nodeToSplit->type=twoNode;
		node->parent=nodeToSplit->parent;
		//update the size value
		nodeToSplit->size=1;
		if (nodeToSplit->leftChild!=NULL)
			nodeToSplit->size+=nodeToSplit->leftChild->size;
		if (nodeToSplit->rightChild!=NULL)
			nodeToSplit->size+=nodeToSplit->rightChild->size;
		*addedNode=node;
		return toUpperLevel;
}




void BalancedTreeNode::addPair(Pair *pair){
    if (type == empty){
        smallKey=pair;
        type=twoNode;
    }
    else if (type == twoNode && pair->timestamp > smallKey->timestamp){
        bigKey=pair;
        type=threeNode;
    }
    else if (type == twoNode && pair->timestamp < smallKey->timestamp){
        bigKey=smallKey;
        smallKey=pair;
        type=threeNode;
    }
    size++;
}

BalancedTreeNode* BalancedTreeNode::searchInsertionLeaf(int64_t timestamp) {
	if (leftChild == NULL && rightChild == NULL)  //spostato dopo il while
			return this;

	BalancedTreeNode *searchNode;
	if (type==twoNode){
		if (timestamp > smallKey->timestamp)
			searchNode=rightChild;
		else
			searchNode=leftChild;
	} else {  //three-node
		if (timestamp > bigKey->timestamp)
			searchNode=rightChild;
		else if (timestamp > smallKey->timestamp)
			searchNode=middleChild;
		else
			searchNode=leftChild;
	}
	return searchNode->searchInsertionLeaf(timestamp);

}

BalancedTreeNode* BalancedTreeNode::search(int64_t timestamp,int64_t *index) {
	if (timestamp == smallKey->timestamp){
		if (leftChild!=NULL)
			*index=leftChild->size;
		else
			*index=0;
		return this;
	} else if (type==threeNode && timestamp == bigKey->timestamp){
		*index=1;
		if (leftChild!=NULL)
			*index+=leftChild->size;
		if (middleChild!=NULL)
			*index+=middleChild->size;
		return this;
	}

	//if not found in a leaf, then return NULL: the key doesn't exist!
	if (leftChild == NULL && rightChild == NULL)
		return NULL;

	if (timestamp < smallKey->timestamp)
		return leftChild->search(timestamp,index);
	else if (type==threeNode && timestamp > bigKey->timestamp ||
			 type==twoNode && timestamp > smallKey->timestamp){
		int64_t rightSubTreeIndex;
		BalancedTreeNode* node=rightChild->search(timestamp,&rightSubTreeIndex);
		*index=rightSubTreeIndex+1; //there is at least the smallKey at the left of rightSubTree
		if (leftChild!=NULL)
			*index+=leftChild->size;
		if (middleChild!=NULL)
			*index+=middleChild->size;
		if (bigKey!=NULL)
			*index+=1;
		return node;
	}
	else{
		int64_t middleSubTreeIndex;
		BalancedTreeNode* node=middleChild->search(timestamp,&middleSubTreeIndex);
		*index=middleSubTreeIndex+1; //there is at least the smallKey at the left of middleSubTree
		if (leftChild!=NULL)
			*index+=leftChild->size;
		return node;

	}

}

uint32_t BalancedTreeNode::select(int64_t position,BalancedTreeNode **node,Pair **pair){
	int64_t ls;  //left child size
	if (leftChild!=NULL)
		ls=leftChild->size;
	else
		ls=0;
	if (position == ls){
		*node=this;
		*pair=smallKey;
		return smallKey->character;
	}
	else if (position < ls)
		return leftChild->select(position,node,pair);
	else {
		if (type==twoNode)
			//to compute the position within the right child we must subtract from the original position the number
			//of the elements preceding the right child elements (i.e. the left child elements and smallkey)
			return rightChild->select(position-(ls+1),node,pair);
		else{
			int64_t ms;
			if (middleChild!=NULL)
				ms=middleChild->size;
			else
				ms=0;
			if (position < ls + 1 + ms)
				middleChild->select(position - (ls+1),node,pair);
			else if (position == ls + 1 +ms){
				*node=this;
				*pair=bigKey;
				return bigKey->character;
			}
			else
				rightChild->select(position-(ls+ms+2),node,pair);
		}
	}
}



uint32_t BalancedTreeMTFList::getFirstCharacter() {
	return leftmostLeaf->smallKey->character;
}

int64_t BalancedTreeMTFList::indexOf(uint32_t character, BalancedTreeNode** node, Pair **pair) {
	int64_t timestamp=timestamps[character];

	int64_t index;
	*node=searchInTree(timestamp,pair,&index);
	return index;
}









int64_t BalancedTreeMTFList::moveToFront(uint32_t character, int64_t newTimestamp,BalancedTreeNode* node,Pair *pair) {
	//std::chrono::time_point<std::chrono::system_clock> startTime,endTime;
	//startTime=std::chrono::system_clock::now();
	deleteFromTree(pair->timestamp,node,pair);
	//endTime=std::chrono::system_clock::now();
	//double deleteTime=(double)std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()/(double)1000;
	//cout << "delete time: " <<deleteTime << endl;

	//startTime=std::chrono::system_clock::now();
	insertIntoTree(newTimestamp,character);
	//endTime=std::chrono::system_clock::now();
	//double insertTime=(double)std::chrono::duration_cast<std::chrono::nanoseconds>(endTime - startTime).count()/(double)1000;
	//cout << "insert time: " <<insertTime << endl;

	timestamps[character]=newTimestamp;

}







/*
 * Returns the character in a given position of the list
 */
uint32_t BalancedTreeMTFList::getCharacter(int64_t position,BalancedTreeNode **node,Pair **pair){
	return root->select(position,node,pair);
}



void BalancedTreeMTFList::dumpTree() {
	if (root != NULL)
		root->dumpSubTree(0);
	cout <<endl;
}


void BalancedTreeNode::dumpSubTree(int64_t level){
	cout << std::string(3*level, ' ');
	if (type==twoNode){
		cout << "[" << smallKey->timestamp << "," << smallKey->character << "]" << "(" << size<< ")" << endl;
		if (leftChild!=NULL)
			leftChild->dumpSubTree(level+1);

		if (rightChild!=NULL)
			rightChild->dumpSubTree(level+1);
	} else{
		cout << "[" << smallKey->timestamp << "," << smallKey->character << "]"  ;
		cout << "[" << bigKey->timestamp << "," << bigKey->character << "]" << "(" << size<< ")" <<endl;
		if (leftChild!=NULL)
					leftChild->dumpSubTree(level+1);
		if (middleChild!=NULL)
			middleChild->dumpSubTree(level+1);
		if (rightChild!=NULL)
			rightChild->dumpSubTree(level+1);
	}

}


void BalancedTreeMTFList::deleteFromTree(int64_t timestamp,BalancedTreeNode *node,Pair *pair) {
	BalancedTreeNode *nextNodeToFix=NULL;
	if (node->leftChild!=NULL)
			node=swapWithSuccessor(node,pair);

		if (node->type == threeNode){ //simply remove from leaf its small key
			if (timestamp == node->smallKey->timestamp){
				node->smallKey=node->bigKey;
				node->bigKey=NULL;

			} else
				node->bigKey=NULL;
			node->type=twoNode;
		} else{
			bool fixed=false;
			do{
				if (node == root){ //the node to delete is the tree root
					root=node->leftChild;
					if (root!=NULL)
						root->parent=NULL;
					fixed=true;
				} else {
					if (tryLocalRotation(node) || tryLocalMerge(node))
						fixed=true;
					else if (tryGlobalMerge(node,&nextNodeToFix)){
						node=nextNodeToFix;
					}

				}
			} while (!fixed);
		}

		//update all the ancestor node sizes
		do{
			node->updateNodeSize();
			node=node->parent;
		}while (node!=NULL);
}

void BalancedTreeMTFList::deleteFromTree(int64_t timestamp) {
	Pair *pair;
	int64_t index;
	BalancedTreeNode *node=searchInTree(timestamp,&pair,&index);

	deleteFromTree(timestamp,node,pair);

}



bool BalancedTreeMTFList::tryLocalRotation(BalancedTreeNode *pNode){  //successful if the node has a sibling with two keys
	BalancedTreeNode* pParentNode = NULL;
	BalancedTreeNode* pRightNode = NULL;
	BalancedTreeNode* pMiddleNode = NULL;
	BalancedTreeNode* pLeftNode = NULL;
	bool fixed = false;

	pParentNode = pNode->parent;
	pRightNode = pParentNode->rightChild;
	pMiddleNode = pParentNode->middleChild;
	pLeftNode = pParentNode->leftChild;
	/* Case 1: Take from right. */
	if(pNode == pLeftNode &&
	   pParentNode->type == twoNode &&
	   pRightNode->type == threeNode)
	{
		pLeftNode->smallKey=pParentNode->smallKey;
		pLeftNode->rightChild=pRightNode->leftChild;
		if (pLeftNode->rightChild!=NULL)
			pLeftNode->rightChild->parent=pLeftNode;
		pLeftNode->type=twoNode;
		pParentNode->smallKey=pRightNode->smallKey;
		pRightNode->smallKey=pRightNode->bigKey;
		pRightNode->bigKey=NULL;
		pRightNode->leftChild= pRightNode->middleChild;
		pRightNode->middleChild=NULL;
		pRightNode->type=twoNode;
		fixed = true;
		pLeftNode->updateNodeSize();
		pRightNode->updateNodeSize();
	}
	/* Case 2: Take from middle first. */
	else if(pNode == pLeftNode &&
			pParentNode->type == threeNode &&
			pMiddleNode->type == threeNode)
	{
		pLeftNode->smallKey=pParentNode->smallKey;
		pLeftNode->rightChild= pMiddleNode->leftChild;
		if (pLeftNode->rightChild!=NULL)
			pLeftNode->rightChild->parent=pLeftNode;
		pLeftNode->type=twoNode;
		pParentNode->smallKey=pMiddleNode->smallKey;
		pMiddleNode->smallKey=pMiddleNode->bigKey;
		pMiddleNode->bigKey=NULL;
		pMiddleNode->leftChild= pMiddleNode->middleChild;
		pMiddleNode->middleChild=NULL;
		pMiddleNode->type=twoNode;
		pLeftNode->updateNodeSize();
		pMiddleNode->updateNodeSize();
		fixed = true;
	}
	/* Case 3: Take from right. */
	else if(pNode == pLeftNode &&
			pParentNode->type == threeNode &&
			pRightNode->type == threeNode)
	{
		pLeftNode->smallKey=pParentNode->smallKey;
		pLeftNode->rightChild= pMiddleNode->leftChild;
		if (pLeftNode->rightChild!=NULL)
			pLeftNode->rightChild->parent=pLeftNode;
		pLeftNode->type=twoNode;
		pParentNode->smallKey=pMiddleNode->smallKey;
		pMiddleNode->smallKey=pParentNode->bigKey;
		pMiddleNode->leftChild= pMiddleNode->rightChild;

		pMiddleNode->rightChild= pRightNode->leftChild;
		if (pMiddleNode->rightChild!=NULL)
			pMiddleNode->rightChild->parent=pMiddleNode;
		pParentNode->bigKey=pRightNode->smallKey;
		pRightNode->smallKey=pRightNode->bigKey;
		pRightNode->bigKey=NULL;
		pRightNode->leftChild= pRightNode->middleChild;
		pRightNode->middleChild=NULL;
		pRightNode->type=twoNode;
		pLeftNode->updateNodeSize();
		pMiddleNode->updateNodeSize();
		pRightNode->updateNodeSize();
		fixed = true;
	}
	/* Case 4: Take from right first. */
	else if(pNode == pMiddleNode &&
			pRightNode->type == threeNode)
	{
		pMiddleNode->smallKey=pParentNode->bigKey;
		pMiddleNode->rightChild=pRightNode->leftChild;
		if (pMiddleNode->rightChild!=NULL)
			pMiddleNode->rightChild->parent=pMiddleNode;
		pMiddleNode->type=twoNode;
		pParentNode->bigKey=pRightNode->smallKey;
		pRightNode->smallKey=pRightNode->bigKey;
		pRightNode->bigKey=NULL;
		pRightNode->leftChild=pRightNode->middleChild;
		pRightNode->middleChild=NULL;
		pRightNode->type=twoNode;
		pMiddleNode->updateNodeSize();
		pRightNode->updateNodeSize();
		fixed = true;
	}
	/* Case 5: Take from left. */
	else if (pNode == pMiddleNode &&
			 pLeftNode->type == threeNode)
	{
		pMiddleNode->smallKey=pParentNode->smallKey;
		pMiddleNode->rightChild=pMiddleNode->leftChild;
		pMiddleNode->leftChild= pLeftNode->rightChild;
		if (pMiddleNode->leftChild!=NULL)
			pMiddleNode->leftChild->parent=pMiddleNode;
		pMiddleNode->type=twoNode;
		pParentNode->smallKey=pLeftNode->bigKey;
		pLeftNode->bigKey=NULL;
		pLeftNode->rightChild= pLeftNode->middleChild;
		pLeftNode->middleChild=NULL;
		pLeftNode->type=twoNode;
		pMiddleNode->updateNodeSize();
		pLeftNode->updateNodeSize();
		fixed = true;
	}
	/* Case 6: Take from left. */
	if(pNode == pRightNode &&
		pParentNode->type == twoNode &&
		pLeftNode->type == threeNode)
	{
		pRightNode->smallKey=pParentNode->smallKey;
		pRightNode->rightChild= pRightNode->leftChild;
		pRightNode->leftChild= pLeftNode->rightChild;
		if (pRightNode->leftChild!=NULL)
			pRightNode->leftChild->parent=pRightNode;
		pRightNode->type=twoNode;
		pParentNode->smallKey=pLeftNode->bigKey;
		pLeftNode->bigKey=NULL;
		pLeftNode->rightChild= pLeftNode->middleChild;
		pLeftNode->middleChild=NULL;
		pLeftNode->type=twoNode;
		pRightNode->updateNodeSize();
		pLeftNode->updateNodeSize();
		fixed = true;
	}
	/* Case 7: Take from middle. */
	else if (pNode == pRightNode &&
			 pParentNode->type == threeNode &&
			 pMiddleNode->type ==  threeNode)
	{
		pRightNode->smallKey=pParentNode->bigKey;
		pRightNode->rightChild=pRightNode->leftChild;
		pRightNode->leftChild= pMiddleNode->rightChild;
		if (pRightNode->leftChild!=NULL)
			pRightNode->leftChild->parent=pRightNode;
		pRightNode->type=twoNode;
		pParentNode->bigKey=pMiddleNode->bigKey;
		pMiddleNode->bigKey=NULL;
		pMiddleNode->rightChild= pMiddleNode->middleChild;
		pMiddleNode->middleChild=NULL;
		pMiddleNode->type=twoNode;
		pRightNode->updateNodeSize();
		pMiddleNode->updateNodeSize();
		fixed = true;
	}
	/* Case 8: Take from left. */
	else if (pNode == pRightNode &&
			 pParentNode->type == threeNode &&
			 pLeftNode->type ==  threeNode )
	{
		pRightNode->smallKey=pParentNode->bigKey;
		pRightNode->rightChild=pRightNode->leftChild;
		pRightNode->leftChild= pMiddleNode->rightChild;
		if (pRightNode->leftChild!=NULL)
			pRightNode->leftChild->parent=pRightNode;
		pRightNode->type=twoNode;
		pParentNode->bigKey=pMiddleNode->smallKey;
		pMiddleNode->smallKey=pParentNode->smallKey;
		pMiddleNode->rightChild= pMiddleNode->leftChild;
		pMiddleNode->leftChild=pLeftNode->rightChild;
		if (pMiddleNode->leftChild!=NULL)
			pMiddleNode->leftChild->parent=pMiddleNode;
		pParentNode->smallKey=pLeftNode->bigKey;
		pLeftNode->bigKey=NULL;
		pLeftNode->rightChild=pLeftNode->middleChild;
		pLeftNode->middleChild=NULL;
		pLeftNode->type=twoNode;
		pRightNode->updateNodeSize();
		pMiddleNode->updateNodeSize();
		pLeftNode->updateNodeSize();
		fixed = true;
	}
	return fixed;
}

bool BalancedTreeMTFList::tryLocalMerge(BalancedTreeNode *pNode){  //successful if the node has two siblings
	BalancedTreeNode* pParentNode = NULL;
	BalancedTreeNode* pRightNode = NULL;
	BalancedTreeNode* pMiddleNode = NULL;
	BalancedTreeNode* pLeftNode = NULL;
	bool fixed = false;

	pParentNode = pNode->parent;
	pRightNode = pParentNode->rightChild;
	pMiddleNode = pParentNode->middleChild;
	pLeftNode = pParentNode->leftChild;

	if(pNode == pLeftNode &&
	   pParentNode->type == threeNode) /* Case 15: Merge to the left. */
	{
		pLeftNode->smallKey=pParentNode->smallKey;
		pLeftNode->bigKey=pMiddleNode->smallKey;
		pLeftNode->middleChild= pMiddleNode->leftChild;
		if (pLeftNode->middleChild!=NULL)
			pLeftNode->middleChild->parent=pLeftNode;
		pLeftNode->rightChild= pMiddleNode->rightChild;
		if (pLeftNode->rightChild!=NULL)
			pLeftNode->rightChild->parent=pLeftNode;
		pLeftNode->type=threeNode;
		pParentNode->smallKey=pParentNode->bigKey;
		pParentNode->bigKey=NULL;
		pParentNode->middleChild=NULL;
		pParentNode->type=twoNode;

		pLeftNode->updateNodeSize();
		pParentNode->updateNodeSize();

		delete pMiddleNode;
		fixed = true;
	}
	else if(pNode == pMiddleNode &&
			pParentNode->type == threeNode) /* Case 16: Merge to the left. */
	{
		pLeftNode->bigKey=pParentNode->smallKey;
		pLeftNode->middleChild= pLeftNode->rightChild;
		pLeftNode->rightChild= pMiddleNode->leftChild;
		if (pLeftNode->rightChild!=NULL)
			pLeftNode->rightChild->parent=pLeftNode;
		pLeftNode->type=threeNode;
		pParentNode->smallKey=pParentNode->bigKey;
		pParentNode->bigKey=NULL;
		pParentNode->middleChild=NULL;
		pParentNode->type=twoNode;
		pLeftNode->updateNodeSize();
		pParentNode->updateNodeSize();

		delete pMiddleNode;
		fixed = true;
	}
	else if(pNode == pRightNode &&
			pParentNode->type == threeNode) /* Case 17: Merge to the right. */
	{
		pRightNode->smallKey=pMiddleNode->smallKey;
		pRightNode->bigKey=pParentNode->bigKey;
		pRightNode->rightChild= pRightNode->leftChild;
		pRightNode->middleChild= pMiddleNode->rightChild;
		if (pRightNode->middleChild!=NULL)
			pRightNode->middleChild->parent=pRightNode;
		pRightNode->leftChild= pMiddleNode->leftChild;
		if (pRightNode->leftChild!=NULL)
			pRightNode->leftChild->parent=pRightNode;
		pRightNode->type=threeNode;
		pParentNode->bigKey=NULL;
		pParentNode->middleChild=NULL;
		pParentNode->type=twoNode;

		pRightNode->updateNodeSize();
		pParentNode->updateNodeSize();

		delete pMiddleNode;
		fixed = true;
	}
	return fixed;
}

bool BalancedTreeMTFList::tryGlobalMerge(BalancedTreeNode *pNode,BalancedTreeNode **ppNextNodeToFix){
	BalancedTreeNode* pParentNode = pNode->parent;
	BalancedTreeNode* pRightNode = pParentNode->rightChild;
	BalancedTreeNode* pLeftNode = pParentNode->leftChild;

	bool fixed = false;
	*ppNextNodeToFix=NULL;
	/* Case 18: Node merger propagates. */
	if(pNode == pLeftNode){
		pLeftNode->smallKey=pParentNode->smallKey;
		pLeftNode->bigKey=pRightNode->smallKey;
		pLeftNode->middleChild=pRightNode->leftChild;
		if (pLeftNode->middleChild!=NULL)
			pLeftNode->middleChild->parent=pLeftNode;
		pLeftNode->rightChild=pRightNode->rightChild;
		if (pLeftNode->rightChild!=NULL)
			pLeftNode->rightChild->parent=pLeftNode;
		pLeftNode->type=threeNode;
		pParentNode->smallKey=NULL;
		pParentNode->rightChild=NULL;
		pParentNode->type=empty;
		delete pRightNode;
		*ppNextNodeToFix = pParentNode;

		pLeftNode->updateNodeSize();
		pParentNode->updateNodeSize();
		fixed=true;
	}
	/* Case 19: Node merger propagates. */
	else if(pNode == pRightNode){
		pLeftNode->bigKey=pParentNode->smallKey;
		pLeftNode->middleChild=pLeftNode->rightChild;
		if (pLeftNode->middleChild!=NULL)
			pLeftNode->middleChild->parent=pLeftNode;
		pLeftNode->rightChild= pRightNode->leftChild;
		if (pLeftNode->rightChild!=NULL)
			pLeftNode->rightChild->parent=pLeftNode;
		pLeftNode->type=threeNode;
		pParentNode->rightChild=NULL;
		pParentNode->smallKey=NULL;
		pParentNode->type=empty;
		delete pRightNode;
		*ppNextNodeToFix = pParentNode;

		pLeftNode->updateNodeSize();
		pParentNode->updateNodeSize();
		fixed=true;
	}
	return fixed;
}




BalancedTreeNode * BalancedTreeMTFList::swapWithSuccessor(BalancedTreeNode* node,
		Pair* pair) {
	BalancedTreeNode *successorLeaf=getSuccessorLeaf(node,pair);
	if (pair->timestamp == node->smallKey->timestamp){
		node->smallKey=successorLeaf->smallKey;
		successorLeaf->smallKey=pair;
	} else if (pair->timestamp == node->bigKey->timestamp){
		node->bigKey=successorLeaf->smallKey;
		successorLeaf->smallKey=pair;
	}
	return successorLeaf;

}

BalancedTreeNode* BalancedTreeMTFList::getSuccessorLeaf(BalancedTreeNode* node,
		Pair *pKey) {
	//reach the leftmost leaf of the edge immediately to the right of the key to delete
	BalancedTreeNode *successorNode;
	if (pKey->timestamp==node->smallKey->timestamp){
		if (node->type==twoNode)
			successorNode=node->rightChild;
		else
			successorNode=node->middleChild;
	} else if (pKey->timestamp==node->bigKey->timestamp)
		successorNode=node->rightChild;

	while (successorNode->leftChild !=NULL)
		successorNode=successorNode->leftChild;

	return successorNode;
}

BalancedTreeNode* BalancedTreeMTFList::searchInTree(int64_t timestamp,Pair **pair,int64_t *index) {
	BalancedTreeNode *foundNode=root->search(timestamp,index);
	if (foundNode!=NULL){
		if (foundNode->smallKey->timestamp==timestamp)
			*pair=foundNode->smallKey;
		else
			*pair=foundNode->bigKey;
	}
	return foundNode;
}





} /* namespace std */
