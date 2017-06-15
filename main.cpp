#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <sstream>
#include <string>

using namespace std;
using std::stringstream;
using std::string;

/////////////////////////////////////////////// BST STARTS ///////////////////////////////////////////////

template<class T>
class BinarySearchTree
{
	private:
		template<class U>
		class Node
		{
			friend class BinarySearchTree;
			U value;
			Node<U> *right;
			Node<U> *left;
			Node(U value, Node<U>* right = NULL, Node<U>* left = NULL)
			{
				this->value = value;
				this->right = right;
				this->left = left;
			}
		};
		Node<T>* root;
		void insert(Node<T> *&, T);
		void destroySubtree(Node<T> *);
		void remove (Node<T> *&, T);
		void performDeletion(Node<T> *&);
		void displayInOrder(Node<T> *);
		void displayInPreOrder(Node<T> *);
		void displayInPostOrder(Node<T> *);
	
	public:
		BinarySearchTree() { root = NULL; }
		~BinarySearchTree() { destroySubtree(root); }
		void insert(T value) { insert(root, value); }
		void remove(T value) { remove(root, value); }
		void makeDeletion(Node<T> *& node) { performDeletion(node); }
		void showInOrder() { displayInOrder(root); }
		void showInPostOrder() { displayInPostOrder(root); }
		void showInPreOrder() { displayInPreOrder(root); }
		bool search(T value);
};

/*
 * Insert a new node in to the tree. 
 */
template <class T>
void BinarySearchTree<T>::insert(Node<T> *&tree, T value)
{
	if(tree == NULL)						// Check to see if tree is null
	{
		tree = new Node<T>(value);			// If so, make the new node the root node
		return;
	}
	
	if(tree->value == value)				// Does the node already exist in the tree?
	{
		return;								// If so, return to prevent duplicate entries
	}

	
	if(tree->value > value)					// If the value is greater than a root node
	{
		insert(tree->left, value);			// Recursively call insert to the left subtree
	}
	else
	{
		insert(tree->right, value);			// Recursively call insert to the right subtree
	}
}

/*
 * Remove a specific node from the tree
 */
template <class T>
void BinarySearchTree<T>::remove(Node<T> *&tree, T value)
{
	if (tree == NULL)						// Make sure tree isn't empty
	{
		return;
	}

	if (tree->value > value)				// Find where the value is
	{
		remove(tree->left, value);
	}
	else if (tree->value < value)
	{
		remove(tree->right, value);
	}
	else
	{
		performDeletion(tree);				// Delete the node from the tree
	}
}

template <class T>
void BinarySearchTree<T>::performDeletion(Node<T> *&tree)
{

	Node<T> *nodeToRemove = tree;			// Hold node to be deleted
	Node<T> *nodeLink;						// used to find left subtree link attachment
	
	if (tree->right == NULL)				// if right child is null
	{				
		tree = tree->left;					// replace the tree with the left subtree
	}
	else if (tree->left == NULL)			// If left child is null
	{
		tree = tree->right;					// replace the tree with the right subtree
	}
	else
	{
		nodeLink = tree->right;				// We have two children from the root node, attach the right side of the tree
		while (nodeLink->left != NULL)		// Find the smallest node in the right-left subtree
		{
			nodeLink = nodeLink->left;		// assign left subtree to nodeLink
		}
		nodeLink->left = tree->left;		// attach left subtree of oringinal tree as left subtree of smallest node in right subtree
		tree = tree->right;					
	}
	delete nodeToRemove;					// Remove the node we want to delete
}

template <class T>
void BinarySearchTree<T>::displayInOrder(Node<T> *tree)
{
	if (tree)
	{
		displayInOrder(tree->left);
		std::cout << tree->value << " ";
		displayInOrder(tree->right);
	}
	std::cout << std::endl;
}

template <class T>
void BinarySearchTree<T>::displayInPostOrder(Node<T> *tree)
{
	if(tree)
	{
		cout << tree->value << " ";
		displayInPostOrder(tree->left);
		displayInPostOrder(tree->right);
	}
	std::cout << std::endl;
}

template <class T>
void BinarySearchTree<T>::displayInPreOrder(Node<T> *tree)
{
	if(tree)
	{
		displayInPreOrder(tree->left);
		displayInPreOrder(tree->right);
		cout << tree->value << " ";
	}
	std::cout << std::endl;
}


template <class T>
bool BinarySearchTree<T>::search(T value)
{
	Node<T>* tree = root;

	while (tree)
	{
		if (tree->value)
		{
			return true;
		}
		else if (tree->value > value)
		{
			tree = tree->left;
		}
		else
		{
			tree = tree->right;
		}
	}
	return false;
}


template <class T>
void BinarySearchTree<T>::destroySubtree(Node<T>* tree)
{
	if (!tree)
		return;
	destroySubtree(tree->left);
	destroySubtree(tree->right);
	delete tree;
}

/////////////////////////////////////////////// HEAP STARTS ///////////////////////////////////////////////

/*
 * Class Declaration
 */
template <class T>
class BinaryHeap
{
    private:
        vector <T> heap;
        int left(T parent);
        int right(T parent);
        int parent(T child);
        void heapifyup(int index);
        void heapifydown(int index);
    public:
        BinaryHeap()
        {}
        void Insert(T element);
        void DeleteMin();
        int ExtractMin();
        void DisplayHeap();
        int Size();
};
/*
 * Return Heap Size
 */
template <class T>
int BinaryHeap<T>::Size()
{
    return heap.size();
}
 
/*
 * Insert Element into a Heap
 */
template <class T>
void BinaryHeap<T>::Insert(T element)
{
    heap.push_back(element);
    heapifyup(heap.size() -1);
}
/*
 * Delete Minimum Element
 */
template <class T>
void BinaryHeap<T>::DeleteMin()
{
    if (heap.size() == 0)
    {
        cout<<"Heap is Empty"<<endl;
        return;
    }
    heap[0] = heap.at(heap.size() - 1);
    heap.pop_back();
    heapifydown(0);
    //cout<<"Element Deleted"<<endl;
}
 
/*
 * Extract Minimum Element
 */
 template <class T>
int BinaryHeap<T>::ExtractMin()
{
    if (heap.size() == 0)
    {
        return -1;
    }
    else
        return heap.front();
}
 
/*
 * Display Heap
 */
template <class T>
void BinaryHeap<T>::DisplayHeap()
{
    vector <int>::iterator pos = heap.begin();
    cout<<"Heap -->  ";
    while (pos != heap.end())
    {
        cout<<*pos<<" ";
        pos++;
    }
    cout<<endl;
}
 
/*
 * Return Left Child
 */
template <class T>
int BinaryHeap<T>::left(T parent)
{
    int l = 2 * parent + 1;
    if (l < heap.size())
        return l;
    else
        return -1;
}
 
/*
 * Return Right Child
 */
template <class T>
int BinaryHeap<T>::right(T parent)
{
    int r = 2 * parent + 2;
    if (r < heap.size())
        return r;
    else
        return -1;
}
 
/*
 * Return Parent
 */
template <class T>
int BinaryHeap<T>::parent(T child)
{
    int p = (child - 1)/2;
    if (child == 0)
        return -1;
    else
        return p;
}
 
/*
 * Heapify- Maintain Heap Structure bottom up
 */
template <class T>
void BinaryHeap<T>::heapifyup(int in)
{
    if (in >= 0 && parent(in) >= 0 && heap[parent(in)] > heap[in])
    {
        int temp = heap[in];
        heap[in] = heap[parent(in)];
        heap[parent(in)] = temp;
        heapifyup(parent(in));
    }
}
 
/*
 * Heapify- Maintain Heap Structure top down
 */
template <class T> 
void BinaryHeap<T>::heapifydown(int in)
{
 
    int child = left(in);
    int child1 = right(in);
    if (child >= 0 && child1 >= 0 && heap[child] > heap[child1])
    {
       child = child1;
    }
    if (child > 0 && heap[in] > heap[child])
    {
        int temp = heap[in];
        heap[in] = heap[child];
        heap[child] = temp;
        heapifydown(child);
    }
}

/////////////////////////////////////////////// RED-BLACK STARTS ///////////////////////////////////////////


#include <iostream>
#include <queue>

using namespace std;


enum Color{BLACK = 0, RED = 1};
template <class T>
struct rbtNode{rbtNode<T> *left, *right, *parent; T data; Color color;}; // 0 = black, 1 = red

template <class T>
class RBT{

public:
    RBT();
    ~RBT();

    void insert(T key);
    void bfs();
    void DrawTree();
    bool remove(T key);
    bool search(T key);
    bool isEmpty();


private:
    rbtNode<T> *root;
    int DrawTreePrivate(rbtNode<T> *tree, int is_left, int offset, int depth, char s[20][255]);

    /* helper functions */
    rbtNode<T> *createLeaf(T data);
    rbtNode<T> *findGrandparent(rbtNode<T>* node);
    rbtNode<T> *findUncle(rbtNode<T>* node);
    rbtNode<T> *findSibling(rbtNode<T>* node);
    rbtNode<T> *searchForNode(T key);
    rbtNode<T> *inorderPredecessor(rbtNode<T> *removenode);

    /* restructuring & recoloring methods */
    void recolor(rbtNode<T>* x);
    void fixDoubleblack(rbtNode<T>* current);
    void swapColors(rbtNode<T>* a, rbtNode<T> *b);
    void restructure(rbtNode<T> *x);
    void fixRoot();
    void fixTree(rbtNode<T> *start);
    void singleRotateRight(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent);
    void singleRotateLeft(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent);
    void rotateRightRight(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent);
    void rotateLeftLeft(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent);
    void rotateLeftRight(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent);
    void rotateRightLeft(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent);
    bool isLeaf(rbtNode<T> *current) const;
    bool NodeHasRedChild(rbtNode<T> *current) const;
    void removeLeaf(rbtNode<T> *leafnode);
    void removeRoot(rbtNode<T> *rootnode);

};

template <class T>
RBT<T>::RBT() { root = NULL; }

template <class T>
RBT<T>::~RBT() { }

template <class T>
void RBT<T>::fixRoot() {
    if(root->parent){
        while(root->parent)
            root = root->parent;
    }
}

template <class T>
void RBT<T>::insert(T key) {
    if(isEmpty()) {
        root = createLeaf(key);
        root->color = BLACK;
    }
    else{
        rbtNode<T> *curr{root}, *trail{NULL};
        while(curr){
            trail = curr;
            (key < curr->data) ? curr = curr->left : curr = curr->right;
        }

        // insert the node
        rbtNode<T> *newNode{NULL};
        if (key < trail->data){
            trail->left = createLeaf(key); trail->left->parent = trail;
            newNode = trail->left;
        }
        else{
            trail->right = createLeaf(key); trail->right->parent = trail;
            newNode = trail->right;
        }

        fixTree(newNode);
    }
}

template <class T>
void RBT<T>::fixTree(rbtNode<T> *start) {

    if(start->parent){
        if(start->parent->color == 1 && start->color == 1 )
        {
            // verify node has an uncle
            rbtNode<T>*uncle = findUncle(start);
            if(uncle){
                if(uncle->color == 1) // if the node's uncle is red, recolor the tree
                    recolor(start);
                else { // if the node's uncle is black, restructure the tree
                    restructure(start);
                    fixRoot();
                }
            }
            else { // NULL uncle nodes still require restructuring
                restructure(start);
                fixRoot();
            }
        }
    }
}

template <class T>
void RBT<T>::restructure(rbtNode<T> *x) {
    rbtNode<T> *grandparent = findGrandparent(x);
    if(grandparent){
        rbtNode<T> *parent = x->parent;
        if(grandparent->left == parent){
            if(parent->left == x) {
                rotateLeftLeft(x, parent, grandparent);
                fixTree(parent);
            }
            else {
                rotateLeftRight(x, parent, grandparent);
                fixTree(x); // come back here
            }
        }
        else{
            if(parent->right == x) {
                rotateRightRight(x, parent, grandparent);
                fixTree(parent);
            }
            else {
                rotateRightLeft(x, parent, grandparent);
                fixTree(x);
            }
        }
    }
}

template <class T>
void RBT<T>::singleRotateLeft(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent) {
    parent->right = x->left;
    if(x->left)
        x->left->parent = parent;
    x->left = parent;
    parent->parent = x;
    x->parent = grandparent;
    grandparent->left = x;
}

template <class T>
void RBT<T>::singleRotateRight(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent) {
    parent->left = x->right;
    if(x->right)
        x->right->parent = parent;
    x->right = parent;
    parent->parent = x;
    x->parent = grandparent;
    grandparent->right = x;
}

template <class T>
void RBT<T>::rotateLeftLeft(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent) {
    rbtNode<T> *greatGrandparent = grandparent->parent;
    grandparent->left = parent->right;
    if(parent->right)
        parent->right->parent = grandparent;
    grandparent->parent = parent;
    parent->right = grandparent;
    parent->parent = NULL;
    swapColors(parent, grandparent);

    if(greatGrandparent) {
        parent->parent = greatGrandparent;
        greatGrandparent->right == grandparent ?
                greatGrandparent->right = parent : greatGrandparent->left = parent;
    }
}

template <class T>
void RBT<T>::rotateRightRight(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent) {
    rbtNode<T> *greatGrandparent = grandparent->parent;
    grandparent->right = parent->left;
    if(parent->left)
        parent->left->parent = grandparent;
    grandparent->parent = parent;
    parent->left = grandparent;
    parent->parent = NULL;
    swapColors(parent, grandparent);

    if(greatGrandparent) {
        parent->parent = greatGrandparent;
        greatGrandparent->right == grandparent ?
                greatGrandparent->right = parent : greatGrandparent->left = parent;
    }
}

template <class T>
void RBT<T>::rotateLeftRight(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent){
    singleRotateLeft(x,parent,grandparent);
    rotateLeftLeft(parent,x,grandparent);
}

template <class T>
void RBT<T>::rotateRightLeft(rbtNode<T> *x, rbtNode<T> *parent, rbtNode<T> *grandparent){
    singleRotateRight(x,parent,grandparent);
    rotateRightRight(parent,x,grandparent);
}


template <class T>
void RBT<T>::recolor(rbtNode<T> *x) {
        x->parent->color = BLACK;
        rbtNode<T> *uncle = findUncle(x);
        if(uncle)
            uncle->color = BLACK;
        rbtNode<T> *grandparent = findGrandparent(x);
        if(grandparent) {
            grandparent->color = RED;
            if(root == grandparent)
                grandparent->color = BLACK;
        }
        fixTree(grandparent); // move up tree
}

template <class T>
void RBT<T>::swapColors(rbtNode<T> *a, rbtNode<T> *b) {
    Color temp = a->color; a->color = b->color; b->color = temp;
}


template <class T>
bool RBT<T>::remove(T key) {
    rbtNode<T> *removenode = searchForNode(key);

    if(!removenode)
        return false;
    else if(isLeaf(removenode))
        removeLeaf(removenode);
    else
        removeRoot(removenode);
    return true;
}

template <class T>
void RBT<T>::fixDoubleblack(rbtNode<T>* current){
    rbtNode<T>*sibling = findSibling(current);
    if(sibling){
        if(sibling->color == BLACK){
            if(NodeHasRedChild(sibling)){

            }
            else{ // node must have black children

            }
        }
        else{ // sibling must be red

        }
    }
    else{
        // treat as if node has a black sibling
    }
}

template <class T>
bool RBT<T>::NodeHasRedChild(rbtNode<T> *current) const{
    if(!current)
        return false;
    if(current->left)
        if(current->left->color == RED)
            return true;
    if(current->right)
        if(current->right->color == RED)
            return true;
    return false;
}


template <class T>
bool RBT<T>::search(T key) {
    rbtNode<T> *searchNode = searchForNode(key);
    if(searchNode)
        return true;
    return false;
}

template <class T>
bool RBT<T>::isEmpty() {
    if(!root)
        return true;
    return false;
}

template <class T>
rbtNode<T>*RBT<T>::createLeaf(T data) {
    rbtNode<T>* newNode = new rbtNode<T>;
    newNode->data   = data;
    newNode->left   = NULL;
    newNode->right  = NULL;
    newNode->parent = NULL;
    newNode->color  = RED;
    return newNode;
}

template <class T>
rbtNode<T>*RBT<T>::searchForNode(T key) {
    if(isEmpty())
        return NULL;
    rbtNode<T> *curr = root;
    bool found = false;
    while(curr && !found){
        if(key > curr->data)
            curr = curr->right;
        else if(key < curr->data)
            curr = curr->left;
        else
            found = true;
    }
    if(found)
        return curr;
    return NULL;
}

// finds grandparent of a node
template <class T>
rbtNode<T>* RBT<T>::findGrandparent(rbtNode<T> *node) {
    if(node->parent)
            return node->parent->parent;
    return NULL;
}

// finds uncle of a node
template <class T>
rbtNode<T>* RBT<T>::findUncle(rbtNode<T> *node) {
    rbtNode<T> *g = findGrandparent(node);
    if(!g)
        return NULL;
    if(node->parent == g->left)
        return g->right;
    else
        return g->left;
}

template <class T>
rbtNode<T>* RBT<T>::findSibling(rbtNode<T>* node){
    if(!node || !node->parent)
        return NULL;
    if(node->parent->left == node)
        return node->parent->right;
    else
        return node->parent->left;
}

template <class T>
void RBT<T>::bfs() {
    if(!isEmpty()) {
        std::queue<rbtNode<T> *> queue;
        queue.push(root);
        while (!queue.empty()) {
            rbtNode<T> *temp = queue.front();
            queue.pop();

            if (temp->left)
                queue.push(temp->left);
            if (temp->right)
                queue.push(temp->right);

            std::cout << temp->data;
            temp->color == BLACK ? std::cout << "(black) " : std::cout << "(red) ";
        }
    }
}


// Drawing functions
template <class T>
void RBT<T>::DrawTree() {
    char s[20][255];
    for (int i = 0; i < 20; i++)
        sprintf(s[i], "%80s", " ");

    DrawTreePrivate(root, 0, 0, 0, s);

    for (int i = 0; i < 20; i++)
        printf("%s\n", s[i]);

}

template <class T>
int RBT<T>::DrawTreePrivate(rbtNode<T> *tree, int is_left, int offset, int depth, char s[20][255])
{
    char b[20];
    int width = 5;

    if (!tree) return 0;

    sprintf(b, "(%03d)", tree->color);

    int left  = DrawTreePrivate(tree->left,  1, offset,                depth + 1, s);
    int right = DrawTreePrivate(tree->right, 0, offset + left + width, depth + 1, s);


    for (int i = 0; i < width; i++)
        s[2 * depth][offset + left + i] = b[i];

    if (depth && is_left) {

        for (int i = 0; i < width + right; i++)
            s[2 * depth - 1][offset + left + width/2 + i] = '-';

        s[2 * depth - 1][offset + left + width/2] = '+';
        s[2 * depth - 1][offset + left + width + right + width/2] = '+';

    } else if (depth && !is_left) {

        for (int i = 0; i < left + width; i++)
            s[2 * depth - 1][offset - width/2 + i] = '-';

        s[2 * depth - 1][offset + left + width/2] = '+';
        s[2 * depth - 1][offset - width/2 - 1] = '+';
    }

    return left + width + right;
}

template <class T>
rbtNode<T> *RBT<T>::inorderPredecessor(rbtNode<T> *removenode) {
    if(isEmpty())
        return NULL;
    rbtNode<T> *curr = removenode->left;
    if(!curr)
        return NULL;

    while(curr->right){
        curr = curr->right;
    }
    return curr;
}

template <class T>
bool RBT<T>::isLeaf(rbtNode<T> *current) const{
    if(!current->left && !current->right)
        return true;
    return false;
}

template <class T>
void RBT<T>::removeRoot(rbtNode<T> *rootnode) {
    // if the node only has a right child only ** needs fixing **
    if(rootnode->right && !rootnode->left){
        rbtNode<T> *parent = rootnode->parent;

        /* If the node is red, delete node like a normal BST.
         * Else if the node is black, correct double blackness.*/
        if(rootnode->color == RED){
            if(parent){
                parent->right == rootnode ? parent->right = rootnode->right : parent->left = rootnode->right;
                rootnode->right->parent = parent;
                rootnode->right = NULL;
            }
            else{
                root = root->right;
                root->parent = NULL;
                rootnode->right = NULL;
            }
        }
        else{ // Node with only 1 right child is black
            fixDoubleblack(rootnode);
        }

        rootnode->parent = NULL;
        delete rootnode;
    }

        // if the node has a left child
    else {
        rbtNode<T> *replacement = inorderPredecessor(rootnode);
        rootnode->data = replacement->data;

        // normal BST delete
        if(replacement->color == RED) {
            if (isLeaf(replacement))
                removeLeaf(replacement);
            else {
                rbtNode<T> *parent = replacement->parent;
                parent->right == replacement ? parent->right = replacement->left : parent->left = replacement->left;
                replacement->left->parent = parent;
                replacement->left = NULL;
                replacement->parent = NULL;
            }
        }
        else
            fixDoubleblack(replacement);

        delete replacement;
    }
}

template <class T>
void RBT<T>::removeLeaf(rbtNode<T> *leafnode) {
    if(leafnode->parent){
        rbtNode<T> *parent = leafnode->parent;
        if(parent->left == leafnode)
            parent->left = NULL;
        else
            parent->right = NULL;

        /* If leaf is red, simply disconnect from tree and delete.
         * If leaf is black, fix blackness then delete. */
        if(leafnode->color == RED)
            leafnode->parent = NULL;
        else
            fixDoubleblack(leafnode);
        delete leafnode;
    }
    else{
        delete root;
        root = NULL;
    }
}

/////////////////////////////////////////////// AVL STARTS ///////////////////////////////////////////////

template<class T>
class AVLTree;

//The maximum imbalance that is allowed in the AVL Tree.
const static int MAX_IMBAL = 2;

/* AVL Node class, holds data and references to adjacent nodes. */
template<class T>
class AVLNode {

	public:
		friend class AVLTree<T>;
		AVLNode();
		AVLNode(const T & newData);
		AVLNode(const T & element, int theHeight);
		~AVLNode();
	private:
		AVLNode<T> * left;
		AVLNode<T> * right;
		const T & data;
		int leftHeight;
		int rightHeight;
		int height();
};

/* AVL Tree class, holds all of the functions etc for the AVL Tree. */
template<class T>
class AVLTree {

	public:
		friend class AVLNode<T>;
		AVLTree();
		AVLTree(const T & element);
		~AVLTree();
		void insert(const T & element);
		bool remove(const T & element);
		const T * find(const T & element);
		const T * findMax();
		const T * findMin();
		int size();
		const string inOrderTraversal();
		void clear();

	private:
		AVLNode<T> * root;
		int treeSize;

		void delete_traversal(AVLNode<T> * subRoot);
		AVLNode<T> * insert(AVLNode<T> * subRoot, const T & element);
		AVLNode<T> * remove(AVLNode<T> * subRoot, const T & element, bool * val);
		AVLNode<T> * balance(AVLNode<T> * subRoot);
		AVLNode<T> * singleRotateRightChild(AVLNode<T> * subRoot);
		AVLNode<T> * singleRotateLeftChild(AVLNode<T> * subRoot);
		const T * find(AVLNode<T> * subRoot, const T & element);
		AVLNode<T> * findMin(AVLNode<T> * subRoot);
		AVLNode<T> * findMax(AVLNode<T> * subRoot);
		void inOrderTraversal(AVLNode<T> * subRoot, string & output, int depth);
};

//~AVLTree Functions--------------------------------------------------------------------
/* AVLTree constructor, sets root to NULL. */
template<class T>
AVLTree<T>::AVLTree() {

	root = NULL;
	treeSize = 0;
}

/* AVLTree constructor, sets root to element. */
template<class T>
AVLTree<T>::AVLTree(const T & element) {

	root = new AVLNode<T>(element);
	treeSize = 1;
}

/* AVLTree destructor, deletes all of the nodes left in the tree. */
template<class T>
AVLTree<T>::~AVLTree() {

	delete_traversal(root);
}

/* Deletes all the nodes in the subTree that has subRoot as the root. */
template<class T>
void AVLTree<T>::delete_traversal(AVLNode<T> * subRoot) {

	if (!subRoot)
		return;

	delete_traversal(subRoot->left);
	delete_traversal(subRoot->right);
	delete subRoot;
}

/* Inserts the passed element into the AVLTree. */
template<class T>
void AVLTree<T>::insert(const T& element) {

	root = insert(root, element);
	treeSize++;
}

/* Helper method for inserting the passed element into the AVLTree. */
template<class T>
AVLNode<T> * AVLTree<T>::insert(AVLNode<T> * subRoot, const T & element) {

	AVLNode<T> * returnValue = subRoot;

	if (!subRoot) {
	
		returnValue = new AVLNode<T>(element);
	}
	else if (element >= subRoot->data) {
	
		subRoot->right = insert(subRoot->right, element);
		subRoot->rightHeight++;
	}
	else {
	
		subRoot->left = insert(subRoot->left, element);
		subRoot->leftHeight++;
	}

	return balance(returnValue);
}

/* Removes the passed element from the AVLTree (if it is present in the AVLTree). */
template<class T>
bool AVLTree<T>::remove(const T & element) {

	if (!root) {
	
		return false;
	}
	else {

		bool val = true;
		root = remove(root, element, &val);
		return val;
	}
}

/* Helper method to remove the passed element from the AVLTree (if it is present). */
template<class T>
AVLNode<T> * AVLTree<T>::remove(AVLNode<T> * subRoot, const T & element, bool * val) {

	AVLNode<T> * returnValue;

	if (!subRoot) {

		*val = false;
		returnValue = NULL;
	}
	else if (subRoot->data == element){

		treeSize--;

		//No children of the found node, remove it.
		if (!subRoot->left && !subRoot->right) {

			delete subRoot;
			returnValue = NULL;
		}
		//Two children of the found node. Remove it, complexly.
		else if (subRoot->left && subRoot->right) {

			AVLNode<T> * maxLeft = findMax(subRoot->left);
			AVLNode<T> * newSubRoot = new AVLNode<T>(maxLeft->data);
			newSubRoot->left = remove(subRoot->left, maxLeft->data, val);
			newSubRoot->right = subRoot->right;
			newSubRoot->leftHeight = subRoot->leftHeight;
			newSubRoot->rightHeight = subRoot->rightHeight;
			//Subtract from leftHeight to account for removal call above.
			newSubRoot->leftHeight--;
			//Increase the treeSize to make up for the remove called above.
			treeSize++;

			delete subRoot;
			returnValue = newSubRoot;
		}
		//One child of the found node. Remove it, and promote child.
		else {

			AVLNode<T> * temp = (subRoot->right) ? subRoot->right : subRoot->left;

			delete subRoot;
			returnValue = temp;
		}
	}
	else if (element >= subRoot->data) {

		subRoot->right = remove(subRoot->right, element, val);
		//If *val is true, then removal succeeded. Subtract from this node's rightHeight.
		if (*val)
			subRoot->rightHeight--;
		returnValue = subRoot;
	}
	else {

		subRoot->left = remove(subRoot->left, element, val);
		//If *val is true, then removal succeeded. Subtract from this node's leftHeight.
		if (*val)
			subRoot->leftHeight--;
		returnValue = subRoot;
	}

	return balance(returnValue);
}

/* Check the balance of the AVLTree that has subRoot as its root. */
template<class T>
AVLNode<T> * AVLTree<T>::balance(AVLNode<T> * subRoot) {

	AVLNode<T> * returnValue = subRoot;

	if (!subRoot) {

		returnValue = subRoot;
	}
	//If the left subTree is too large.
	else if (subRoot->leftHeight - subRoot->rightHeight > MAX_IMBAL) {

		//Single rotation is necessary if the left-subtree is larger than the right.
		if(subRoot->left->leftHeight >= subRoot->left->rightHeight) {

			returnValue = singleRotateLeftChild(subRoot);
		}
		//Double rotation is necessary if the right-subtree is larger than the left.
		else {

			subRoot->left = singleRotateRightChild(subRoot->left);
			returnValue = singleRotateLeftChild(subRoot);
		}
	}
	//If the right subTree is too large.
	else if (subRoot->rightHeight - subRoot->leftHeight > MAX_IMBAL) {

		//Single rotation is necessary if the right-subtree is larger than the left.
		if (subRoot->right->rightHeight >= subRoot->right->leftHeight) {

			returnValue = singleRotateRightChild(subRoot);
		}
		//Double rotation is necessary if the left sub-tree is larger than the right.
		else {

			subRoot->right = singleRotateLeftChild(subRoot->right);
			returnValue = singleRotateRightChild(subRoot);
		}
	}

	return returnValue;
}

/* Rotate the left subTree over to the root. */
template<class T>
AVLNode<T> * AVLTree<T>::singleRotateLeftChild(AVLNode<T> * subRoot) {

	AVLNode<T> * temp = subRoot->left;
	subRoot->left = temp->right;
	subRoot->leftHeight = temp->rightHeight;
	temp->right = subRoot;
	temp->rightHeight = subRoot->height() + 1;

	return temp;
}

/* Rotate the right subTree over to the root. */
template<class T>
AVLNode<T> * AVLTree<T>::singleRotateRightChild(AVLNode<T> * subRoot) {

	AVLNode<T> * temp = subRoot->right;
	subRoot->right = temp->left;
	subRoot->rightHeight = temp->leftHeight;
	temp->left = subRoot;
	temp->leftHeight = subRoot->height() + 1;

	return temp;
}

/* Searches the AVLTree for the passed in element.
 * Returns NULL if the element cannot be found.
 */
template<class T>
const T * AVLTree<T>::find(const T & element) {

	return find(root, element);
}

/* Helper method for finding the passed in element in the subTree with subRoot as the root. */
template<class T>
const T * AVLTree<T>::find(AVLNode<T> * subRoot, const T & element) {

	if (!subRoot)
		return NULL;
	else if (subRoot->data == element)
		return &subRoot->data;
	else if (subRoot->data > element)
		return find(subRoot->left, element);
	else 
		return find(subRoot->right, element);
}

/* Finds the minimum element in the tree. */
template<class T>
const T * AVLTree<T>::findMin() {

	if (!root)
		return NULL;
	return findMin(root)->data;
}

/* Helper function for finding the minimum element in the tree. */
template<class T>
AVLNode<T> * AVLTree<T>::findMin(AVLNode<T> * subRoot) {

	if (subRoot->left)
		return findMin(subRoot->left);
	else
		return subRoot;
}

/* Finds the maximum element in the AVLTree. */
template<class T>
const T * AVLTree<T>::findMax() {

	if (!root)
		return NULL;
	return findMax(root)->data;
}

/* Helper function for finding the maximum element in the tree. */
template<class T>
AVLNode<T> * AVLTree<T>::findMax(AVLNode<T> * subRoot) {

	if (subRoot->right)
		return findMax(subRoot->right);
	else
		return subRoot;
}

/* Performs an inOrderTraversal of the AVLTree. */
template<class T>
const string AVLTree<T>::inOrderTraversal() {

	string output;

	if (root)
		inOrderTraversal(root, output, 0);
	else
		output += "--|";

	return output;
}

/* Helper function for inOrderTraversal. Performs an inOrderTraversal of the subTree with subRoot as its root. */
template<class T>
void AVLTree<T>::inOrderTraversal(AVLNode<T> * subRoot, string & output, int depth) {

	stringstream ss;

	if (!subRoot)
		return;

	if (subRoot->left)
		inOrderTraversal(subRoot->right, output, depth + 1);

	for (int i = 0; i < depth; i++)
		output += "----";
	output += "|";
	ss << subRoot->data;// << "|||leftHeight = " << subRoot->leftHeight << "|||rightHeight = " << subRoot->rightHeight << "|||addr: " << subRoot;
	output += ss.str();
	output += "\n";
	
	if (subRoot->right)
		inOrderTraversal(subRoot->left, output, depth + 1);
}

/* Gets the size of the AVLTree. */
template<class T>
int AVLTree<T>::size() {

	return treeSize;
}

/* Empties out the AVLTree. */
template<class T>
void AVLTree<T>::clear() {

	delete_traversal(root);
	root = NULL;
	treeSize = 0;
}

//~AVLNode functions---------------------------------------------------------------------
/* Constructor for AVLNode, sets the node's data to element. */
template<class T>
AVLNode<T>::AVLNode(const T & element) : data(element) {

	leftHeight = 0;
	rightHeight = 0;
	left = NULL;
	right = NULL;
}

/* Destructor for AVLNode. */
template<class T>
AVLNode<T>::~AVLNode() {

	left = NULL;
	right = NULL;
}

/* Gets the total height of the node. Adds together leftHeight & rightHeight. */
template<class T>
int AVLNode<T>::height() {

	return leftHeight + rightHeight;
}


/////////////////////////////////////////////// MAIN STARTS ///////////////////////////////////////////////

int main () {
  int number;
  vector<int> insert_numbers;
  ifstream myfile ("insertion_data-500000.txt");
  //Taking insertion numbers to an integer array from txt file 
  if (myfile.is_open())
  {
    while ( myfile >> number )
    {
      insert_numbers.push_back(number);
    }
    myfile.close();
  }
  else cout << "Unable to open file";
  
  vector<int> delete_numbers;
  ifstream myfile2 ("deletion_data-50000.txt");
  //Taking deletion numbers to an integer array from txt file 
  if (myfile2.is_open())
  {
    while ( myfile2 >> number )
    {
      delete_numbers.push_back(number);
    }
    myfile2.close();
  }
  else cout << "Unable to open file"; 
  
  vector<int> search_numbers;
  ifstream myfile3 ("search_data-50000.txt");
  //Taking searching numbers to an integer array from txt file 
  if (myfile3.is_open())
  {
    while ( myfile3 >> number )
    {
      search_numbers.push_back(number);
    }
    myfile3.close();
  }
  else cout << "Unable to open file"; 
  

  clock_t start;
  double duration;
  
  cout<<"INSERT :"<<endl;
  
  start = clock();
  	BinarySearchTree<int> bst;
  	for (int i = 0; i < insert_numbers.size(); i++)	//BST Insertion
  	{
		bst.insert(insert_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Binary Search Tree Insertion Time = "<< duration <<" sec"<<'\n';  
  
  start = clock(); 
  BinaryHeap<int> bheap;
  	for (int i = 0; i < insert_numbers.size(); i++)	//Heap Insertion
  	{
		bheap.Insert(insert_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Binary Heap Insertion Time = "<< duration <<" sec"<<'\n';
  
  start = clock();
  RBT<int> rbt;
  	for (int i = 0; i < insert_numbers.size(); i++)	//Red-Black Tree Insertion
  	{
		rbt.insert(insert_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Red-Black Tree Insertion Time = "<< duration <<" sec"<<'\n';
  
  start = clock();
  AVLTree<int> avl;
  	for (int i = 0; i < insert_numbers.size(); i++)	//AVL Tree Insertion
  	{
		avl.insert(insert_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"AVL Tree Insertion Time = "<< duration <<" sec"<<'\n';
  
  cout<<"SEARCH :"<<endl;
  
  start = clock();
  	for (int i = 0; i < search_numbers.size(); i++)	//BST Search
  	{
		bst.search(search_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Binary Search Tree Search Time = "<< duration <<" sec"<<'\n';
  
  cout<<"Binary Heap Search Time = N/A"<<'\n';
  
  start = clock();
  	for (int i = 0; i < search_numbers.size(); i++)	//Red-Black Tree Search
  	{
		rbt.search(search_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Red-Black Tree Search Time = "<< duration <<" sec"<<'\n';
  
  start = clock();
  	for (int i = 0; i < search_numbers.size(); i++)	//AVL Tree Search
  	{
		avl.find(search_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"AVL Tree Search Time = "<< duration <<" sec"<<'\n';
  
  cout<<"DELETE :"<<endl;
  
  start = clock();
  	for (int i = 0; i < delete_numbers.size(); i++)	//BST Delete
  	{
		bst.remove(delete_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Binary Search Tree Delete Time = "<< duration <<" sec"<<'\n';
  
  start = clock();
  	for (int i = 0; i < delete_numbers.size(); i++)	//Heap Delete
  	{
		bheap.DeleteMin();
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Binary Heap Delete(DeleteMin) Time = "<< duration <<" sec"<<'\n';
  
  start = clock();
  	for (int i = 0; i < delete_numbers.size(); i++)	//Red-Black Tree Delete
  	{
		rbt.remove(delete_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"Red-Black Tree Delete Time = "<< duration <<" sec"<<'\n';
  
  start = clock();
  	for (int i = 0; i < delete_numbers.size(); i++)	//AVL Tree Delete
  	{
		avl.remove(delete_numbers[i]);
  	}
  duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
  cout<<"AVL Tree Delete Time = "<< duration <<" sec"<<'\n';
  return 0;
}
