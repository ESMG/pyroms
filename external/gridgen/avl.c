/******************************************************************************
 *
 * File:        avl.c
 *
 * Created      03/01/89
 *
 * Author:      Brad Appleton
 *
 * Description: Implementation of an AVL tree.
 *
 * Modifications:
 *
 *   Fri Jul 14 1989, Rev 1.0, brad(0165)
 *
 *   Wed Oct  8 14:11:07 EST 2003 -- modified by Pavel Sakov.
 *     Downloaded original code from 
 *     http://www.cmcrossroads.com/bradapp/ftp/src/libs/C++/libavl.tar.gz
 *     Reshuffled style and names.
 *     Modified to hold a shallow copy of data only.
 *     Simplified walk action format.
 *     Fixed avlfind() (node type should be updated in the process).
 *     Added avltree_findnextitem() and avltree_findprevitem().
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include "avl.h"                /* public types for avl trees */

typedef enum { LEFT = 0, RIGHT = 1 } DIRECTION;

#define OPPOSITE(x) (1 - (x))

/* return codes used by avlinsert(), avldelete(), and avlbalance() */
#define  HEIGHT_UNCHANGED 0
#define  HEIGHT_CHANGED   1

/*
 * IS_TREE     --  both subtrees are non-empty
 * IS_LBRANCH  --  left subtree is non-empty; right is empty
 * IS_RBRANCH  --  right subtree is non-empty; left is empty
 * IS_LEAF     --  both subtrees are empty
 * IS_NULL     --  given tree is empty
 */
typedef enum { IS_TREE, IS_LBRANCH, IS_RBRANCH, IS_LEAF, IS_NULL } NODETYPE;

typedef int (*avlcompare) ( /* void* item1, void* item2, NODETYPE node */ );

struct avlnode;
typedef struct avlnode avlnode;

struct avlnode {
    void* item;                 /* pointer to item */
    int bal;                    /* balance factor */
    avlnode* subtree[2];        /* LEFT and RIGHT subtrees */
};

struct avltree {
    avlnode* root;              /* pointer to the root node of the tree */
    avl_compare compare;        /* function used to compare keys */
    int count;                  /* number of nodes in the tree */
};

/* Constructor. Creates a new node.
 */
static avlnode* avlnode_create(void* item)
{
    avlnode* node;

    node = (avlnode*) malloc(sizeof(avlnode));
    node->item = item;          /* shallow copy only */
    node->bal = 0;
    node->subtree[LEFT] = node->subtree[RIGHT] = NULL;

    return node;
}

/* Destructor.
 * Resets the node pointer to NULL. Does not affect the data.
 */
static void avlnode_destroy(avlnode** node)
{
    free(*node);
    *node = NULL;
}

/* Determines the node type.
 */
static NODETYPE node_gettype(avlnode* tree)
{
    if (tree == NULL)
        return IS_NULL;
    else if ((tree->subtree[LEFT] != NULL) && (tree->subtree[RIGHT] != NULL))
        return IS_TREE;
    else if (tree->subtree[LEFT] != NULL)
        return IS_LBRANCH;
    else if (tree->subtree[RIGHT] != NULL)
        return IS_RBRANCH;
    else
        return IS_LEAF;
}

/* Compare function used to find the minimal element in a tree.
 */
static int avlmin(void* item1, void* item2, NODETYPE ntype)
{
    if (ntype == IS_RBRANCH || ntype == IS_LEAF)
        return 0;               /* left subtree is empty -- this is the
                                 * minimum */
    else
        return -1;              /* keep going left */
}

/* Compare function used to find the maximal element in a tree.
 */
static int avlmax(void* item1, void* item2, NODETYPE ntype)
{
    if (ntype == IS_LBRANCH || ntype == IS_LEAF)
        return 0;               /* right subtree is empty -- this is the
                                 * maximum */
    else
        return 1;               /* keep going right */
}

/* Rotates a node in specified direction to restore the balance of a tree.
 */
static int avlrotate1(avlnode** rootp, DIRECTION dir)
{
    DIRECTION other_dir = OPPOSITE(dir);
    avlnode* old_root = *rootp;
    int ht_unchanged = (*rootp)->subtree[other_dir]->bal == 0;;

    /*
     * assign new root 
     */
    *rootp = old_root->subtree[other_dir];

    /*
     * new-root exchanges it's "dir" subtree for it's parent 
     */
    old_root->subtree[other_dir] = (*rootp)->subtree[dir];
    (*rootp)->subtree[dir] = old_root;

    /*
     * update balances 
     */
    if (dir == LEFT)
        --((*rootp)->bal);
    else
        ++((*rootp)->bal);
    old_root->bal = -(*rootp)->bal;

    return ht_unchanged;
}

/* Rotates a node in specified direction and then in the opposite direction
 * to restore the balance of a tree.
 */
static void avlrotate2(avlnode** rootp, DIRECTION dir)
{
    DIRECTION other_dir = OPPOSITE(dir);
    avlnode* old_root = *rootp;
    avlnode* old_other_dir_subtree = (*rootp)->subtree[other_dir];

    /*
     * assign new root 
     */
    *rootp = old_root->subtree[other_dir]->subtree[dir];

    /*
     * new-root exchanges it's "dir" subtree for it's grandparent 
     */
    old_root->subtree[other_dir] = (*rootp)->subtree[dir];
    (*rootp)->subtree[dir] = old_root;

    /*
     * new-root exchanges it's "other-dir" subtree for it's parent 
     */
    old_other_dir_subtree->subtree[dir] = (*rootp)->subtree[other_dir];
    (*rootp)->subtree[other_dir] = old_other_dir_subtree;

    /*
     * update balances 
     */
    {
        int bal = (*rootp)->bal;

        (*rootp)->subtree[LEFT]->bal = (bal > 0) ? -bal : 0;
        (*rootp)->subtree[RIGHT]->bal = (bal > 0) ? 0 : -bal;
        (*rootp)->bal = 0;
    }
}

/* Determines and performs the  sequence of rotations needed to restore the
 * balance of the tree.
 * Returns 1 if tree height changed due to rotation; 0 otherwise.
 */
static int avlbalance(avlnode** rootp)
{
    int special_case = 0;

    if ((*rootp)->bal < -1) {   /* needs a right rotation */
        if ((*rootp)->subtree[LEFT]->bal == 1)
            avlrotate2(rootp, RIGHT);   /* needs a double RL rotation */
        else                    /* needs a single RR rotation */
            special_case = avlrotate1(rootp, RIGHT);
    } else if ((*rootp)->bal > 1) {     /* needs a left rotation */
        if ((*rootp)->subtree[RIGHT]->bal == -1)
            avlrotate2(rootp, LEFT);    /* needs a double LR rotation */
        else                    /* needs a single LL rotation */
            special_case = avlrotate1(rootp, LEFT);
    } else
        return HEIGHT_UNCHANGED;        /* no rotation made */

    return (special_case) ? HEIGHT_UNCHANGED : HEIGHT_CHANGED;
}

/* Finds an item in the tree.
 *
 * @param item Image of the item to find
 * @param tree A pointer to AVL tree
 * @param compare The compare function
 * @return A pointer to the item or NULL if not found.
 */
static void* avlfind(avlnode* tree, void* item, avlcompare compare)
{
    NODETYPE ntype = node_gettype(tree);
    int cmp;

    while (tree != NULL && (cmp = compare(item, tree->item, ntype)) != 0) {
        tree = tree->subtree[(cmp < 0) ? LEFT : RIGHT];
        ntype = node_gettype(tree);
    }

    return (tree == NULL) ? NULL : tree->item;
}

/* Inserts an item into the given tree.
 *
 * @param rootp A pointer to AVL tree
 * @param item A pointer to a pointer to the item to insert.
 * @param compare The compare function
 * @return Whether the tree height has changed (?).
 *
 * On return, *item is NULL if insertion succeeded, otherwise unchanged
 */
static int avlinsert(avlnode** rootp, void** item, avlcompare compare)
{
    int increase;
    int cmp;

    if (*rootp == NULL) {       /* insert new node here */
        *rootp = avlnode_create(*item);
        *item = NULL;           /* set return value in item */
        return HEIGHT_CHANGED;
    }

    cmp = compare(*item, (*rootp)->item, 0);

    if (cmp < 0) {              /* insert into the left subtree */
        increase = -avlinsert(&(*rootp)->subtree[LEFT], item, compare);
        if (*item != NULL)
            return HEIGHT_UNCHANGED;
    } else if (cmp > 0) {       /* insert into the right subtree */
        increase = avlinsert(&(*rootp)->subtree[RIGHT], item, compare);
        if (*item != NULL)
            return HEIGHT_UNCHANGED;
    } else {                    /* item already exists */
        *item = (*rootp)->item; /* set return value in item */
        return HEIGHT_UNCHANGED;
    }

    (*rootp)->bal += increase;  /* update balance factor */

    /*
     * re-balance if needed -- height of current tree increases only if its
     * subtree height increases and the current tree needs no rotation
     */
    if (increase && (*rootp)->bal)
        return 1 - avlbalance(rootp);

    return HEIGHT_UNCHANGED;
}

/* Delete an item from the tree.
 *
 * @param rootp A pointer to AVL tree
 * @param item A pointer to a pointer to the item to insert.
 * @Param size Item size in bytes.
 * @param compare The compare function
 * @return Whether the tree height has changed (?).
 *
 * On return, *item points to the deleted item item (or NULL if deletion
 * failed).
 */
static int avldelete(avlnode** rootp, void** item, avlcompare compare)
{
    int decrease = 0;
    int cmp;
    avlnode* old_root = *rootp;
    NODETYPE ntype = node_gettype(*rootp);
    DIRECTION dir = (ntype == IS_LBRANCH) ? LEFT : RIGHT;

    if (*rootp == NULL) {       /* item not found */
        *item = NULL;           /* set return value in item */
        return HEIGHT_UNCHANGED;
    }

    cmp = compare(*item, (*rootp)->item, ntype);        /* compare item items 
                                                         */

    if (cmp < 0) {              /* delete from left subtree */
        decrease = -avldelete(&(*rootp)->subtree[LEFT], item, compare);
        if (*item == NULL)
            return HEIGHT_UNCHANGED;
    } else if (cmp > 0) {       /* delete from right subtree */
        decrease = avldelete(&(*rootp)->subtree[RIGHT], item, compare);
        if (*item == NULL)
            return HEIGHT_UNCHANGED;
    } else {                    /* cmp == 0 */
        *item = (*rootp)->item; /* set return value in item */

        /*
         *  At this point, we know that "cmp" is zero and "*rootp" points to
         *  the node that we need to delete.  There are three possible cases:
         *
         *     1) The node is a leaf.  Remove it and return.
         *
         *     2) The node is a branch (has only 1 child). Make "*rootp"
         *        (the pointer to this node) point to the child.
         *
         *     3) The node has two children. We swap item with the successor of
         *        "*rootp" (the smallest item in its right subtree) and delete
         *        the successor from the right subtree of "*rootp".  The
         *        identifier "decrease" should be reset if the subtree height
         *        decreased due to the deletion of the successor of "rootp".
         */
        switch (ntype) {        /* what kind of node are we removing? */
        case IS_LEAF:
            avlnode_destroy(rootp);     /* free the leaf, its height */
            return HEIGHT_CHANGED;      /* changes from 1 to 0, return 1 */

        case IS_RBRANCH:       /* only child becomes new root */
        case IS_LBRANCH:
            *rootp = (*rootp)->subtree[dir];
            avlnode_destroy(&old_root); /* free the deleted node */
            return HEIGHT_CHANGED;      /* we just intened the "dir" subtree */

        case IS_TREE:
            decrease = avldelete(&(*rootp)->subtree[RIGHT], &((*rootp)->item), avlmin);
            break;
        case IS_NULL:
            assert(ntype != IS_NULL);
        }
    }

    (*rootp)->bal -= decrease;  /* update balance factor */

    /*
     * rebalance if necessary -- the height of current tree changes if one
     * of two things happens: (1) a rotation was performed which changed
     * the height of the subtree (2) the subtree height decreased and now
     * matches the height of its other subtree (so the current tree now
     * has a zero balance when it previously did not)
     */
    if (decrease && (*rootp)->bal)
        return avlbalance(rootp);       /* rebalance and see if height
                                         * changed */
    else if (decrease && !(*rootp)->bal)
        return HEIGHT_CHANGED;  /* balanced because subtree decreased */
    else
        return HEIGHT_UNCHANGED;
}

/* Traverses the given tree performing "action" upon each item item
 * encountered.
 */
static void avlwalk(avlnode* tree, avl_action action, SIBLING_ORDER order, int level)
{
    DIRECTION dir1 = (order == MIN_TO_MAX) ? LEFT : RIGHT;
    DIRECTION dir2 = OPPOSITE(dir1);

    if (tree == NULL || action == NULL)
        return;

    if (tree->subtree[dir1] != NULL)
        avlwalk(tree->subtree[dir1], action, order, level + 1);

    action(tree->item);

    if (tree->subtree[dir2] != NULL)
        avlwalk(tree->subtree[dir2], action, order, level + 1);
}

/* Frees up space for all nodes in a given tree performing "action" upon each
 * item item encountered (only perform "action" if it is a non-null function).
 */
static void avlfree(avlnode** rootp, avl_action action, SIBLING_ORDER order, int level)
{
    DIRECTION dir1 = (order == MIN_TO_MAX) ? LEFT : RIGHT;
    DIRECTION dir2 = OPPOSITE(dir1);

    if (*rootp != NULL) {
        if ((*rootp)->subtree[dir1] != NULL)
            avlfree(&((*rootp)->subtree[dir1]), action, order, level + 1);
        if (action != NULL)
            action((*rootp)->item);
        if ((*rootp)->subtree[dir2] != NULL)
            avlfree(&((*rootp)->subtree[dir2]), action, order, level + 1);
        free(*rootp);
    }
}

/* Constructor.
 */
avltree* avltree_create(avl_compare compare)
{
    avltree* tree;

    tree = malloc(sizeof(avltree));
    tree->root = NULL;
    tree->compare = compare;
    tree->count = 0;

    return tree;
}

void avltree_destroy(avltree* tree, int deep)
{
    avlfree(&(tree->root), (deep) ? free : NULL, MIN_TO_MAX, 1);
    free(tree);
}

/* Traverses the tree and performs the specified action on each item item in
 * the tree.
 */
void avltree_walk(avltree* tree, avl_action action, SIBLING_ORDER order)
{
    avlwalk(tree->root, action, order, 1);
}

/* Counts the number of nodes in the tree.
 */
int avltree_getnitems(avltree* tree)
{
    return tree->count;
}

/* Inserts the item into the tree.
 */
void* avltree_insertitem(avltree* tree, void* item)
{
    avlinsert(&tree->root, &item, tree->compare);
    if (item == NULL)
        ++(tree->count);

    return item;
}

/* Deletes an item from the tree.
 */
void* avltree_deleteitem(avltree* tree, void* item)
{
    avldelete(&tree->root, &item, tree->compare);
    if (item != NULL)
        --(tree->count);

    return item;
}

/* Finds an item in the tree. 
 * Returns its address (NULL if not found).
 */
void* avltree_finditem(avltree* tree, void* item)
{
    return avlfind(tree->root, item, tree->compare);
}

/* Finds the next item in the tree. 
 * Returns its address (NULL if not found).
 */
void* avltree_findnextitem(avltree* tree, void* item)
{
    avlnode* root = tree->root;
    avlnode* parent = NULL;
    int cmp;

    while (root != NULL && (cmp = tree->compare(item, root->item)) != 0) {
        if (cmp < 0)
            parent = root;
        root = root->subtree[(cmp < 0) ? LEFT : RIGHT];
    }

    if (root == NULL)
        return NULL;
    if (root->subtree[RIGHT] != NULL)
        return avlfind(root->subtree[RIGHT], NULL, avlmin);
    if (parent != NULL)
        return parent->item;
    return NULL;
}

/* Finds the previous item in the tree. 
 * Returns its address (NULL if not found).
 */
void* avltree_findprevitem(avltree* tree, void* item)
{
    avlnode* root = tree->root;
    avlnode* parent = NULL;
    int cmp;

    while (root != NULL && (cmp = tree->compare(item, root->item)) != 0) {
        if (cmp > 0)
            parent = root;
        root = root->subtree[(cmp < 0) ? LEFT : RIGHT];
    }

    if (root == NULL)
        return NULL;
    if (root->subtree[LEFT] != NULL)
        return avlfind(root->subtree[LEFT], NULL, avlmax);
    if (parent != NULL)
        return parent->item;
    return NULL;
}

/* Deletes the minimal item from the tree.
 */
void* avltree_deleteminitem(avltree* tree)
{
    void* item;

    avldelete(&tree->root, &item, avlmin);
    if (item != NULL)
        --(tree->count);

    return item;
}

/* Finds the minimal item in the tree.
 * Returns its address (NULL if not found).
 */
void* avltree_findminitem(avltree* tree)
{
    return avlfind(tree->root, NULL, avlmin);
}

/* Deletes the maximal item from the tree.
 */
void* avltree_deletemaxitem(avltree* tree)
{
    void* item;

    avldelete(&tree->root, &item, avlmax);
    if (item != NULL)
        --(tree->count);

    return item;
}

/* Finds the maximal item in the tree.
 * Return its address (NULL if not found).
 */
void* avltree_findmaxitem(avltree* tree)
{
    return avlfind(tree->root, NULL, avlmax);
}

#if defined (TEST_AVL)

#define NVALS 50

static void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "error: avltree: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    exit(1);
}

static int* generate_shuffled_numbers(int n)
{
    int* v = malloc(n * sizeof(int));
    int i;

    for (i = 0; i < n; ++i)
        v[i] = i - n / 2;

    for (i = 0; i < n; ++i) {
        int nn = (int) ((double) n * (double) rand() / ((double) RAND_MAX + 1.0));
        int tmp = v[i];

        v[i] = v[nn];
        v[nn] = tmp;
    }

    return v;
}

static void print_item(int n, int v[])
{
    int i;

    printf("  item =");
    for (i = 0; i < n; ++i)
        printf(" %d", v[i]);
    printf("\n");
}

static int compare_int(void* i1, void* i2)
{
    return (*(int*) i1 - *(int*) i2);
}

static void tree_printitem(void* item)
{
    printf(" %d", *(int*) item);
}

static void tree_printall(avltree* tree, SIBLING_ORDER order)
{
    printf("  data =");
    avltree_walk(tree, tree_printitem, order);
    printf("\n");
}

static void searchitem(avltree* tree, int n, int* v, int i)
{
    int* item;

    if (i >= n)
        return;

    printf(" item %d:\n", i);
    item = avltree_finditem(tree, &v[i]);
    printf("    item %d = %d\n", i, *item);
    item = avltree_findnextitem(tree, &v[i]);
    if (item != NULL)
        printf("    next(item %d) = %d\n", i, *item);
    else
        printf("    next(item %d) = <none>\n", i);
    item = avltree_findprevitem(tree, &v[i]);
    if (item != NULL)
        printf("    previous(item %d) = %d\n", i, *item);
    else
        printf("    previous(item %d) = <none>\n", i);
}

int main(int argc, char* argv[])
{
    avltree* tree = NULL;
    int n = NVALS;
    int* v = NULL;
    int* p;
    int i;

    if (argc > 1)
        n = atoi(argv[1]);
    if (n <= 0)
        quit("test: nothing to do (n = 0)\n");

    printf("generating %d shuffled integers:\n", n);
    v = generate_shuffled_numbers(n);
    print_item(n, v);

    printf("generating a tree:\n");
    tree = avltree_create(compare_int);

    printf("inserting the numbers into the tree:\n");
    for (i = 0; i < n; i++)
        avltree_insertitem(tree, &v[i]);

    printf("printing the tree MIN_TO_MAX:\n");
    tree_printall(tree, MIN_TO_MAX);
    printf("printing the tree MAX_TO_MIN:\n");
    tree_printall(tree, MAX_TO_MIN);

    printf("searching for maximal and minimal items:\n");
    p = avltree_findminitem(tree);
    printf("  min = %d\n", *(int*) p);
    p = avltree_findmaxitem(tree);
    printf("  max = %d\n", *(int*) p);

    printf("searching for a specified, next and previous items:\n");
    searchitem(tree, n, v, 2);
    searchitem(tree, n, v, 3);
    searchitem(tree, n, v, 4);

    printf("deleting every third item starting from the first:\n");
    for (i = 0; i < n; i += 3)
        avltree_deleteitem(tree, &v[i]);
    printf("printing the tree:\n");
    tree_printall(tree, MIN_TO_MAX);

    printf("deleting every third item starting from the second:\n");
    for (i = 1; i < n; i += 3)
        avltree_deleteitem(tree, &v[i]);
    printf("printing the tree:\n");
    tree_printall(tree, MIN_TO_MAX);

    printf("deleting minimal and maximal items:\n");
    avltree_deleteminitem(tree);
    avltree_deletemaxitem(tree);
    printf("printing the tree:\n");
    tree_printall(tree, MIN_TO_MAX);

    printf("destroying the tree:\n");
    avltree_destroy(tree, 0);

    free(v);

    return 0;
}

#endif                          /* TEST_AVL */
