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
 *     Added avltree_findnextitem() and avltree_findprevitem().
 *
******************************************************************************/

#if !defined(_AVL_H)
#define _AVL_H

typedef enum { MIN_TO_MAX, MAX_TO_MIN } SIBLING_ORDER;

struct avltree;
typedef struct avltree avltree;

typedef int (*avl_compare) (void* item1, void* item2);
typedef void (*avl_action) (void* item);

avltree* avltree_create(avl_compare compare);
void avltree_destroy(avltree* tree, int deep);

void avltree_walk(avltree*, avl_action action, SIBLING_ORDER order);
int avltree_getnitems(avltree* tree);

/* returns NULL on success , `item' otherwise */
void* avltree_insertitem(avltree* tree, void* item);

/* returns `item' on success, NULL otherwise */
void* avltree_deleteitem(avltree* tree, void* item);

/* returns the sought item on success, NULL otherwise */
void* avltree_finditem(avltree* tree, void* image);

/* returns the sought item on success, NULL otherwise */
void* avltree_findnextitem(avltree* tree, void* item);

/* returns the sought item on success, NULL otherwise */
void* avltree_findprevitem(avltree* tree, void* item);

/* returns the sought item on success, NULL otherwise */
void* avltree_findminitem(avltree*);

/* returns the sought item on success, NULL otherwise */
void* avltree_findmaxitem(avltree* tree);

/* returns the item on success, NULL otherwise */
void* avltree_deleteminitem(avltree* tree);

/* returns the deleted item on success, NULL otherwise */
void* avltree_deletemaxitem(avltree* tree);

#endif
