#ifndef HEADER_hash
#define HEADER_hash

typedef struct _hashelem
{
  char             *name;
  int               index;
  struct _hashelem *next;
  struct _hashelem *nextelem;
} hashelem;

typedef struct _hashtable
{
  hashelem         **table;
  int              size;
  int              base;
  int              count;
  struct _hashelem *first;
  struct _hashelem *last;
} hashtable;

#ifdef __cplusplus
extern "C" {
#endif

hashtable *create_hash_table(int size, int base);
void      free_hash_table(hashtable *ht);
hashelem  *findhash(const char *name, hashtable *ht);
hashelem  *puthash(const char *name, int index, hashelem **list, hashtable *ht);
void      drophash(const char *name, hashelem **list, hashtable *ht);
void      free_hash_item(hashelem **hp);
hashtable *copy_hash_table(hashtable *ht, hashelem **list, int newsize);

#ifdef __cplusplus
 }
#endif

#endif /* HEADER_lp_hash */
