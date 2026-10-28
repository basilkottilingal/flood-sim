#ifndef _GEOJSON_H_
#define _GEOJSON_H_

  #include "ast.h"
  #include "json.tab.h"

  typedef struct GJHashEntry
  {
    const char *key;
    const Ast * geometry;
    struct GJHashEntry *next;
  } GJHashEntry;

  typedef struct
  { 
    GJHashEntry **b;
    size_t nb, n, mask;
  } GJHashTable;

  typedef struct
  {
    GJHashTable hashtable;
    Ast * tree;
  } GeoJSON;

  /* API */
  int   geojson (const char * json, GeoJSON * gj, const char * uniqueTag);
  const Ast * geojson_locate (GeoJSON * gj, const char * byname);

#endif
