#ifndef _GEOJSON_H_
#define _GEOJSON_H_

  #include "ast.h"

  typedef struct
  {
  }

  typedef struct
  {
    Ast ** geometry;
    Ast * json;
  } GeoJSON;

  /* API */
  int   geojson (const char * json, GeoJSON * g);
  Ast * geojson_locate (const char * byname);
  

#endif
