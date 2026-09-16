#ifndef _GEOTIFF_FILEMAP_H_
#define _GEOTIFF_FILEMAP_H_

  enum
  {
    FILEMAP_TIFF,     /* geotiff file mapped via mmap ()            */
    FILEMAP_PIXELS,   /* 2D pixel data (decoded from geotiff file)  */
  } ;

  /* api */
  char * filemap_address   ( int type );
  size_t filemap_size      ( int type );
  int    filemap_tiff      ( const char * tiff );
  int    filemap_pixels    ( size_t npixels );
  void   filemap_close_all ( const char * error_if_any );
         /* in case error_if_any != NULL, it will exit (-1) */
  
#endif
