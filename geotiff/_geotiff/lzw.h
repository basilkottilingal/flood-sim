#ifndef _GEOTIFF_LZWD_H_
#define _GEOTIFF_LZWD_H_

  /*
  .. APIs.
  .. (a) decode compressed lzw encoding.
  ..     Return value = 0 if successful (in which case return value is the
  ..     length of the decoded data) and return value < 0 in case of error
  .. (b) error handling in case return value of lzw_decode () < 0.
  ..     Use the return val of lzw_decode () as the argument.
  */
  int  lzw_decode (const void * encoding, size_t len, void * buff, size_t outlen);
  void lzw_error  (int e);

#endif 
