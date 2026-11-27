#ifndef _WEBGL_ASSETS_H_
#define _WEBGL_ASSETS_H_

  /*
  .. The files webgl_scene_write () puts in its HTML, compiled in by
  .. gen-assets.sh (xxd -i) as webgl-assets.c: the page's template, the
  .. viewer, three.js and its OrbitControls (vendor/, MIT, revision 160).
  */

  extern unsigned char viewer_viewer_html [];
  extern unsigned int  viewer_viewer_html_len;
  extern unsigned char viewer_viewer_js [];
  extern unsigned int  viewer_viewer_js_len;
  extern unsigned char vendor_three_module_min_js [];
  extern unsigned int  vendor_three_module_min_js_len;
  extern unsigned char vendor_OrbitControls_js [];
  extern unsigned int  vendor_OrbitControls_js_len;

#endif
