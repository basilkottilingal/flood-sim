#!/bin/sh
# webgl-assets.c: the viewer's template and script, and three.js, as C
# arrays (xxd -i), so that a program writes its HTML with no file beside
# it. Run from libwebgl: sh gen-assets.sh [out]. The names follow the
# paths: viewer/viewer.html -> viewer_viewer_html, viewer_viewer_html_len.
#
# three.js (vendor/three.module.min.js, vendor/OrbitControls.js) is
# embedded when it is there (sh vendor/fetch.sh), and the pages work
# offline; when it is not, or with WEBGL_ONLINE=1, its arrays are empty
# and the pages load it from jsDelivr (webgl_offline () tells which).
set -e
out=${1:-webgl-assets.c}
{
  echo "/* made by gen-assets.sh from viewer/ and vendor/: do not edit */"
  for f in viewer/viewer.html viewer/viewer.js; do
    xxd -i "$f"
  done
  for f in vendor/three.module.min.js vendor/OrbitControls.js; do
    name=$(echo "$f" | sed 's/[^A-Za-z0-9]/_/g')
    if [ -z "$WEBGL_ONLINE" ] && [ -f "$f" ]; then
      xxd -i "$f"
    else
      echo "unsigned char $name [] = { 0 };"
      echo "unsigned int ${name}_len = 0;"
    fi
  done
} > "$out.tmp"
mv "$out.tmp" "$out"
