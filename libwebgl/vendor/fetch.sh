#!/bin/sh
# three.js r160 for libwebgl, from jsDelivr (the npm package three@0.160.0),
# checked against the sha256 sums in README.md; then rebuild (make, or sh
# tests/run.sh) and the pages work offline. Run from libwebgl:
#   sh vendor/fetch.sh
set -e
url=https://cdn.jsdelivr.net/npm/three@0.160.0
cd "$(dirname "$0")"
curl -fsSL -o three.module.min.js "$url/build/three.module.min.js"
curl -fsSL -o OrbitControls.js "$url/examples/jsm/controls/OrbitControls.js"
sha256sum -c <<SUMS
3e690ac7d180b0aadf0891bea39eec643e29e2d3e75c99b18689518665f69ba6  three.module.min.js
5a44a9e86a2a0fb11933eed69bc2cd33c76a496854c1aed6ed776efa87d7b064  OrbitControls.js
SUMS
