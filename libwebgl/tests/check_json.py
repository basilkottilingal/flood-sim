# The JSON in a written HTML file, parsed by Python's own parser: the
# manifest (its layers, fields and blob ids against the blobs in the file)
# and the import map. Usage: python3 -I tests/check_json.py file.html
import json, re, sys
html = open (sys.argv[1], encoding = 'utf-8').read ()
m = re.search (r'<script id="manifest" type="application/json">(.*?)</script>',
               html, re.S)
man = json.loads (m.group (1))
blobs = set (re.findall (r'<script id="(b\d+)" type="application/octet-stream">',
                         html))
ids = []
for L in man['layers']:
    for k in ('box', 'xy', 'e'):
        if k in L:
            ids.append (L[k])
    for f in L['fields']:
        ids.append (f['data'])
        assert f['where'] in ('cell', 'block', 'point', 'element'), f
assert sorted (ids) == sorted (blobs), (sorted (ids), sorted (blobs))
im = re.search (r'<script type="importmap">(.*?)</script>', html, re.S)
imports = json.loads (im.group (1))['imports']
assert set (imports) == {'three', 'three/addons/controls/OrbitControls.js'}
assert all (v.startswith ('data:text/javascript;base64,') or
            v.startswith ('https://cdn.jsdelivr.net/npm/three@0.160.0/')
            for v in imports.values ()), imports
print ('json ok: %d layers, %d blobs, title %r'
       % (len (man['layers']), len (blobs), man['title']))
