# libwebgl: webgl.c (plain C) and the adapters webgl-tree.c (libb-amr) and
# webgl-hmesh.c (libhmesh), built with cwrap as the other libraries are.
# webgl-assets.c is made from viewer/ and vendor/ by gen-assets.sh.
DEPS := b-amr hmesh utils
CFLAGS += -D_GNU_SOURCE -dimension=2
CPPFLAGS += -I../libb-amr -I../libhmesh
include ../lib.mk

# three.js embedded if vendor/ has it (sh vendor/fetch.sh), else loaded
# from jsDelivr by the pages
webgl-assets.c: gen-assets.sh viewer/viewer.html viewer/viewer.js \
                $(wildcard vendor/*.js)
	sh gen-assets.sh $@
