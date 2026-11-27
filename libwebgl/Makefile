# libwebgl: webgl.c (plain C) and the adapters webgl-tree.c (libb-amr) and
# webgl-hmesh.c (libhmesh), built with cwrap as the other libraries are.
# webgl-assets.c is made from viewer/ and vendor/ by gen-assets.sh.
DEPS := b-amr hmesh utils
CFLAGS += -D_GNU_SOURCE -dimension=2
CPPFLAGS += -I../libb-amr -I../libhmesh
include ../lib.mk

# the command: bamr2html FILE.bamr [OUT.html] (plain C, no libb-amr)
bamr2html: tools/bamr2html.c webgl.c webgl-bamr.c webgl-assets.c webgl.h \
           webgl-bamr.h
	gcc -O2 -std=gnu11 -D_GNU_SOURCE -I. tools/bamr2html.c webgl.c \
	  webgl-bamr.c webgl-assets.c -lm -o $@

# three.js embedded if vendor/ has it (sh vendor/fetch.sh), else loaded
# from jsDelivr by the pages
webgl-assets.c: gen-assets.sh viewer/viewer.html viewer/viewer.js \
                $(wildcard vendor/*.js)
	sh gen-assets.sh $@
