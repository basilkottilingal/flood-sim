# The red-black Gauss-Seidel design note: ../red-black-gauss-seidel.pdf
# (not tracked: see ../.gitignore). Needs python3 with numpy, matplotlib,
# reportlab, and cwrap built (../../../libc-pp).
CWRAP ?= ../../../libc-pp/build/cwrap
AMR   := ../../../libb-amr
SRC   := $(AMR)/tree.c $(AMR)/grid-scalar.c $(AMR)/memindex.c $(AMR)/voxel-memindex.c \
         $(AMR)/nano-vdb.c $(AMR)/pool.c $(AMR)/block-ghosts.c $(AMR)/adapt.c \
         $(AMR)/../libutils/bits.c
OUT   := ../red-black-gauss-seidel.pdf

$(OUT): build_pdf.py figures.py lfa.json cost2.json cost3.json
	python3 figures.py
	python3 build_pdf.py $(OUT)

lfa.json: lfa.py
	python3 lfa.py > $@

cost%.json: cost.c
	$(CWRAP) -dimension=$* -D_GNU_SOURCE -std=gnu11 -O2 -I$(AMR) -I$(AMR)/.. cost.c $(SRC) -lm -o cost$*
	./cost$* > $@

clean:
	$(RM) fig*.png lfa.json cost2.json cost3.json cost2 cost3 $(OUT)

.PHONY: clean
