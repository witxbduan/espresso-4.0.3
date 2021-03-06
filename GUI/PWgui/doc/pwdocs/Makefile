PW_DOC_DIR = ../../../../Doc
LD1_DOC_DIR = ../../../../atomic_doc

INPUT_FILES_HTML = \
        INPUT_PW.html \
        INPUT_PH.html \
        INPUT_PP.html \
        INPUT_PROJWFC.html \
        INPUT_LD1.html \
        INPUT_D3.html

INPUT_FILES_TXT = $(INPUT_FILES_HTML:.html=.txt)

IMAGES     = democritos.png pwscf.png
GIF_IMAGES = $(IMAGES:.png=.gif)
LINK_FILES = $(INPUT_FILES_TXT) 


%.pdf : %.tex
	if test ! -f $@ ; then pdflatex $<; fi
	pdflatex $<

%.gif : %.png
	convert $< $@


all: input_html compile_user_guide


links:
	( cd $(PW_DOC_DIR); $(MAKE) )
	for file in $(LINK_FILES); do ln -sf $(PW_DOC_DIR)/$$file . ; done
#	ln -sf $(LD1_DOC_DIR)/INPUT_LD1 .


clean_links:
	for file in $(LINK_FILES); do \
		if test -f $$file ; then rm -f $$file ; fi ; \
	done


input_html: links 
	./htmlise.sh PW      pw.x      INPUT_PW.txt      > INPUT_PW.html
	./htmlise.sh PH      ph.x      INPUT_PH.txt      > INPUT_PH.html
	./htmlise.sh PP      pw.x      INPUT_PP.txt      > INPUT_PP.html
	./htmlise.sh ProjWFC projwfc.x INPUT_PROJWFC.txt > INPUT_PROJWFC.html
	./htmlise.sh D3      d3.x      INPUT_D3.txt      > INPUT_D3.html
	./htmlise.sh LD1     ld1.x     INPUT_LD1.txt     > INPUT_LD1.html

compile_user_guide: $(GIF_IMAGES) user_guide



$(GIF_IMAGES): $(IMAGES)


#user_guide.tex: user_guide.html
#	./html2latex.sh $<

user_guide: user_guide.tex
	if test ! -f user_guide.pdf ; then pdflatex user_guide; fi
	pdflatex user_guide
	rm -rf user_guide/	
	latex2html \
		-t "User's Guide for Quantum-ESPRESSO" \
		-html_version 3.2,math \
		-toc_depth 5 -split 5 -toc_stars -show_section_numbers \
		-local_icons \
		user_guide.tex
	cd user_guide; \
	for file in *.png; do \
		convert $$file $${file%.png}.gif; \
		rm -f $$file; \
	done; \
	for file in *.html; do \
		cp $$file /tmp/$$file; \
		cat /tmp/$$file | sed 's/\.png/\.gif/g' - | sed 's/HREF="http/NAME="http/g' | sed 's/mathend000#//g' | sed 's/&nbsp;/ /g' - > $$file; \
		rm -f /tmp/$$file; \
	done


clean_user_guide:
	-rm -f user_guide.aux user_guide.log user_guide.out user_guide.toc *.gif


clean: clean_links clean_user_guide
	-rm -f INPUT_*.html *~
	-rm -rf user_guide.pdf user_guide/	

distclean: clean
veryclean: clean