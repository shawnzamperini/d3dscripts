#  VERSION : 21.10.2000 00:50
# Modified by D. P. Stotler ca. 2006

SHELL=/bin/sh
SYSTEM:=$(shell uname)

VPATH = src:src/shared
INCLUDE =
MSCL=../../mscl/$(MACHINE)

ifeq ($(SYSTEM),SunOS)
  MACHINE=SUN
  FC=f90
  CPP=/usr/ccs/lib/cpp
  FFLAGS=-g
  LDLIBS=-L$(MSCL)
  LDFLAGS=-lmscl -lf77compat
endif
ifeq ($(SYSTEM),AIX)
  MACHINE=IBM
# Kukushkin has flag settings, but rather complicated and may be out of date
endif
ifeq ($(SYSTEM),OSF1)
  MACHINE=ALPHA
endif
ifeq ($(SYSTEM),HP-UX)
  MACHINE=HP
  FC=f77
  CPP=/lib/cpp
  FFLAGS=-g -C -R8 +FPOZV
endif
ifeq ($(SYSTEM),IRIX)
  MACHINE=SGI
endif
ifeq ($(SYSTEM),Linux)
  MACHINE=LINUX
  CPP=/lib/cpp -C -traditional
  LDLIBS=-L$(MSCL)
  LDFLAGS=-lmscl
#
#  See about setting the compiler via an environment variable
#  COMPILER=PATHSCALE
#  COMPILER=PGROUP
#
# Note: have detected an incompatibility between PathScale and NCAR. 
# Suggest compiling without NCAR when using PathScale Fortran.
#
  ifeq ($(COMPILER),PATHSCALE)
    FC=pathf90
    FFLAGS=-g 
  endif
  ifeq ($(COMPILER),PGROUP)
    FC=pgf90
    FFLAGS=-g
  endif 
endif

ALL = cropequ dg2dg dg2ef ef2dg vr2dg nk2dg dg2vr risepsi prinequ jt2dg pb2dg pt2dg

creq =  cropequ.o rdeqdg.o wreqdg.o rdeqlh.o
dgdg =	dg2dg.o rdeqdg.o wreqdg.o rdeqlh.o
dgminuspsidg =	dg_minuspsi_dg.o rdeqdg.o wreqdg.o rdeqlh.o
dgef =	dg2ef.o rdeqdg.o wrefit.o rdeqlh.o
efdg =	ef2dg.o wreqdg.o rdefit.o
vrdg =	vr2dg.o wreqdg.o rdeqvr.o
nkdg =	nk2dg.o wreqdg.o rdeqnk.o
jtdg =	jt2dg.o rdeqjt.o wreqdg.o rdeqjh.o
dgvr =	dg2vr.o rdeqdg.o wreqvr.o rdeqlh.o
rsps =	risepsi.o rdeqdg.o wreqdg.o rdeqlh.o
preq =	prinequ.o rdeqdg.o rdeqlh.o
pbdg =	pb2dg.o rdeqpb.o wreqdg.o
ptdg =	pt2dg.o
culham = culham.o rdefit.o wreqdg.o

include LISTOBJ

CREQ = $(creq:%.o=$(MACHINE)/%.o)
DGDG = $(dgdg:%.o=$(MACHINE)/%.o)
DGMINUSPSIDG = $(dgminuspsidg:%.o=$(MACHINE)/%.o)
DGEF = $(dgef:%.o=$(MACHINE)/%.o)
DGVR = $(dgvr:%.o=$(MACHINE)/%.o)
EFDG = $(efdg:%.o=$(MACHINE)/%.o)
VRDG = $(vrdg:%.o=$(MACHINE)/%.o)
NKDG = $(nkdg:%.o=$(MACHINE)/%.o)
RSPS = $(rsps:%.o=$(MACHINE)/%.o)
PREQ = $(preq:%.o=$(MACHINE)/%.o)
JTDG = $(jtdg:%.o=$(MACHINE)/%.o)
PBDG = $(pbdg:%.o=$(MACHINE)/%.o)
PTDG = $(ptdg:%.o=$(MACHINE)/%.o)
CULHAM = $(culham:%.o=$(MACHINE)/%.o)

COMPILE=${FC} -c ${FFLAGS}
$(MACHINE)/%.o : %.F
	${CPP} -P -C ${INCLUDE} $< >${MACHINE}/$*.f
	$(COMPILE) -o ${MACHINE}/$*.o ${MACHINE}/$*.f
	if [ -f $*.o ]; then /bin/mv $*.o ${MACHINE}; fi

all: $(ALL:%=$(MACHINE)/%)

${MACHINE}/cropequ: $(CREQ)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/dg2dg: $(DGDG)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/culham: $(CULHAM)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/dg2ef: $(DGEF)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/dg2vr: $(DGVR)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/ef2dg: $(EFDG)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/vr2dg: $(VRDG)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/nk2dg: $(NKDG)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/jt2dg: $(JTDG)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/risepsi: $(RSPS)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/prinequ: $(PREQ)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/pb2dg: $(PBDG)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/pt2dg: $(PTDG)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

${MACHINE}/dg_minuspsi_dg: $(DGMINUSPSIDG)
	${FC} ${FFLAGS} -o $@ $^ ${LDLIBS} $(LDFLAGS) 

clean:
	cd ${MACHINE}; rm -rf *.o $(ALL)

tags:
	rm -f TAGS ; etags src/*.F src/shared/*.F

depend: ${OBJS:.o=.f}
	echo ${OBJS}
	echo ${OBJS:.o=.f}
	makedepend -f  ${MACHINE}/dependencies ${INCLUDE} $^
	mv ${MACHINE}/dependencies ${MACHINE}/dependencies.bak
	sed -e '3,$$s/^/${MACHINE}\//' ${MACHINE}/dependencies.bak > ${MACHINE}/dependencies

# listobj:
#	@cd src ; echo "OBJS =" *.f | sed -e 's/ [^ /]*\// /g' -e 's/\.f/.o/g' > ../LISTOBJ

listobj:
	@rm -f LISTOBJ; touch LISTOBJ; l="OBJS ="; \
	for d in `echo "$(VPATH)" | tr : \ `; do \
		l="$$l `(cd $$d > /dev/null; echo *.F)`"; \
	done; \
	E="-e 's/\.F/\.o/g'" ; \
	echo "$$l" | eval sed "$$E" > LISTOBJ

LISTOBJ: listobj


include ${MACHINE}/dependencies

