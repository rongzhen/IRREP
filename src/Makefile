include ../make.inc

MAKE = make

#-------------------------------------------------------------------------------
# Suffix rules
#-------------------------------------------------------------------------------
.SUFFIXES: .o .f90
.f90.o:
	$(F90) $(F90_OPTS) -c $<

#-------------------------------------------------------------------------------
# Source files
#-------------------------------------------------------------------------------
SRC_modules = modmain.f90 modxcifc.f90

SRC_main = main.f90

SRC_exirrep = exirrep.f90 \
 irchar.f90 irchartab.f90 irclasses.f90 \
 irlkgrpg0.f90 irlkgrp.f90 irkpt.f90 \
 irpgrpname.f90 irrotaxis.f90  \
 irsymops.f90 irtimerev.f90  \
 irwrite.f90 irwritesym.f90  \
 irfindmmi.f90 irvint.f90 irminv.f90 irsu2.f90 irwritedate.f90  \
 irmm.f90   

SRC_routines = \
 autoradmt.f90 findprim.f90 findsym.f90 findsymcrys.f90 genppts.f90 rfinp.f90 \
 rvfcross.f90 seceqn.f90 symrf.f90 symrfir.f90 wavefmt.f90

SRC_routines_flux = \
 elfplot.f90 rhonorm.f90 energy.f90 spinchar.f90 \
 zpotclmt.f90 writegeom.f90 nfftifc.f90 zfftifc.f90 allatoms.f90 gridsize.f90 \
 poteff.f90 genrmesh.f90 readfermi.f90 potcoul.f90 gensfacgp.f90 checkmt.f90 \
 zfinp.f90 match.f90 force.f90 forcek.f90 writeefg.f90 packeff.f90 \
 bandchar.f90 findsymlat.f90 \
 genlofr.f90 atom.f90 writefermi.f90 writekpts.f90 fsmfield.f90 mossbauer.f90 \
 occupy.f90 writelinen.f90 writeinfo.f90 readinput.f90 charge.f90 moment.f90 \
 writesym.f90 genidxlo.f90 gencore.f90 addrhocr.f90 gengvec.f90 \
 genshtmat.f90 plot1d.f90 plot2d.f90 plot3d.f90 updatpos.f90 writeiad.f90 \
 symvect.f90 vecplot.f90 genylmg.f90 linengy.f90 init0.f90 init1.f90 \
 gengpvec.f90 genpmat.f90 ggamt.f90 ggair.f90 genveffig.f90 gencfun.f90 \
 genapwfr.f90 seceqnfv.f90 seceqnsv.f90 seceqnss.f90 getngkmax.f90 \
 dos.f90 rhoinit.f90 potplot.f90 writestate.f90 \
 potxc.f90 zpotcoul.f90 gndstate.f90 rhovalk.f90 readstate.f90 bandstr.f90 \
 writeeval.f90 rhoplot.f90 gwf2val.f90 gwf2cr.f90 \
 rfarray.f90 fermisurf.f90 symrvfir.f90 \
 symrvf.f90 symrfmt.f90 wfplot.f90 hmlrad.f90 olprad.f90 olpistl.f90 olpaa.f90 \
 olpalo.f90 olplolo.f90 hmlistl.f90 hmlaa.f90 hmlalo.f90 hmllolo.f90 init2.f90 \
 writelsj.f90 wavefcr.f90 gdft.f90 gradrf.f90 dbxcplot.f90 effmass.f90 \
 projsbf.f90 seceqnhf.f90 writelat.f90 energynn.f90 findigp0.f90 \
 writeengy.f90 writechg.f90 writeforce.f90 genjlgpr.f90 findkpt.f90 \
 putevalfv.f90 getevalfv.f90 putevalsv.f90 getevalsv.f90 putoccsv.f90 \
 getoccsv.f90 putevecfv.f90 getevecfv.f90 putevecsv.f90 getevecsv.f90 \
 delevec.f90 genvmatk.f90 genvmat.f90 genwfsv.f90 \
 findsymsite.f90 geomplot.f90 hartfock.f90 genkinmat.f90 energykncr.f90

SRC_phonon = \
 phonon.f90 dyntask.f90 phcell.f90 phdisp.f90 readdyn.f90 dynqtor.f90 \
 dynrtoq.f90 dynsymapp.f90 dynsym.f90 dyndiag.f90 phdos.f90 sumrule.f90 \
 writephn.f90

SRC_optics = \
 writepmat.f90 linopt.f90 linoptk.f90 moke.f90

SRC_advanced = \
 vnlrho.f90 vnlrhomt.f90 genwiq2.f90 \
 exxengy.f90 exxengyk.f90 oepmain.f90 oepresk.f90 oepvnl.f90 \
 oepvnlk.f90 oepmagmt.f90 oepmag.f90 zrhoqint.f90

SRC_rdmft =
# rdmft.f90 rdmoccinit.f90 rdmdkdc.f90 rdmwriteengy.f90 rdmwritederiv.f90 \
# rdmvnlwf.f90 rdmdedc.f90 rdmdexcdc.f90 rdmvarywf.f90 rdmengy.f90 \
# rdmengyxc.f90 rdmminwf.f90 rdmdedn.f90 rdmdexcdn.f90 rdmvnl.f90 rdmvnlk.f90 \
# rdmvaryn.f90 rdmminn.f90

SRC_lib = \
 euler.f90 wigner3j.f90 gaunt.f90 gauntyry.f90 r3mm.f90 r3mtm.f90 r3mmt.f90 \
 r3mv.f90 r3mtv.f90 r3cross.f90 r3dist.f90 r3taxi.f90 r3dot.f90 r3minv.f90 \
 r3mdet.f90 r3frac.f90 i3mdet.f90 factnm.f90 factr.f90 hermite.f90 \
 brzint.f90 sphcrd.f90 sphcover.f90 erf.f90 clebgor.f90 \
 sbessel.f90 sbesseldm.f90 genylm.f90 genrlm.f90 radmesh.f90 zmatinp.f90 \
 lopzflm.f90 sortidx.f90 gcd.f90 zfmtinp.f90 rfmtinp.f90 findband.f90 \
 gradzfmt.f90 gradrfmt.f90 ztorflm.f90 rtozflm.f90 zflmconj.f90 rotzflm.f90 \
 polynom.f90 sdelta.f90 stheta.f90 sdelta_mp.f90 stheta_mp.f90 sdelta_fd.f90 \
 stheta_fd.f90 rdiracint.f90 rdiracdme.f90 rdirac.f90 rschrodint.f90 \
 rschroddme.f90 rschrodapp.f90 reciplat.f90 \
 connect.f90 \
 flushifc.f90 spline.f90 writewiq2.f90 rfinterp.f90 rfmtctof.f90 fderiv.f90 \
 mixer.f90 fsmooth.f90 rotaxang.f90 i3minv.f90

SRC_lib_flux =

SRC_xc = \
  xc_pzca.f90 xc_pwca.f90 xc_pbe.f90 xc_am05.f90 xc_xalpha.f90 xc_wc06.f90 \
  x_wc06.f90 x_pbe.f90 c_pbe.f90 c_pbe_gcor.f90

SRC = $(SRC_modules) $(SRC_main) $(SRC_routines) $(SRC_routines_flux) \
 $(SRC_lib) $(SRC_lib_flux) $(SRC_advanced) $(SRC_phonon) $(SRC_optics) \
 $(SRC_xc) $(SRC_rdmft) $(SRC_exirrep) 

OBJ = $(SRC:.f90=.o)
EXE = exciting

exciting:	$(OBJ)
	$(F90) $(F90_OPTS) -o $(EXE) $(OBJ) $(LIB_LPK) $(LIB_FFT) $(LIB_SYS)

blas:
	cd BLAS; $(MAKE); cp blas.a ..

lapack:
	cd LAPACK; $(MAKE); cp lapack.a ..

fft:
	cd fftlib; $(MAKE); cp fftlib.a ..

xc:
	cd xclib; $(MAKE); cp xclib.a ..

all:	blas lapack fft exciting

clean:
	rm -f *.o *.mod *~ fort.* ifc* *.gcno gmon.out *.aux *.dvi *.log *.pdf \
         *.tex *.toc $(EXE)

cleanall:
	cd BLAS; $(MAKE) clean
	cd LAPACK; $(MAKE) clean
	cd fftlib; $(MAKE) clean
	$(MAKE) clean

doc:
	./protex -s $(SRC_main) $(SRC_modules) $(SRC_routines) \
         $(SRC_routines_flux) $(SRC_lib) $(SRC_lib_flux) $(SRC_xc) \
         $(SRC_advanced) $(SRC_phonon) > exciting.tex
	 latex exciting;bibtex exciting;latex exciting;latex exciting
	 dvipdf exciting.dvi
	./protex -s exirrepdoc.f90  $(SRC_exirrep) > exirrep.tex
	 latex exirrep;bibtex exirrep;latex exirrep;latex exirrep   
	 dvipdf exirrep.dvi  

doctemp:
	./protex -s $(SRC_main) $(SRC_routines) $(SRC_advanced)  > exciting.tex
	latex exciting

backup:
	tar -czf ex.tgz $(SRC) fftlib eos spacegroup species junk BLAS LAPACK \
         Makefile notes.txt docs docedit examples protex COPYING README \
         exciting.bib depend ../release

VERSION = $$(awk -F"/" '/data version/ {print $$2}' modmain.f90 | sed 's/ //g;s/,/./g')$
RELEASE = ../release/exciting

release:
	rm -rf $(RELEASE)
	mkdir $(RELEASE)
	cd BLAS; $(MAKE) clean
	cd LAPACK; $(MAKE) clean
	cd fftlib; $(MAKE) clean
	cd eos; $(MAKE) clean
	cd spacegroup; $(MAKE) clean
	cd species; $(MAKE) clean
	mkdir $(RELEASE)/src
	cp -r BLAS $(RELEASE)/src
	cp -r LAPACK $(RELEASE)/src
	cp -r fftlib $(RELEASE)/src
	cp -r eos $(RELEASE)/src
	cp -r spacegroup $(RELEASE)/src
	cp -r species $(RELEASE)/src
	rm -f $(RELEASE)/src/species/*.in
	cp $(SRC) Makefile protex exciting.bib $(RELEASE)/src
	cd examples; find . -type f -name *.OUT -exec rm -f {} \;
	cd examples; find . -type f -name *~ -exec rm -f {} \;
	cp -r examples $(RELEASE)
	mkdir $(RELEASE)/species
	cp species/*.in $(RELEASE)/species
	cp COPYING $(RELEASE)
	cp README $(RELEASE)
	cp ../Makefile $(RELEASE)
	cp ../setup $(RELEASE)
	mkdir $(RELEASE)/utilities
	cp -r NEdit $(RELEASE)/utilities
	$(MAKE) doc
	cp exciting.pdf docs/exciting
	cd spacegroup;$(MAKE) doc;cp spacegroup.pdf ../docs/spacegroup;$(MAKE) clean
	cp -r docs $(RELEASE)
	tar -C ../release -czf ../release/exciting-$(VERSION).tgz exciting
	cp docs/exciting/exciting.pdf ../release
	cp docs/spacegroup/spacegroup.pdf ../release

lines:
	cat $(SRC) | wc -l

