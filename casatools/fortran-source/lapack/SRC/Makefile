include ../make.inc

clean:
	rm -f *.o liblapack*

FC=gfortran
OPT=-O2
FFLAGS = -fPIC $(OPT)

OBJS=	cbdsqr.o cgbbrd.o cgbcon.o cgbequb.o cgbequ.o cgbrfs.o \
	cgbsv.o cgbsvx.o cgbtf2.o cgbtrf.o cgbtrs.o cgebak.o \
	cgebal.o cgebd2.o cgebrd.o cgecon.o cgeequb.o cgeequ.o cgees.o \
	cgeesx.o cgeev.o cgeevx.o cgegs.o cgegv.o cgehd2.o cgehrd.o cgelq2.o \
	cgelqf.o cgelsd.o cgels.o cgelss.o cgelsx.o cgelsy.o cgeql2.o \
	cgeqlf.o cgeqp3.o cgeqpf.o cgeqr2.o cgeqrf.o cgerfs.o \
	cgerq2.o cgerqf.o cgesc2.o cgesdd.o cgesvd.o cgesv.o cgesvx.o \
	cgetc2.o cgetf2.o cgetrf.o cgetri.o cgetrs.o cggbak.o \
	cggbal.o cgges.o cggesx.o cggev.o cggevx.o cggglm.o cgghrd.o \
	cgglse.o cggqrf.o cggrqf.o cggsvd.o cggsvp.o cgtcon.o cgtrfs.o \
	cgtsv.o cgtsvx.o cgttrf.o cgttrs.o cgtts2.o chbevd.o chbev.o \
	chbevx.o chbgst.o chbgvd.o chbgv.o chbgvx.o chbtrd.o checon.o \
	cheequb.o cheevd.o cheev.o cheevr.o cheevx.o chegs2.o chegst.o \
	chegvd.o chegv.o chegvx.o cherfs.o chesv.o chesvx.o \
	chetd2.o chetf2.o chetrd.o chetrf.o chetri.o chetrs.o \
	chfrk.o chgeqz.o chla_transtype.o chpcon.o chpevd.o chpev.o \
	chpevx.o chpgst.o chpgvd.o chpgv.o chpgvx.o chprfs.o chpsv.o \
	chpsvx.o chptrd.o chptrf.o chptri.o chptrs.o chsein.o chseqr.o \
	clabrd.o clacgv.o clacn2.o clacon.o clacp2.o clacpy.o clacrm.o \
	clacrt.o cladiv.o claed0.o claed7.o claed8.o claein.o claesy.o \
	claev2.o clag2z.o \
	clags2.o clagtm.o \
	clahef.o \
	clahqr.o clahr2.o clahrd.o claic1.o \
	clals0.o clalsa.o clalsd.o clangb.o clange.o clangt.o clanhb.o \
	clanhe.o clanhf.o clanhp.o clanhs.o clanht.o clansb.o clansp.o \
	clansy.o clantb.o clantp.o clantr.o clapll.o clapmt.o \
	claqgb.o claqge.o claqhb.o claqhe.o claqhp.o \
	claqp2.o claqps.o claqr0.o claqr1.o claqr2.o claqr3.o claqr4.o \
	claqr5.o claqsb.o claqsp.o claqsy.o clar1v.o clar2v.o clarcm.o \
	clarfb.o clarf.o clarfg.o clarfp.o clarft.o clarfx.o clargv.o \
	clarnv.o clarrv.o clartg.o clartv.o \
	clarzb.o clarz.o clarzt.o clascl.o claset.o clasr.o \
	classq.o claswp.o clasyf.o \
	clatbs.o \
	clatdf.o clatps.o clatrd.o clatrs.o clatrz.o clatzm.o clauu2.o \
	clauum.o cpbcon.o cpbequ.o cpbrfs.o cpbstf.o cpbsv.o \
	cpbsvx.o cpbtf2.o cpbtrf.o cpbtrs.o cpftrf.o cpftri.o cpftrs.o \
	cpocon.o cpoequb.o cpoequ.o cporfs.o cposv.o cposvx.o \
	cpotf2.o cpotrf.o cpotri.o cpotrs.o cppcon.o cppequ.o \
	cpprfs.o cppsv.o cppsvx.o cpptrf.o cpptri.o cpptrs.o cpstf2.o \
	cpstrf.o cptcon.o cpteqr.o cptrfs.o cptsv.o cptsvx.o cpttrf.o \
	cpttrs.o cptts2.o crot.o cspcon.o cspmv.o cspr.o csprfs.o cspsv.o \
	cspsvx.o csptrf.o csptri.o csptrs.o csrscl.o cstedc.o cstegr.o \
	cstein.o cstemr.o csteqr.o csycon.o csyequb.o csymv.o csyr.o \
	csyrfs.o csysv.o csysvx.o csytf2.o csytrf.o \
	csytri.o csytrs.o ctbcon.o ctbrfs.o ctbtrs.o ctfsm.o ctftri.o \
	ctfttp.o ctfttr.o ctgevc.o ctgex2.o ctgexc.o ctgsen.o ctgsja.o \
	ctgsna.o ctgsy2.o ctgsyl.o ctpcon.o ctprfs.o ctptri.o ctptrs.o \
	ctpttf.o ctpttr.o ctrcon.o ctrevc.o ctrexc.o ctrrfs.o ctrsen.o \
	ctrsna.o ctrsyl.o ctrti2.o ctrtri.o ctrtrs.o ctrttf.o ctrttp.o \
	ctzrqf.o ctzrzf.o cung2l.o cung2r.o cungbr.o cunghr.o cungl2.o \
	cunglq.o cungql.o cungqr.o cungr2.o cungrq.o cungtr.o cunm2l.o \
	cunm2r.o cunmbr.o cunmhr.o cunml2.o cunmlq.o cunmql.o cunmqr.o \
	cunmr2.o cunmr3.o cunmrq.o cunmrz.o cunmtr.o cupgtr.o cupmtr.o \
	dbdsdc.o dbdsqr.o ddisna.o dgbbrd.o dgbcon.o dgbequb.o dgbequ.o \
	dgbrfs.o dgbsv.o dgbsvx.o dgbtf2.o dgbtrf.o \
	dgbtrs.o dgebak.o dgebal.o dgebd2.o dgebrd.o dgecon.o dgeequb.o \
	dgeequ.o dgees.o dgeesx.o dgeev.o dgeevx.o dgegs.o dgegv.o \
	dgehd2.o dgehrd.o dgejsv.o dgelq2.o dgelqf.o dgelsd.o dgels.o \
	dgelss.o dgelsx.o dgelsy.o dgeql2.o dgeqlf.o dgeqp3.o dgeqpf.o \
	dgeqr2.o dgeqrf.o dgerfs.o dgerq2.o dgerqf.o dgesc2.o \
	dgesdd.o dgesvd.o dgesv.o dgesvj.o dgesvx.o dgetc2.o \
	dgetf2.o dgetrf.o dgetri.o dgetrs.o dggbak.o dggbal.o dgges.o \
	dggesx.o dggev.o dggevx.o dggglm.o dgghrd.o dgglse.o dggqrf.o \
	dggrqf.o dggsvd.o dggsvp.o dgsvj0.o dgsvj1.o dgtcon.o dgtrfs.o \
	dgtsv.o dgtsvx.o dgttrf.o dgttrs.o dgtts2.o dhgeqz.o dhsein.o \
	dhseqr.o disnan.o dlabad.o dlabrd.o dlacn2.o dlacon.o dlacpy.o \
	dladiv.o dlae2.o dlaebz.o dlaed0.o dlaed1.o dlaed2.o dlaed3.o \
	dlaed4.o dlaed5.o dlaed6.o dlaed7.o dlaed8.o dlaed9.o dlaeda.o \
	dlaein.o dlaev2.o dlaexc.o dlag2.o dlag2s.o \
	dlags2.o dlagtf.o dlagtm.o \
	dlagts.o dlagv2.o dlahqr.o dlahr2.o dlahrd.o dlaic1.o \
	dlaisnan.o dlaln2.o dlals0.o dlalsa.o dlalsd.o \
	dlamrg.o dlaneg.o dlangb.o dlange.o dlangt.o dlanhs.o dlansb.o \
	dlansf.o dlansp.o dlanst.o dlansy.o dlantb.o dlantp.o dlantr.o \
	dlanv2.o dlapll.o dlapmt.o \
	dlapy2.o dlapy3.o dlaqgb.o dlaqge.o dlaqp2.o \
	dlaqps.o dlaqr0.o dlaqr1.o dlaqr2.o dlaqr3.o dlaqr4.o dlaqr5.o \
	dlaqsb.o dlaqsp.o dlaqsy.o dlaqtr.o dlar1v.o dlar2v.o dlarfb.o \
	dlarf.o dlarfg.o dlarfp.o dlarft.o dlarfx.o dlargv.o dlarnv.o \
	dlarra.o dlarrb.o dlarrc.o dlarrd.o dlarre.o \
	dlarrf.o dlarrj.o dlarrk.o dlarrr.o dlarrv.o \
	dlartg.o dlartv.o dlaruv.o dlarzb.o dlarz.o dlarzt.o dlas2.o \
	dlascl.o dlasd0.o dlasd1.o dlasd2.o dlasd3.o \
	dlasd4.o dlasd5.o dlasd6.o dlasd7.o dlasd8.o dlasda.o dlasdq.o \
	dlasdt.o dlaset.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o dlasq5.o \
	dlasq6.o dlasr.o dlasrt.o dlassq.o dlasv2.o dlaswp.o dlasy2.o \
	dlasyf.o \
	dlat2s.o dlatbs.o dlatdf.o dlatps.o dlatrd.o \
	dlatrs.o dlatrz.o dlatzm.o dlauu2.o dlauum.o \
	dopgtr.o dopmtr.o dorg2l.o dorg2r.o dorgbr.o dorghr.o dorgl2.o \
	dorglq.o dorgql.o dorgqr.o dorgr2.o dorgrq.o dorgtr.o dorm2l.o \
	dorm2r.o dormbr.o dormhr.o dorml2.o dormlq.o dormql.o dormqr.o \
	dormr2.o dormr3.o dormrq.o dormrz.o dormtr.o dpbcon.o dpbequ.o \
	dpbrfs.o dpbstf.o dpbsv.o dpbsvx.o dpbtf2.o dpbtrf.o dpbtrs.o \
	dpftrf.o dpftri.o dpftrs.o dpocon.o dpoequb.o dpoequ.o \
	dporfs.o dposv.o dposvx.o dpotf2.o dpotrf.o \
	dpotri.o dpotrs.o dppcon.o dppequ.o dpprfs.o dppsv.o dppsvx.o \
	dpptrf.o dpptri.o dpptrs.o dpstf2.o dpstrf.o dptcon.o dpteqr.o \
	dptrfs.o dptsv.o dptsvx.o dpttrf.o dpttrs.o dptts2.o drscl.o \
	dsbevd.o dsbev.o dsbevx.o dsbgst.o dsbgvd.o dsbgv.o dsbgvx.o \
	dsbtrd.o dsfrk.o dsgesv.o dspcon.o dspevd.o dspev.o dspevx.o \
	dspgst.o dspgvd.o dspgv.o dspgvx.o dsposv.o dsprfs.o dspsv.o \
	dspsvx.o dsptrd.o dsptrf.o dsptri.o dsptrs.o dstebz.o dstedc.o \
	dstegr.o dstein.o dstemr.o dsteqr.o dsterf.o dstevd.o dstev.o \
	dstevr.o dstevx.o dsycon.o dsyequb.o dsyevd.o dsyev.o dsyevr.o \
	dsyevx.o dsygs2.o dsygst.o dsygvd.o dsygv.o dsygvx.o dsyrfs.o \
	dsysv.o dsysvx.o dsytd2.o dsytf2.o dsytrd.o \
	dsytrf.o dsytri.o dsytrs.o dtbcon.o dtbrfs.o dtbtrs.o dtfsm.o \
	dtftri.o dtfttp.o dtfttr.o dtgevc.o dtgex2.o dtgexc.o dtgsen.o \
	dtgsja.o dtgsna.o dtgsy2.o dtgsyl.o dtpcon.o dtprfs.o dtptri.o \
	dtptrs.o dtpttf.o dtpttr.o dtrcon.o dtrevc.o dtrexc.o dtrrfs.o \
	dtrsen.o dtrsna.o dtrsyl.o dtrti2.o dtrtri.o dtrtrs.o dtrttf.o \
	dtrttp.o dtzrqf.o dtzrzf.o dzsum1.o icmax1.o ieeeck.o ilaclc.o \
	ilaclr.o iladiag.o iladlc.o iladlr.o ilaenv.o ilaprec.o \
	ilaslc.o ilaslr.o ilatrans.o ilauplo.o ilaver.o ilazlc.o \
	ilazlr.o iparmq.o izmax1.o lsamen.o sbdsdc.o sbdsqr.o scsum1.o \
	sdisna.o sgbbrd.o sgbcon.o sgbequb.o sgbequ.o sgbrfs.o \
	sgbsv.o sgbsvx.o sgbtf2.o sgbtrf.o sgbtrs.o \
	sgebak.o sgebal.o sgebd2.o sgebrd.o sgecon.o sgeequb.o sgeequ.o \
	sgees.o sgeesx.o sgeev.o sgeevx.o sgegs.o sgegv.o sgehd2.o \
	sgehrd.o sgejsv.o sgelq2.o sgelqf.o sgelsd.o sgels.o sgelss.o \
	sgelsx.o sgelsy.o sgeql2.o sgeqlf.o sgeqp3.o sgeqpf.o sgeqr2.o \
	sgeqrf.o sgerfs.o sgerq2.o sgerqf.o sgesc2.o sgesdd.o \
	sgesvd.o sgesv.o sgesvj.o sgesvx.o sgetc2.o sgetf2.o \
	sgetrf.o sgetri.o sgetrs.o sggbak.o sggbal.o sgges.o sggesx.o \
	sggev.o sggevx.o sggglm.o sgghrd.o sgglse.o sggqrf.o sggrqf.o \
	sggsvd.o sggsvp.o sgsvj0.o sgsvj1.o sgtcon.o sgtrfs.o sgtsv.o \
	sgtsvx.o sgttrf.o sgttrs.o sgtts2.o shgeqz.o shsein.o shseqr.o \
	sisnan.o slabad.o slabrd.o slacn2.o slacon.o slacpy.o sladiv.o \
	slae2.o slaebz.o slaed0.o slaed1.o slaed2.o slaed3.o slaed4.o \
	slaed5.o slaed6.o slaed7.o slaed8.o slaed9.o slaeda.o slaein.o \
	slaev2.o slaexc.o slag2d.o slag2.o \
	slags2.o slagtf.o slagtm.o slagts.o \
	slagv2.o slahqr.o slahr2.o slahrd.o slaic1.o slaisnan.o \
	slaln2.o slals0.o slalsa.o slalsd.o slamrg.o \
	slaneg.o slangb.o slange.o slangt.o slanhs.o slansb.o slansf.o \
	slansp.o slanst.o slansy.o slantb.o slantp.o slantr.o slanv2.o \
	slapll.o slapmt.o \
	slapy2.o slapy3.o slaqgb.o slaqge.o slaqp2.o \
	slaqps.o slaqr0.o slaqr1.o slaqr2.o slaqr3.o slaqr4.o slaqr5.o \
	slaqsb.o slaqsp.o slaqsy.o slaqtr.o slar1v.o slar2v.o slarfb.o \
	slarf.o slarfg.o slarfp.o slarft.o slarfx.o slargv.o slarnv.o \
	slarra.o slarrb.o slarrc.o slarrd.o slarre.o \
	slarrf.o slarrj.o slarrk.o slarrr.o slarrv.o \
	slartg.o slartv.o slaruv.o slarzb.o slarz.o slarzt.o slas2.o \
	slascl.o slasd0.o slasd1.o slasd2.o slasd3.o \
	slasd4.o slasd5.o slasd6.o slasd7.o slasd8.o slasda.o slasdq.o \
	slasdt.o slaset.o slasq1.o slasq2.o slasq3.o slasq4.o slasq5.o \
	slasq6.o slasr.o slasrt.o slassq.o slasv2.o slaswp.o slasy2.o \
	slasyf.o slatbs.o slatdf.o slatps.o slatrd.o slatrs.o \
	slatrz.o slatzm.o slauu2.o slauum.o sopgtr.o \
	sopmtr.o sorg2l.o sorg2r.o sorgbr.o sorghr.o sorgl2.o sorglq.o \
	sorgql.o sorgqr.o sorgr2.o sorgrq.o sorgtr.o sorm2l.o sorm2r.o \
	sormbr.o sormhr.o sorml2.o sormlq.o sormql.o sormqr.o sormr2.o \
	sormr3.o sormrq.o sormrz.o sormtr.o spbcon.o spbequ.o spbrfs.o \
	spbstf.o spbsv.o spbsvx.o spbtf2.o spbtrf.o spbtrs.o spftrf.o \
	spftri.o spftrs.o spocon.o spoequb.o spoequ.o sporfs.o \
	sposv.o sposvx.o spotf2.o spotrf.o \
	spotri.o spotrs.o sppcon.o sppequ.o spprfs.o sppsv.o sppsvx.o \
	spptrf.o spptri.o spptrs.o spstf2.o spstrf.o sptcon.o spteqr.o \
	sptrfs.o sptsv.o sptsvx.o spttrf.o spttrs.o sptts2.o srscl.o \
	ssbevd.o ssbev.o ssbevx.o ssbgst.o ssbgvd.o ssbgv.o ssbgvx.o \
	ssbtrd.o ssfrk.o sspcon.o sspevd.o sspev.o sspevx.o sspgst.o \
	sspgvd.o sspgv.o sspgvx.o ssprfs.o sspsv.o sspsvx.o ssptrd.o \
	ssptrf.o ssptri.o ssptrs.o sstebz.o sstedc.o sstegr.o sstein.o \
	sstemr.o ssteqr.o ssterf.o sstevd.o sstev.o sstevr.o sstevx.o \
	ssycon.o ssyequb.o ssyevd.o ssyev.o ssyevr.o ssyevx.o ssygs2.o \
	ssygst.o ssygvd.o ssygv.o ssygvx.o ssyrfs.o ssysv.o \
	ssysvx.o ssytd2.o ssytf2.o ssytrd.o ssytrf.o \
	ssytri.o ssytrs.o stbcon.o stbrfs.o stbtrs.o stfsm.o stftri.o \
	stfttp.o stfttr.o stgevc.o stgex2.o stgexc.o stgsen.o stgsja.o \
	stgsna.o stgsy2.o stgsyl.o stpcon.o stprfs.o stptri.o stptrs.o \
	stpttf.o stpttr.o strcon.o strevc.o strexc.o strrfs.o strsen.o \
	strsna.o strsyl.o strti2.o strtri.o strtrs.o strttf.o strttp.o \
	stzrqf.o stzrzf.o xerbla_array.o xerbla.o zbdsqr.o zcgesv.o \
	zcposv.o zdrscl.o zgbbrd.o zgbcon.o zgbequb.o zgbequ.o \
	zgbrfs.o zgbsv.o zgbsvx.o zgbtf2.o \
	zgbtrf.o zgbtrs.o zgebak.o zgebal.o zgebd2.o zgebrd.o zgecon.o \
	zgeequb.o zgeequ.o zgees.o zgeesx.o zgeev.o zgeevx.o zgegs.o \
	zgegv.o zgehd2.o zgehrd.o zgelq2.o zgelqf.o zgelsd.o zgels.o \
	zgelss.o zgelsx.o zgelsy.o zgeql2.o zgeqlf.o zgeqp3.o \
	zgeqpf.o zgeqr2.o zgeqrf.o zgerfs.o zgerq2.o \
	zgerqf.o zgesc2.o zgesdd.o zgesvd.o zgesv.o zgesvx.o \
	zgetc2.o zgetf2.o zgetrf.o zgetri.o zgetrs.o zggbak.o zggbal.o \
	zgges.o zggesx.o zggev.o zggevx.o zggglm.o zgghrd.o zgglse.o \
	zggqrf.o zggrqf.o zggsvd.o zggsvp.o zgtcon.o zgtrfs.o zgtsv.o \
	zgtsvx.o zgttrf.o zgttrs.o zgtts2.o zhbevd.o zhbev.o zhbevx.o \
	zhbgst.o zhbgvd.o zhbgv.o zhbgvx.o zhbtrd.o zhecon.o zheequb.o \
	zheevd.o zheev.o zheevr.o zheevx.o zhegs2.o zhegst.o zhegvd.o \
	zhegv.o zhegvx.o zherfs.o zhesv.o zhesvx.o \
	zhetd2.o zhetf2.o zhetrd.o zhetrf.o zhetri.o zhetrs.o zhfrk.o \
	zhgeqz.o zhpcon.o zhpevd.o zhpev.o zhpevx.o zhpgst.o zhpgvd.o \
	zhpgv.o zhpgvx.o zhprfs.o zhpsv.o zhpsvx.o zhptrd.o zhptrf.o \
	zhptri.o zhptrs.o zhsein.o zhseqr.o zlabrd.o zlacgv.o zlacn2.o \
	zlacon.o zlacp2.o zlacpy.o zlacrm.o zlacrt.o zladiv.o zlaed0.o \
	zlaed7.o zlaed8.o zlaein.o zlaesy.o zlaev2.o zlag2c.o \
	zlags2.o \
	zlagtm.o zlahef.o \
	zlahqr.o zlahr2.o \
	zlahrd.o zlaic1.o zlals0.o zlalsa.o zlalsd.o \
	zlangb.o zlange.o zlangt.o zlanhb.o zlanhe.o zlanhf.o zlanhp.o \
	zlanhs.o zlanht.o zlansb.o zlansp.o zlansy.o zlantb.o zlantp.o \
	zlantr.o zlapll.o zlapmt.o \
	zlaqgb.o zlaqge.o \
	zlaqhb.o zlaqhe.o zlaqhp.o zlaqp2.o zlaqps.o zlaqr0.o zlaqr1.o \
	zlaqr2.o zlaqr3.o zlaqr4.o zlaqr5.o zlaqsb.o zlaqsp.o zlaqsy.o \
	zlar1v.o zlar2v.o zlarcm.o zlarfb.o zlarf.o zlarfg.o zlarfp.o \
	zlarft.o zlarfx.o zlargv.o zlarnv.o zlarrv.o \
	zlartg.o zlartv.o zlarzb.o zlarz.o zlarzt.o \
	zlascl.o zlaset.o zlasr.o zlassq.o zlaswp.o \
	zlasyf.o \
	zlat2c.o zlatbs.o \
	zlatdf.o zlatps.o zlatrd.o zlatrs.o zlatrz.o zlatzm.o zlauu2.o \
	zlauum.o zpbcon.o zpbequ.o zpbrfs.o zpbstf.o \
	zpbsv.o zpbsvx.o zpbtf2.o zpbtrf.o zpbtrs.o zpftrf.o zpftri.o \
	zpftrs.o zpocon.o zpoequb.o zpoequ.o zporfs.o \
	zposv.o zposvx.o zpotf2.o zpotrf.o zpotri.o zpotrs.o \
	zppcon.o zppequ.o zpprfs.o zppsv.o zppsvx.o zpptrf.o zpptri.o \
	zpptrs.o zpstf2.o zpstrf.o zptcon.o zpteqr.o zptrfs.o zptsv.o \
	zptsvx.o zpttrf.o zpttrs.o zptts2.o zrot.o zspcon.o zspmv.o \
	zspr.o zsprfs.o zspsv.o zspsvx.o zsptrf.o zsptri.o zsptrs.o \
	zstedc.o zstegr.o zstein.o zstemr.o zsteqr.o zsycon.o \
	zsyequb.o zsymv.o zsyr.o zsyrfs.o zsysv.o zsysvx.o \
	zsytf2.o zsytrf.o zsytri.o zsytrs.o ztbcon.o \
	ztbrfs.o ztbtrs.o ztfsm.o ztftri.o ztfttp.o ztfttr.o ztgevc.o \
	ztgex2.o ztgexc.o ztgsen.o ztgsja.o ztgsna.o ztgsy2.o ztgsyl.o \
	ztpcon.o ztprfs.o ztptri.o ztptrs.o ztpttf.o ztpttr.o ztrcon.o \
	ztrevc.o ztrexc.o ztrrfs.o ztrsen.o ztrsna.o ztrsyl.o ztrti2.o \
	ztrtri.o ztrtrs.o ztrttf.o ztrttp.o ztzrqf.o ztzrzf.o zung2l.o \
	zung2r.o zungbr.o zunghr.o zungl2.o zunglq.o zungql.o zungqr.o \
	zungr2.o zungrq.o zungtr.o zunm2l.o zunm2r.o zunmbr.o zunmhr.o \
	zunml2.o zunmlq.o zunmql.o zunmqr.o zunmr2.o zunmr3.o zunmrq.o \
	zunmrz.o zunmtr.o zupgtr.o zupmtr.o ../INSTALL/dlamch.o \
	../INSTALL/dsecnd_INT_ETIME.o ../INSTALL/lsame.o \
	../INSTALL/second_INT_ETIME.o ../INSTALL/slamch.o

SXLASRC = sgesvxx.o sgerfsx.o sla_gerfsx_extended.o sla_geamv.o		\
   sla_gercond.o sla_rpvgrw.o ssysvxx.o ssyrfsx.o			\
   sla_syrfsx_extended.o sla_syamv.o sla_syrcond.o sla_syrpvgrw.o	\
   sposvxx.o sporfsx.o sla_porfsx_extended.o sla_porcond.o		\
   sla_porpvgrw.o sgbsvxx.o sgbrfsx.o sla_gbrfsx_extended.o		\
   sla_gbamv.o sla_gbrcond.o sla_gbrpvgrw.o sla_lin_berr.o slarscl2.o	\
   slascl2.o sla_wwaddw.o

DXLASRC = dgesvxx.o dgerfsx.o dla_gerfsx_extended.o dla_geamv.o		\
   dla_gercond.o dla_rpvgrw.o dsysvxx.o dsyrfsx.o			\
   dla_syrfsx_extended.o dla_syamv.o dla_syrcond.o dla_syrpvgrw.o	\
   dposvxx.o dporfsx.o dla_porfsx_extended.o dla_porcond.o		\
   dla_porpvgrw.o dgbsvxx.o dgbrfsx.o dla_gbrfsx_extended.o		\
   dla_gbamv.o dla_gbrcond.o dla_gbrpvgrw.o dla_lin_berr.o dlarscl2.o	\
   dlascl2.o dla_wwaddw.o

CXLASRC = cgesvxx.o cgerfsx.o cla_gerfsx_extended.o cla_geamv.o		\
   cla_gercond_c.o cla_gercond_x.o cla_rpvgrw.o				\
   csysvxx.o csyrfsx.o cla_syrfsx_extended.o cla_syamv.o		\
   cla_syrcond_c.o cla_syrcond_x.o cla_syrpvgrw.o			\
   cposvxx.o cporfsx.o cla_porfsx_extended.o				\
   cla_porcond_c.o cla_porcond_x.o cla_porpvgrw.o			\
   cgbsvxx.o cgbrfsx.o cla_gbrfsx_extended.o cla_gbamv.o		\
   cla_gbrcond_c.o cla_gbrcond_x.o cla_gbrpvgrw.o			\
   chesvxx.o cherfsx.o cla_herfsx_extended.o cla_heamv.o		\
   cla_hercond_c.o cla_hercond_x.o cla_herpvgrw.o			\
   cla_lin_berr.o clarscl2.o clascl2.o cla_wwaddw.o

ZXLASRC = zgesvxx.o zgerfsx.o zla_gerfsx_extended.o zla_geamv.o		\
   zla_gercond_c.o zla_gercond_x.o zla_rpvgrw.o zsysvxx.o zsyrfsx.o	\
   zla_syrfsx_extended.o zla_syamv.o zla_syrcond_c.o zla_syrcond_x.o	\
   zla_syrpvgrw.o zposvxx.o zporfsx.o zla_porfsx_extended.o		\
   zla_porcond_c.o zla_porcond_x.o zla_porpvgrw.o zgbsvxx.o zgbrfsx.o	\
   zla_gbrfsx_extended.o zla_gbamv.o zla_gbrcond_c.o zla_gbrcond_x.o	\
   zla_gbrpvgrw.o zhesvxx.o zherfsx.o zla_herfsx_extended.o		\
   zla_heamv.o zla_hercond_c.o zla_hercond_x.o zla_herpvgrw.o		\
   zla_lin_berr.o zlarscl2.o zlascl2.o zla_wwaddw.o

ALLXAUX =

ifdef USEXBLAS
ALLXOBJ=$(SXLASRC) $(DXLASRC) $(CXLASRC) $(ZXLASRC) $(ALLXAUX)
endif

static: $(OBJS) $(ALLXOBJ)
	ar ruv liblapack.a $(OBJS)
	ranlib liblapack.a

shared: $(OBJS) $(ALLXOBJ)
	cc -L$(PREFIX)/lib $(CFLAGS) -shared -Wl,-rpath,some-path-that-does-not-exist-for-the-purpose-of-reserving-a-buffer-for-future-adjustments-to-the-rpath -Wl,-soname,liblapack-casatools.so.3 -o liblapack-casatools.so.3.2.1 $(OBJS) -lblas -lm -lgfortran -lc
	cp liblapack-casatools.so.3.2.1 $(PREFIX)/lib
	cd $(PREFIX)/lib && rm -f liblapack-casatools.so.3 && ln -s liblapack-casatools.so.3.2.1 liblapack-casatools.so.3
	cd $(PREFIX)/lib && rm -f liblapack.so && ln -s liblapack-casatools.so.3.2.1 liblapack.so
