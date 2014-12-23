INSTALLDIR = ./.. # /usr/local/bin/molphy

CC = gcc # cc, gcc, c89
CFLAGS = -O # -O4 or +O
DEFINE = #
LIBS = -lm # -lM
CP = cp # copy command
NUCLE = -DNUC
NJLE = -DNJ
SDLE = -DSD

PRGS = protml protst nucml nucst njdist totalml
# prottpm nuctpm tmlprot protdml nucdml tridist sddist

OBJML = getseq.o seqproc.o mltree.o altree.o qltree.o sltree.o \
	njtree.o njmtree.o prtree.o pstree.o matrixut.o mygetopt.o
OBJPML = protml.o protproc.o dyhfjtt.o mtrev24.o tranprb.o distan.o mlklhd.o
OBJNML = Nucml.o nucproc.o Ntranprb.o Ndistan.o Nmlklhd.o

OBJST = getseq.o seqproc.o matrixut.o mygetopt.o
OBJPST = protst.o seqstat.o protproc.o
OBJNST = Nucst.o Nseqstat.o nucproc.o

OBJD = distproc.o Dpstree.o matrixut.o mygetopt.o
OBJTD = tridist.o triproc.o
OBJND = Njdist.o njproc.o
OBJSD = sddist.o sdproc.o

OBJTO = totalml.o matrixut.o mygetopt.o

OBJPTPM = Prottpm.o tpmproc.o protproc.o dyhfjtt.o mtrev24.o Trnprbt.o \
	distan.o mlklhd.o
OBJNTPM = Nuctpm.o Ntpmproc.o nucproc.o Ntrnprbt.o Ndistan.o Nmlklhd.o

OBJTML = getseq.o seqproc.o triadml.o mltree.o prtree.o matrixut.o mygetopt.o
OBJTMLP = tmlprot.o protproc.o dyhfjtt.o mtrev24.o tranprb.o distan.o mlklhd.o

OBJPDML = protdml.o dmlproc.o protproc.o dyhfjtt.o mtrev24.o tranprb.o \
	distan.o mlklhd.o
OBJNDML = Nucdml.o Ndmlproc.o nucproc.o Ntranprb.o Ndistan.o Nmlklhd.o

all : $(PRGS)

protml : $(OBJML) $(OBJPML)
	$(CC) $(CFLAGS) $(DEFINE) -o $@ $(OBJML) $(OBJPML) $(LIBS)

nucml : $(OBJML) $(OBJNML)
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -o $@ $(OBJML) $(OBJNML) $(LIBS)

protst : $(OBJST) $(OBJPST)
	$(CC) $(CFLAGS) $(DEFINE) -o $@ $(OBJST) $(OBJPST) $(LIBS)

nucst : $(OBJST) $(OBJNST)
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -o $@ $(OBJST) $(OBJNST) $(LIBS)

tridist : $(OBJD) $(OBJTD)
	$(CC) $(CFLAGS) $(DEFINE) -o $@ $(OBJD) $(OBJTD) $(LIBS)

njdist : $(OBJD) $(OBJND)
	$(CC) $(CFLAGS) $(DEFINE) -o $@ $(OBJD) $(OBJND) $(LIBS)

totalml : $(OBJTO)
	$(CC) $(CFLAGS) $(DEFINE) -o $@ $(OBJTO) $(LIBS)

sddist : $(OBJD) $(OBJSD)
	$(CC) $(CFLAGS) $(DEFINE) -o $@ $(OBJD) $(OBJSD) $(LIBS)

prottpm : $(OBJML) $(OBJPTPM)
	$(CC) $(CFLAGS) $(DEFINE) -o $@ $(OBJML) $(OBJPTPM) $(LIBS)

nuctpm : $(OBJML) $(OBJNTPM)
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -o $@ $(OBJML) $(OBJNTPM) $(LIBS)

tmlprot : $(OBJTML) $(OBJTMLP)
	$(CC) $(CFLAGS) $(DEFINE) -o $@ $(OBJTML) $(OBJTMLP) $(LIBS)

protdml : $(OBJML) $(OBJPDML)
	$(CC) $(CFLAGS) $(DEFINE) -o $@ $(OBJML) $(OBJPDML) $(LIBS)

nucdml : $(OBJML) $(OBJNDML)
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -o $@ $(OBJML) $(OBJNDML) $(LIBS)

# protml nucml

protml.o : protml.c protml.h molphy.h 
	$(CC) $(CFLAGS) $(DEFINE) -c protml.c
Nucml.c : protml.c optimtpm.c abratio.c
	$(CP) protml.c Nucml.c
Nucml.o : Nucml.c protml.h molphy.h optimtpm.c abratio.c
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -c Nucml.c
getseq.o : getseq.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c getseq.c
seqproc.o : seqproc.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c seqproc.c
tranprb.o : tranprb.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c tranprb.c
Ntranprb.c : tranprb.c
	$(CP) tranprb.c Ntranprb.c
Ntranprb.o : Ntranprb.c protml.h molphy.h
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -c Ntranprb.c
distan.o : distan.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c distan.c
Ndistan.c : distan.c
	$(CP) distan.c Ndistan.c
Ndistan.o : Ndistan.c protml.h molphy.h
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -c Ndistan.c
mltree.o : mltree.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c mltree.c
altree.o : altree.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c altree.c
qltree.o : qltree.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c qltree.c
sltree.o : sltree.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c sltree.c
njtree.o : njtree.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c njtree.c
njmtree.o : njmtree.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c njmtree.c
mlklhd.o : mlklhd.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c mlklhd.c
Nmlklhd.c : mlklhd.c
	$(CP) mlklhd.c Nmlklhd.c
Nmlklhd.o : Nmlklhd.c protml.h molphy.h
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -c Nmlklhd.c
prtree.o : prtree.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c prtree.c
pstree.o : pstree.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c pstree.c

protproc.o : protproc.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c protproc.c
nucproc.o : nucproc.c protml.h molphy.h
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -c nucproc.c
dyhfjtt.o : dyhfjtt.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c dyhfjtt.c
mtrev24.o : mtrev24.c protml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c mtrev24.c

# protst nucst

protst.o : protst.c protst.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c protst.c
Nucst.c : protst.c
	$(CP) protst.c Nucst.c
Nucst.o : Nucst.c protst.h molphy.h
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -c Nucst.c

seqstat.o : seqstat.c protst.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c seqstat.c
Nseqstat.c : seqstat.c
	$(CP) seqstat.c Nseqstat.c
Nseqstat.o : Nseqstat.c protst.h molphy.h
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -c Nseqstat.c

# tridist njdist

tridist.o : tridist.c tridist.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c tridist.c
triproc.o : triproc.c tridist.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c triproc.c
Njdist.c : tridist.c
	$(CP) tridist.c Njdist.c
Njdist.o : Njdist.c tridist.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) $(NJLE) -c Njdist.c
njproc.o : njproc.c tridist.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c njproc.c
sddist.o : sddist.c sddist.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) $(SDLE) -c sddist.c
sdproc.o : sdproc.c sddist.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c sdproc.c
distproc.o : distproc.c tridist.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c distproc.c
Dpstree.c : pstree.c
	$(CP) pstree.c Dpstree.c
Dpstree.o : Dpstree.c tridist.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) $(NJLE) -c Dpstree.c

# totalml

totalml.o : totalml.c totalml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c totalml.c

# tmlprot

tmlprot.o : tmlprot.c tmlprot.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c tmlprot.c

triadml.o : triadml.c tmlprot.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c triadml.c

# prottpm nuctpm

Prottpm.c : protml.c
	$(CP) protml.c Prottpm.c
Prottpm.o : Prottpm.c protml.h molphy.h
	$(CC) $(CFLAGS) -DTPM $(DEFINE) -c Prottpm.c
Nuctpm.c : protml.c
	$(CP) protml.c Nuctpm.c
Nuctpm.o : Nuctpm.c protml.h molphy.h
	$(CC) $(CFLAGS) -DTPM $(NUCLE) $(DEFINE) -c Nuctpm.c

tpmproc.o : tpmproc.c protml.h molphy.h
	$(CC) $(CFLAGS) -DTPM $(DEFINE) -c tpmproc.c
Ntpmproc.c : tpmproc.c
	$(CP) tpmproc.c Ntpmproc.c
Ntpmproc.o : Ntpmproc.c protml.h molphy.h
	$(CC) $(CFLAGS) -DTPM $(NUCLE) $(DEFINE) -c Ntpmproc.c

Trnprbt.c : tranprb.c
	$(CP) tranprb.c Trnprbt.c
Trnprbt.o : Trnprbt.c protml.h molphy.h
	$(CC) $(CFLAGS) -DTPM $(DEFINE) -c Trnprbt.c
Ntrnprbt.c : tranprb.c
	$(CP) tranprb.c Ntrnprbt.c
Ntrnprbt.o : Ntrnprbt.c protml.h molphy.h
	$(CC) $(CFLAGS) -DTPM $(NUCLE) $(DEFINE) -c Ntrnprbt.c

# protdml nucdml

protdml.o : protdml.c protdml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c protdml.c
Nucdml.c : protdml.c
	$(CP) protdml.c Nucdml.c
Nucdml.o : Nucdml.c protdml.h molphy.h optimtpm.c abratio.c
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -c Nucdml.c
dmlproc.o : dmlproc.c protdml.h molphy.h
	$(CC) $(CFLAGS) $(DEFINE) -c dmlproc.c
Ndmlproc.c : dmlproc.c
	$(CP) dmlproc.c Ndmlproc.c
Ndmlproc.o : Ndmlproc.c protdml.h molphy.h optimtpm.c abratio.c
	$(CC) $(CFLAGS) $(NUCLE) $(DEFINE) -c Ndmlproc.c

# others

matrixut.o : matrixut.c matrixut.h
	$(CC) $(CFLAGS) $(DEFINE) -c matrixut.c
mygetopt.o : mygetopt.c
	$(CC) $(CFLAGS) $(DEFINE) -c mygetopt.c

install : $(PRGS)
	$(CP) $(PRGS) $(INSTALLDIR)
	cd $(INSTALLDIR); \
	chmod 755 $(PRGS)

clobber :
	rm -f *.o $(PRGS)
	cd $(INSTALLDIR); \
	rm $(PRGS)

clean :
	rm -f core *.o $(PRGS)

backup :
	$(CP) [a-z]*.c *.h Makefile ./srcback

public :
	$(CP) [a-z]*.c *.h /users/adachi/work/molphy/molphy-2.3b3/src
	$(CP) Makefile /users/adachi/work/molphy/molphy-2.3b3/src/Makefile.org
	$(CP) ../utl/*.pl /users/adachi/work/molphy/molphy-2.3b3/utl
	$(CP) ../utl/Install /users/adachi/work/molphy/molphy-2.3b3/utl

backupevolu2 :
	rcp [a-z]*.c *.h Makefile evolu2:/users/adachi/work/molphyback

backupsunmh :
	rcp [a-z]*.c *.h Makefile sunmh:/home/sunmh/adachi/work/molphyback

installevolu2 :
	rcp [a-z]*.c *.h Makefile evolu2:/users/adachi/work/molphybeta

installsunmh :
	rcp [a-z]*.c *.h Makefile sunmh:/home/sunmh/adachi/work/molphybeta

RMSRC :
	rm bdate.c compcrit.c dmlproc.c protdml.c protdml.h ptnlkl.c \
		sddist.c sddist.h sdproc.c tmlprot.c tmlprot.h tpgraph.c tpmproc.c \
		triadml.c
