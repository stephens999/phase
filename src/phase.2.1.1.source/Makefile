SHELL   = sh
CXX     = g++
INCLUDE = -I. -I..

#CFLAGS  = -g -O4 -Wall -W -Wno-sign-compare -pedantic -DDEBUG=5   $(INCLUDE)
#CFLAGS  = -g -O4 -Wall -W -Wno-sign-compare -pedantic -DDEBUG=5 -DCHECK_BOUNDARY  $(INCLUDE)
#CFLAGS  = -g -O4 -Wall -W -Wno-sign-compare -pedantic -DNDEBUG  $(INCLUDE)
#CFLAGS = -O2 -DNDEBUG -DTREAT_SNPS_AS_MS $(INCLUDE)
#CFLAGS = -O2 -DBIGDATASETS $(INCLUDE)b
#CFLAGS = -g -DBIGDATASETS $(INCLUDE)
#CFLAGS = -O2 $(INCLUDE)
#CFLAGS = -O2  $(INCLUDE)
#CFLAGS = -g -a -O2 -pg -DNDEBUG $(INCLUDE)

DEFAULTFLAGS =  -O2 -DBIGDATASETS $(INCLUDE)
#DEFAULTFLAGS =  -g -O2 -DBIGDATASETS -DDMALLOC $(INCLUDE)
#DEFAULTFLAGS =  -g -pg -O2 -DBIGDATASETS $(INCLUDE)
BIGFLAGS = -O2 -DBIGDATASETS $(INCLUDE)
SMLFLAGS = -O2 $(INCLUDE)
DEBUGFLAGS = -g -DDEBUG=5 -DCHECK_BOUNDARY $(INCLUDE)

#LDFLAGS = -L../crand/ -lcrand -lm -ldmallocxx
LDFLAGS = -lm
#LDFLAGS = /usr/lib/libm.a -static 

.SUFFIXES: 
.SUFFIXES: .cpp .o .c

.cpp.o:
	$(CXX) $(DEFAULTFLAGS) -c $<

testSNPHAPOBJ = testSNPHaplotype.o SNPHaplotype.o
testHAPOBJ = testHaplotype.o indnode.o utility.o arrayQ.o arrayCC.o arrayDiffProb.o arrayDiploidDiffProb.o arrayDiffCount.o Haplotype.o classpop.o HapList2.o arrayFF.o 
#OBJ = phase.o indnode.o utility.o arrayQ.o arrayFF.o arrayCC.o arrayDiffProb.o arrayDiploidDiffProb.o arrayDiffCount.o HapList2.o SNPHaplotype.o Haplotype.o MSHaplotype.o 
OBJ = phase.o indnode.o utility.o arrayQ.o arrayFF.o arrayCC.o arrayDiffProb.o arrayDiploidDiffProb.o arrayDiffCount.o HapList2.o Haplotype.o classpop.o recom.o HapPairList.o Summary.o mt19937ar.o pacb.lookup.o 
HDR = indnode.hpp classpop.hpp constants.hpp arrayQ.hpp utility.hpp arrayFF.hpp arrayCC.hpp arrayDiffProb.hpp arrayDiploidDiffProb.hpp arrayDiffCount.hpp HapList2.hpp 

PHASE : $(OBJ)
	$(CXX) $(DEFAULTFLAGS) -I. -I.. $(OBJ) -o $@  $(LDFLAGS)
big : clean
	@$(MAKE) PHASE.big "CFLAGS=$(BIGFLAGS)"
small : clean
	@$(MAKE) PHASE.small "CFLAGS=$(SMLFLAGS)"
debug: clean
	@$(MAKE) PHASE.big  "CFLAGS=$(DEBUGFLAGS)"

all: phase test

test: phase test1 test2 test4

test1:
	phase -n -f1 ../test/test2.inp test2.out test2z.out
test2:    
	phase -g1 ../test/test3.inp test3.out test3z.out
test3:    
	phase -H3 ../test/hudson.k1.t4.r40.20.in hudson.out hudsonz.out
test4:    
	phase -D3 ../test/test.inp test1.out test1z.out

testRPLbig:
	./PHASE -D54 -j46 -MR -T -R1 -n test/p2r200t20.100.inp test/p2r200t20.100.RPL2.out

testSPLbig:
	./PHASE -D100 -MS -T -R1 -n test/p2r200t20.100.inp test/p2r200t20.100.SPL.out

testQPLbig:
	./PHASE -D100 -MQ -T -R1 -n test/p2r200t20.100.inp test/p2r200t20.100.QPL.out

testSPLms:
	./PHASE -D100 -MS -D100 -T -R1 -n test/msr4t8.10.inp test/msr4t8.10.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/msr4t8.20.inp test/msr4t8.20.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/msr4t8.30.inp test/msr4t8.30.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/msr4t8.40.inp test/msr4t8.40.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/msr4t8.50.inp test/msr4t8.50.SPL.out

	perl perl/newcheck.pl test/msr4t8.10.inp test/msr4t8.10.SPL.out test/msr4t8.10.SPL.res
	perl perl/newcheck.pl test/msr4t8.20.inp test/msr4t8.20.SPL.out test/msr4t8.20.SPL.res
	perl perl/newcheck.pl test/msr4t8.30.inp test/msr4t8.30.SPL.out test/msr4t8.30.SPL.res
	perl perl/newcheck.pl test/msr4t8.40.inp test/msr4t8.40.SPL.out test/msr4t8.40.SPL.res
	perl perl/newcheck.pl test/msr4t8.50.inp test/msr4t8.50.SPL.out test/msr4t8.50.SPL.res

	perl perl/summarise.pl test/msr4t8.10.SPL.res >  test/SPLms.sum	
	perl perl/summarise.pl test/msr4t8.20.SPL.res >> test/SPLms.sum	
	perl perl/summarise.pl test/msr4t8.30.SPL.res >> test/SPLms.sum	
	perl perl/summarise.pl test/msr4t8.40.SPL.res >> test/SPLms.sum	
	perl perl/summarise.pl test/msr4t8.50.SPL.res >> test/SPLms.sum	

testRPLms:
	./PHASE -D100 -MR -D100 -T -R1 -n test/msr4t8.10.inp test/msr4t8.10.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/msr4t8.20.inp test/msr4t8.20.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/msr4t8.30.inp test/msr4t8.30.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/msr4t8.40.inp test/msr4t8.40.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/msr4t8.50.inp test/msr4t8.50.RPL.out

	perl perl/newcheck.pl test/msr4t8.10.inp test/msr4t8.10.RPL.out test/msr4t8.10.RPL.res
	perl perl/newcheck.pl test/msr4t8.20.inp test/msr4t8.20.RPL.out test/msr4t8.20.RPL.res
	perl perl/newcheck.pl test/msr4t8.30.inp test/msr4t8.30.RPL.out test/msr4t8.30.RPL.res
	perl perl/newcheck.pl test/msr4t8.40.inp test/msr4t8.40.RPL.out test/msr4t8.40.RPL.res
	perl perl/newcheck.pl test/msr4t8.50.inp test/msr4t8.50.RPL.out test/msr4t8.50.RPL.res

	perl perl/summarise.pl test/msr4t8.10.RPL.res >  test/RPLms.sum	
	perl perl/summarise.pl test/msr4t8.20.RPL.res >> test/RPLms.sum	
	perl perl/summarise.pl test/msr4t8.30.RPL.res >> test/RPLms.sum	
	perl perl/summarise.pl test/msr4t8.40.RPL.res >> test/RPLms.sum	
	perl perl/summarise.pl test/msr4t8.50.RPL.res >> test/RPLms.sum	

testQPLms:
	./PHASE -D100 -MQ -D100 -T -R1 -n test/msr4t8.10.inp test/msr4t8.10.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/msr4t8.20.inp test/msr4t8.20.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/msr4t8.30.inp test/msr4t8.30.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/msr4t8.40.inp test/msr4t8.40.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/msr4t8.50.inp test/msr4t8.50.QPL.out

	perl perl/newcheck.pl test/msr4t8.10.inp test/msr4t8.10.QPL.out test/msr4t8.10.QPL.res
	perl perl/newcheck.pl test/msr4t8.20.inp test/msr4t8.20.QPL.out test/msr4t8.20.QPL.res
	perl perl/newcheck.pl test/msr4t8.30.inp test/msr4t8.30.QPL.out test/msr4t8.30.QPL.res
	perl perl/newcheck.pl test/msr4t8.40.inp test/msr4t8.40.QPL.out test/msr4t8.40.QPL.res
	perl perl/newcheck.pl test/msr4t8.50.inp test/msr4t8.50.QPL.out test/msr4t8.50.QPL.res

	perl perl/summarise.pl test/msr4t8.10.QPL.res >  test/QPLms.sum	
	perl perl/summarise.pl test/msr4t8.20.QPL.res >> test/QPLms.sum	
	perl perl/summarise.pl test/msr4t8.30.QPL.res >> test/QPLms.sum	
	perl perl/summarise.pl test/msr4t8.40.QPL.res >> test/QPLms.sum	
	perl perl/summarise.pl test/msr4t8.50.QPL.res >> test/QPLms.sum	

testSPL:
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40.20.inp test/p2r40.20.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40.40.inp test/p2r40.40.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40.60.inp test/p2r40.60.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40.80.inp test/p2r40.80.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40.100.inp test/p2r40.100.SPL.out
	perl perl/newcheck.pl test/p2r40.20.inp test/p2r40.20.SPL.out test/p2r40.20.SPL.res     
	perl perl/newcheck.pl test/p2r40.40.inp test/p2r40.40.SPL.out test/p2r40.40.SPL.res     
	perl perl/newcheck.pl test/p2r40.60.inp test/p2r40.60.SPL.out test/p2r40.60.SPL.res     
	perl perl/newcheck.pl test/p2r40.80.inp test/p2r40.80.SPL.out test/p2r40.80.SPL.res     
	perl perl/newcheck.pl test/p2r40.100.inp test/p2r40.100.SPL.out test/p2r40.100.SPL.res     	
	perl perl/summarise.pl test/p2r40.20.SPL.res >  test/SPL.sum
	perl perl/summarise.pl test/p2r40.40.SPL.res >> test/SPL.sum   
	perl perl/summarise.pl test/p2r40.60.SPL.res >> test/SPL.sum   
	perl perl/summarise.pl test/p2r40.80.SPL.res >> test/SPL.sum   
	perl perl/summarise.pl test/p2r40.100.SPL.res >> test/SPL.sum	

testSPLm:
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40m05.20.inp test/p2r40m05.20.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40m05.40.inp test/p2r40m05.40.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40m05.60.inp test/p2r40m05.60.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40m05.80.inp test/p2r40m05.80.SPL.out
	./PHASE -D100 -MS -D100 -T -R1 -n test/p2r40m05.100.inp test/p2r40m05.100.SPL.out
	perl perl/newcheck.pl test/p2r40m05.20.inp test/p2r40m05.20.SPL.out test/p2r40m05.20.SPL.res     
	perl perl/newcheck.pl test/p2r40m05.40.inp test/p2r40m05.40.SPL.out test/p2r40m05.40.SPL.res     
	perl perl/newcheck.pl test/p2r40m05.60.inp test/p2r40m05.60.SPL.out test/p2r40m05.60.SPL.res     
	perl perl/newcheck.pl test/p2r40m05.80.inp test/p2r40m05.80.SPL.out test/p2r40m05.80.SPL.res     
	perl perl/newcheck.pl test/p2r40m05.100.inp test/p2r40m05.100.SPL.out test/p2r40m05.100.SPL.res     
	perl perl/summarise.pl test/p2r40m05.20.SPL.res > test/SPLm.sum
	perl perl/summarise.pl test/p2r40m05.40.SPL.res >> test/SPLm.sum   
	perl perl/summarise.pl test/p2r40m05.60.SPL.res >> test/SPLm.sum   
	perl perl/summarise.pl test/p2r40m05.80.SPL.res >> test/SPLm.sum   
	perl perl/summarise.pl test/p2r40m05.100.SPL.res >> test/SPLm.sum	

testRPL:
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40.20.inp test/p2r40.20.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40.40.inp test/p2r40.40.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40.60.inp test/p2r40.60.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40.80.inp test/p2r40.80.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40.100.inp test/p2r40.100.RPL.out
	perl perl/newcheck.pl test/p2r40.20.inp test/p2r40.20.RPL.out test/p2r40.20.RPL.res     
	perl perl/newcheck.pl test/p2r40.40.inp test/p2r40.40.RPL.out test/p2r40.40.RPL.res     
	perl perl/newcheck.pl test/p2r40.60.inp test/p2r40.60.RPL.out test/p2r40.60.RPL.res     
	perl perl/newcheck.pl test/p2r40.80.inp test/p2r40.80.RPL.out test/p2r40.80.RPL.res     
	perl perl/newcheck.pl test/p2r40.100.inp test/p2r40.100.RPL.out test/p2r40.100.RPL.res    
	perl perl/summarise.pl test/p2r40.20.RPL.res >  test/RPL.sum
	perl perl/summarise.pl test/p2r40.40.RPL.res >> test/RPL.sum   
	perl perl/summarise.pl test/p2r40.60.RPL.res >> test/RPL.sum   
	perl perl/summarise.pl test/p2r40.80.RPL.res >> test/RPL.sum   
	perl perl/summarise.pl test/p2r40.100.RPL.res >> test/RPL.sum	

testRPLm:
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40m05.20.inp test/p2r40m05.20.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40m05.40.inp test/p2r40m05.40.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40m05.60.inp test/p2r40m05.60.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40m05.80.inp test/p2r40m05.80.RPL.out
	./PHASE -D100 -MR -D100 -T -R1 -n test/p2r40m05.100.inp test/p2r40m05.100.RPL.out
	perl perl/newcheck.pl test/p2r40m05.20.inp test/p2r40m05.20.RPL.out test/p2r40m05.20.RPL.res     
	perl perl/newcheck.pl test/p2r40m05.40.inp test/p2r40m05.40.RPL.out test/p2r40m05.40.RPL.res     
	perl perl/newcheck.pl test/p2r40m05.60.inp test/p2r40m05.60.RPL.out test/p2r40m05.60.RPL.res     
	perl perl/newcheck.pl test/p2r40m05.80.inp test/p2r40m05.80.RPL.out test/p2r40m05.80.RPL.res     
	perl perl/newcheck.pl test/p2r40m05.100.inp test/p2r40m05.100.RPL.out test/p2r40m05.100.RPL.res     
	perl perl/summarise.pl test/p2r40m05.20.RPL.res > test/RPLm.sum
	perl perl/summarise.pl test/p2r40m05.40.RPL.res >> test/RPLm.sum   
	perl perl/summarise.pl test/p2r40m05.60.RPL.res >> test/RPLm.sum   
	perl perl/summarise.pl test/p2r40m05.80.RPL.res >> test/RPLm.sum   
	perl perl/summarise.pl test/p2r40m05.100.RPL.res >> test/RPLm.sum	

testQPL:
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40.20.inp test/p2r40.20.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40.40.inp test/p2r40.40.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40.60.inp test/p2r40.60.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40.80.inp test/p2r40.80.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40.100.inp test/p2r40.100.QPL.out
	perl perl/newcheck.pl test/p2r40.20.inp test/p2r40.20.QPL.out test/p2r40.20.QPL.res     
	perl perl/newcheck.pl test/p2r40.40.inp test/p2r40.40.QPL.out test/p2r40.40.QPL.res     
	perl perl/newcheck.pl test/p2r40.60.inp test/p2r40.60.QPL.out test/p2r40.60.QPL.res     
	perl perl/newcheck.pl test/p2r40.80.inp test/p2r40.80.QPL.out test/p2r40.80.QPL.res     
	perl perl/newcheck.pl test/p2r40.100.inp test/p2r40.100.QPL.out test/p2r40.100.QPL.res     
	perl perl/summarise.pl test/p2r40.20.QPL.res > test/QPL.sum
	perl perl/summarise.pl test/p2r40.40.QPL.res >> test/QPL.sum   
	perl perl/summarise.pl test/p2r40.60.QPL.res >> test/QPL.sum   
	perl perl/summarise.pl test/p2r40.80.QPL.res >> test/QPL.sum   
	perl perl/summarise.pl test/p2r40.100.QPL.res >> test/QPL.sum	


testQPLm:
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40m05.20.inp test/p2r40m05.20.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40m05.40.inp test/p2r40m05.40.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40m05.60.inp test/p2r40m05.60.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40m05.80.inp test/p2r40m05.80.QPL.out
	./PHASE -D100 -MQ -D100 -T -R1 -n test/p2r40m05.100.inp test/p2r40m05.100.QPL.out
	perl perl/newcheck.pl test/p2r40m05.20.inp test/p2r40m05.20.QPL.out test/p2r40m05.20.QPL.res     
	perl perl/newcheck.pl test/p2r40m05.40.inp test/p2r40m05.40.QPL.out test/p2r40m05.40.QPL.res     
	perl perl/newcheck.pl test/p2r40m05.60.inp test/p2r40m05.60.QPL.out test/p2r40m05.60.QPL.res     
	perl perl/newcheck.pl test/p2r40m05.80.inp test/p2r40m05.80.QPL.out test/p2r40m05.80.QPL.res     
	perl perl/newcheck.pl test/p2r40m05.100.inp test/p2r40m05.100.QPL.out test/p2r40m05.100.QPL.res     
	perl perl/summarise.pl test/p2r40m05.20.QPL.res > test/QPLm.sum
	perl perl/summarise.pl test/p2r40m05.40.QPL.res >> test/QPLm.sum   
	perl perl/summarise.pl test/p2r40m05.60.QPL.res >> test/QPLm.sum   
	perl perl/summarise.pl test/p2r40m05.80.QPL.res >> test/QPLm.sum   
	perl perl/summarise.pl test/p2r40m05.100.QPL.res >> test/QPLm.sum	
testSNPHaplotype : $(testSNPHAPOBJ)
	$(CXX) $(DEFAULTFLAGS) -I. -I.. -o $@ $^ $(LDFLAGS)

testHaplotype : $(testHAPOBJ)
	$(CXX) $(DEFAULTFLAGS) -I. -I.. -o $@ $^ $(LDFLAGS)


PHASE.big : $(OBJ)
	$(CXX) $(BIGFLAGS) -I. -I.. -o $@ $^ $(LDFLAGS)

PHASE.small :  $(OBJ)
	$(CXX) $(SMLFLAGS) -o $@ $^ $(LDFLAGS)

postprocess: postprocess.o indnode.o classpop.o utility.o arrayQ.o arrayFF.o arrayCC.o summary.o arrayCC.o arrayDiffProb.o arrayDiploidDiffProb.o arrayDiffCount.o
	$(CXX) $(CFLAGS) -o $@ $^ $(LDFLAGS)

postprocess.o : postprocess.cpp $(HDR)
summary.o: summary.cpp summary.hpp

phase.o   : phase.cpp $(HDR)
indnode.o : indnode.cpp indnode.hpp constants.hpp Haplotype.cpp Haplotype.hpp  
classpop.o : classpop.cpp classpop.hpp constants.hpp arrayQ.hpp arrayCC.cpp arrayCC.hpp arrayDiffProb.hpp arrayDiffProb.cpp indnode.cpp indnode.hpp Haplotype.cpp Haplotype.hpp HapList2.cpp HapList2.hpp utility.cpp utility.hpp
utility.o  : utility.cpp utility.hpp constants.hpp
mt19937ar.o : mt19937ar.c
pacb.lookup.o : pacb.lookup.cpp pacb.lookup.hpp
arrayQ.o   : arrayQ.cpp arrayQ.hpp
arrayFF.o  : arrayFF.cpp arrayFF.hpp arrayQ.hpp arrayCC.cpp arrayCC.hpp arrayDiffProb.hpp arrayDiffProb.cpp Haplotype.cpp Haplotype.hpp
arrayCC.o  : arrayCC.cpp arrayFF.hpp arrayQ.hpp
arrayDiffProb.o : arrayDiffProb.cpp arrayDiffProb.hpp arrayQ.hpp arrayQ.cpp
arrayDiploidDiffProb.o : arrayDiploidDiffProb.cpp arrayDiploidDiffProb.hpp arrayQ.hpp arrayQ.cpp
arrayDiffCount.o : arrayDiffCount.cpp arrayDiffCount.hpp arrayDiffProb.cpp arrayDiffProb.hpp arrayCC.cpp arrayCC.hpp
recom.o : recom.cpp classpop.hpp
Haplotype.o: Haplotype.cpp Haplotype.hpp 
HapList2.o: HapList2.cpp HapList2.hpp Haplotype.cpp Haplotype.hpp pacb.lookup.cpp pacb.lookup.hpp
HapPairList.o: HapPairList.cpp HapPairList.hpp Haplotype.cpp Haplotype.hpp Summary.cpp Summary.hpp
Summary.o: Summary.cpp Summary.hpp Haplotype.cpp Haplotype.hpp


TAGS: *.hpp *.cpp
	etags $^
tags: *.hpp *.cpp
	ctags $^

# Stripping off the symbols can make the executable smaller.
strip: phase
	strip $^

clean:
	rm -rf *.o
ultraclean: clean
	rm -f phase
