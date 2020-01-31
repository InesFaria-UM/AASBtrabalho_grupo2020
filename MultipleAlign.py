
from MyAlign import MyAlign
from MySeq import MySeq
from AlignSeq import AlignSeq
from SubstMatrix import SubstMatrix

class MultipleAlign(object):

    def __init__(self, seqs, alignseq):
        self.seqs = seqs
        self.alignpars = alignseq
    
    def numSeqs(self):
        return len(self.seqs)
    
    def score_SP_col (self, charsCol):
        sc = 0
        for i in range(len(charsCol)-1):
            for j in range (i+1,len(charsCol)):
                sc+=self.alignpars.scorePos(charsCol[i],charsCol[j])       
        return sc
     
    def scoreSP (self, alinhamento):
        spi = 0
        for i in range (len(alinhamento)):
            carsCol= alinhamento.column(i)
            spi=self.score_SP_col(carsCol)
        return spi
    
    def addSeqAlignment (self, alignment, seq):
        res = []
        for i in range(len(alignment.listseqs)+1):
            res.append("")
        cons = MySeq(alignment.consensus(),alignment.tipo)
        self.alignpars.needlemanWunsch(cons, seq)
        align2 = self.alignpars.recoverAlignment()
        orig = 0
        for i in range(len(align2)):
            if align2[0,i]== '-':
                for k in range(len(alignment.listseqs)):
                    res[k] += "-"
            else:
                for k in range(len(alignment.listseqs)):
                    res[k] += alignment[k,orig]
                orig+=1
        res[len(alignment.listseqs)] = align2.listseqs[1]
        return MyAlign(res, alignment.tipo)
    
    def alignConsensus(self):
        self.alignpars.needlemanWunsch(self.seqs[0], self.seqs[1])
        res = self.alignpars.recoverAlignment()

        for i in range(2, len(self.seqs)):
            res = self.addSeqAlignment(res, self.seqs[i])
        return res

def printMat (mat):
    for i in range(0, len(mat)):
        print(mat[i])

def test():  
    s1 = MySeq("PHWAS","protein")
    s2 = MySeq("HWASW","protein")
    s3 = MySeq("HPHWA","protein")
    sm = SubstMatrix()
    sm.loadFromFile("blosum62.mat", "\t")
    aseq = AlignSeq(sm, -8)
    ma = MultipleAlign([s1,s2,s3], aseq)
    alinm = ma.alignConsensus()
    print(ma.scoreSP(alinm))
    print(alinm)
    

# -PHWAS-
# --HWASW
# HPHWA--

def test2():  
    s1 = MySeq("SWSSKLMKKIM","protein")
    s2 = MySeq("SYSLMKLKSWK","protein")
    s3 = MySeq("SWSSLMKLILS","protein")
    s4 = MySeq("SWSLMKLISSW","protein")
    sm = SubstMatrix()
    sm.loadFromFile("blosum62.mat", "\t")
    aseq = AlignSeq(sm, -8)
    ma = MultipleAlign([s1,s2,s3, s4], aseq)
    alinm = ma.alignConsensus()
    print(alinm) 
    print(ma.scoreSP(alinm))
    


#SWSSKLMKKIM--
#SYSLMKLKSWK--
#SWSS-LMKLILS-
#SWS--LMKLISSW
#30

def testBiom():
    #s1 = MySeq("NRCD","protein")   
    #s2 = MySeq("RADC","protein")
    #s3 = MySeq("NRDC","protein")
    
    s1 = MySeq("RADN","protein")   
    s2 = MySeq("ARCD","protein")
    s3 = MySeq("ARDN","protein")
        
    
    sm = SubstMatrix()
    sm.loadFromFile("blosum62.mat", "\t")
    aseq = AlignSeq(sm, -2)
    ma = MultipleAlign([s1,s2,s3], aseq)
    alin = MyAlign(["-RADN", "ARCD-", "AR-DN"], "protein")
    print(alin)
    print(ma.scoreSP(alin))

    alinm = ma.alignConsensus()
    print(alinm) 
    
    print(ma.highQuality(alin))

def testAASB():
    s1 = MySeq("SWSSKLMKKIM","protein")
    s2 = MySeq("SYSLMKLKSWK","protein")
    s3 = MySeq("SWSSLMKLILS","protein")
    s4 = MySeq("SWSLMKLISSW","protein")
    sm = SubstMatrix()
    sm.loadFromFile("blosum62.mat", "\t")
    aseq = AlignSeq(sm, -2)
    ma = MultipleAlign([s1,s2,s3, s4], aseq)

    alin = MyAlign(["SWSSKLMKKIM-", "SYSLMKLKSWK-", "SWS-SLMKLILS", "SWSL-MKLISSW"], "protein")   
    print(alin)
    print(ma.scoreSP(alin))
    
    x = ma.improveAlignment(alin)
    print(x)
    print(ma.scoreSP(x))

def testAASB2017 ():
    s1 = MySeq("ACATATCAT")
    s2 = MySeq("ACTAGATCT")
    s3 = MySeq("AGATATTAG")
    s4 = MySeq("GCATCGATT")
    
    sm = SubstMatrix()
    sm.createFromMatchPars(1,-1,"ACGT")
    aseq = AlignSeq(sm,-1)
    
    print(aseq.needlemanWunsch(s1, s2, True))
    printMat(aseq.S)
    printMat(aseq.T)
    print(aseq.recoverAlignment_with_ties())
    
    ma = MultipleAlign([s1,s2,s3,s4], aseq)
    al = ma.alignConsensus()
    print(al)
    al.distPairs()

def testBook ():
    s1 = MySeq("ATAGC")
    s2 = MySeq("AACC")
    s3 = MySeq("ATGAC")
    
    sm = SubstMatrix()
    sm.createFromMatchPars(1,-1,"ACGT")
    aseq = AlignSeq(sm,-1)
    ma = MultipleAlign([s1,s2,s3], aseq)
    al = ma.alignConsensus()
    print(al)
    
def testines ():
    s1 = MySeq("-AP-SC")
    s2 = MySeq("AACC")
    s3 = MySeq("ATGAC")
    
    sm = SubstMatrix()
    sm.createFromMatchPars(1,-1,"ACGT")
    aseq = AlignSeq(sm,-1)
    ma = MultipleAlign([s1,s2,s3], aseq)
    al = ma.alignConsensus()
    print(al)

if __name__ == "__main__": 
    test()
    #test2()
    #testBiom()
    testAASB2017()
    #test2()
    testBook()
