
from MatrixNum import MatrixNum
from ClustHier import ClustHier
from MySeq import MySeq
from AlignSeq import AlignSeq
from SubstMatrix import SubstMatrix

class UPGMA:


    def __init__(self, seqs, alseq):
        self.seqs = seqs
        self.alseq = alseq
        self.matdist = MatrixNum(len(seqs), len(seqs))
        self.criaMatDists()
        
    def criaMatDists(self):
        for i in range(len(self.seqs)): 
            self.matdist.setValue(i, i, 0.0)
        for i in range(len(self.seqs)):
            for j in range(i, len(self.seqs)):
                s1 = self.seqs[i]
                s2 = self.seqs[j]
                self.alseq.needlemanWunsch(s1,s2)
                alin=self.alseq.recoverAlignment()
                ncd=0
                for k in range(len(alin)):
                    col=alin.column(k)
                    if [col][0] != col[1]:
                        ncd+=1 
                self.matdist.setValue(i,j,ncd)
                self.matdist.setValue(j,i,ncd)

    def run(self):
        ch = ClustHier(self.matdist)
        arv= ch.executeClustering()
        return arv
        

def test1():
    seq1 = MySeq("ATAGCGAT")    
    seq2 = MySeq("ATAGGCCT")    
    seq3 = MySeq("CTAGGCCC")    
    sm = SubstMatrix()    
    sm.createFromMatchPars(3, -1, "ACGT")    
    alseq = AlignSeq(sm, -8)    
    up  = UPGMA([seq1, seq2, seq3], alseq)    
    arv = up.run()    
    arv.printtree() 
    
  
# Root:[2, 1, 0] Dist.: 2.0
#        Left:[2, 1] Dist.: 1.0
#                Left:[2] Dist.: 0
#                Right:[1] Dist.: 0
#        Right:[0] Dist.: 0    
    
    
    
def test2():
    seq1 = MySeq("AGTAG")
    seq2 = MySeq("TATA")
    seq3 = MySeq("TAT")
    seq4 = MySeq("TAGTG")    
    sm = SubstMatrix()
    sm.createFromMatchPars(0, -1, "ACGT")
    alseq = AlignSeq(sm, -1)
    up  = UPGMA([seq1, seq2, seq3, seq4], alseq)
    up.matdist.printmat()
    arv = up.run()
    arv.printtree()

#[0, 3, 3, 2]
#[3, 0, 1, 2]
#[3, 1, 0, 2]
#[2, 2, 2, 0]
#
#Root:[3, 0, 2, 1] Dist.: 1.25
#        Left:[3, 0] Dist.: 1.0
#                Left:[3] Dist.: 0
#                Right:[0] Dist.: 0
#        Right:[2, 1] Dist.: 0.5
#                Left:[2] Dist.: 0
#                Right:[1] Dist.: 0


if __name__ == '__main__': 
    test1()
    print()
    test2()
