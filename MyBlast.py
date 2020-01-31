# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 16:55:47 2015

@author: miguelrocha
"""

class MyBlast:
    '''
    Classe para matrizes de pontos
    '''

    def __init__(self, filename = None, w = 3):
        '''
        Construtor
        '''
        if filename is not None:
            self.readDatabase(filename)
        else:
            self.db = []
        self.w = w
        self.map = None

    def readDatabase(self, db=None,filename=None):
        if filename ==None:
            self.db= db
        else:
            self.db=[]
            f=open(filename)
            for line in f:
                self.db.append(line.rstrip())
            f.close()
        
    def addSequenceDB(self, seq):
        self.db.append(seq)
        
    def buildMap (self, query):
        self.map = {}  
        for i in range(len(query)-self.w+1):
            subseq= query[i:i+self.w]
            if subseq in self.map:
                self.map[subseq].append(i)
            else:
                self.map[subseq]=[i]
                
            
        
    def getHits (self, seq, query):
        if self.map is None: 
            self.buildMap(query)
        res = []
        for i in range(len(query)-self.w+1):
            subseq= query[i:i+self.w]
            if subseq in self.map:
                I= self.map[subseq]
                for elem in I:
                    res.append((elem,i))
        return res
        
    def extendsHit (self, seq, hit, query):
        stq, sts = hit[0], hit[1]
        matfw = 0       
        k=0
        bestk = 0
        while 2*matfw >= k and stq+self.w+k < len(query) and sts+self.w+k < len(seq):
            if query[stq+self.w+k] == seq[sts+self.w+k]: 
                matfw+=1
                bestk = k+1
            k += 1
        size = self.w + bestk
    
        k = 0
        matbw = 0   
        bestk = 0
        while 2*matbw >= k and stq > k and sts > k:
            if query[stq-k-1] == seq[sts-k-1]: 
                matbw+=1
                bestk = k+1
            k+=1       
        size += bestk
        
        return (stq-bestk, sts-bestk, size, self.w+matfw+matbw)
        
    def hitBestScore(self, seq, query):
        hits = self.getHits(seq, query)
        bestScore = -1.0
        best = ()
        for h in hits:
            ext = self.extendsHit(seq, h, query)
            score = ext[3]
            if score > bestScore or (score== bestScore and ext[2] < best[2]):
                bestScore = score
                best = ext
        return best
 
 
    def bestAlignment (self, query):
        self.buildMap(query)
        bestScore = -1.0
        res = (0,0,0,0,0)
        for k in range(0,len(self.db)):
            bestSeq = self.hitBestScore(self.db[k], query)
            if bestSeq != ():
                score = bestSeq[3]  
                if score > bestScore or (score== bestScore and bestSeq[2] < res[2]):
                    bestScore = score
                    res = bestSeq[0], bestSeq[1], bestSeq[2], bestSeq[3], k
        if bestScore < 0: return ()
        else: return res
        
# def test1():
#     mb = MyBlast("seqBlast.txt", 11)
#     query = "gacgcctcgcgctcgcgcgctgaggcaaaaaaaaaaaaaaaaaaaatcggatagctagctgagcgctcgatagcgcgttcgctgcatcgcgtatagcgctgaagctcccggcgagctgtctgtaaatcggatctcatctcgctctatcct"
#     r = mb.bestAlignment(query)
#     print(r)

# # (1, 38, 149, 108, 3)

# def test2():
#     mb = MyBlast("seqBlast.txt", 11)
#     query2 = "cgacgacgacgacgaatgatg"
#     r = mb.bestAlignment(query2)
#     print(r)       

# #(0, 0, 21, 21, 4)

# test1()
# print()
# test2()
