
from MyAlign import MyAlign
from MySeq import MySeq
from SubstMatrix import SubstMatrix
from MatrixNum import MatrixNum

class AlignSeq:

    def __init__(self, sm, g):
        self.g = g
        self.sm = sm
        self.S = None
        self.T = None
        self.seq1 = None
        self.seq2 = None
        
    def scorePos (self, c1, c2):
        if c1 == "-" or c2=="-":
            return self.g
        else:
            return self.sm[c1,c2]
        
    def scoreAlin (self, alin):
        res = 0
        for i in range(len(alin)):
            res += self.scorePos (alin[0][i], alin[1][i])
        return res
    
    def needlemanWunsch (self, seq1, seq2):
        if (seq1.tipo != seq2.tipo): return None
        self.S = MatrixNum(len(seq1)+1, len(seq2)+1)
        self.T = MatrixNum(len(seq1)+1, len(seq2)+1)
        self.seq1 = seq1
        self.seq2 = seq2
        for j in range(1, len(seq2)+1):
            self.S[0, j] = self.g * j
            self.T[0, j] = 3
        for i in range(1, len(seq1)+1):
            self.S[i, 0] = self.g * i
            self.T[i, 0] = 2
        for i in range(0, len(seq1)):
            for j in range(len(seq2)):
                s1 = self.S[i,j] + self.scorePos(seq1[i], seq2[j])
                s2 = self.S[i,j+1] + self.g
                s3 = self.S[i+1,j] + self.g
                self.S[i+1,j+1] = max(s1, s2, s3)
                self.T[i+1,j+1] = max3t(s1, s2, s3)
        return self.S[len(seq1),len(seq2)]   
    
    def needlemanWunschTies (self, seq1, seq2, ties = False):
        if (seq1.tipo != seq2.tipo): return None
        self.S = [[0]]
        self.T = [[0]]
        self.seq1 = seq1
        self.seq2 = seq2
        for j in range(1, len(seq2)+1):
            self.S[0].append(self.g * j)
            if ties: self.T[0].append([3])
            else: self.T[0].append(3)
        for i in range(1, len(seq1)+1):
            self.S.append([self.g * i])
            if ties: self.T.append([[2]])
            else: self.T.append([2])
        for i in range(0, len(seq1)):
            for j in range(len(seq2)):
                s1 = self.S[i][j] + self.scorePos (seq1[i], seq2[j])
                s2 = self.S[i][j+1] + self.g
                s3 = self.S[i+1][j] + self.g
                self.S[i+1].append(max(s1, s2, s3))
                if ties:
                    self.T[i+1].append(max3t_with_ties(s1, s2, s3))
                else:
                    self.T[i+1].append(max3t(s1, s2, s3))
        return self.S[len(seq1)][len(seq2)]
    
    def recoverAlignment (self):
        res = ["", ""]
        i = len(self.seq1)
        j = len(self.seq2)
        while i>0 or j>0:
            if self.T[i,j]==1:
                res[0] = self.seq1[i-1] + res[0]
                res[1] = self.seq2[j-1] + res[1]
                i -= 1
                j -= 1
            elif self.T[i,j] == 3:
                res[0] = "-" + res[0]
                res[1] = self.seq2[j-1] + res[1] 
                j -= 1
            else:
                res[0] = self.seq1[i-1] + res[0]
                res[1] = "-" + res[1]
                i -= 1
        return MyAlign(res, self.seq1.tipo)
    
    def recoverAlignment_with_ties (self):
        i = len(self.seq1)
        j = len(self.seq2)  
        alins = [["", "", i,j]]
        res = []
        while alins:
            al = alins.pop(0)
            i = al[2]
            j = al[3]
            if i==0 and j==0:
                res.append(al[:2])
            else:
                for t in self.T[i][j]:
                    p = []
                    if t==1:
                        p.append(self.seq1[i-1] + al[0])
                        p.append(self.seq2[j-1] + al[1])
                        p.append(i-1)
                        p.append(j-1)
                    elif t == 3:
                        p.append("-" + al[0])
                        p.append(self.seq2[j-1] + al[1])
                        p.append(i)
                        p.append(j-1)
                    else:
                        p.append(self.seq1[i-1] + al[0])
                        p.append("-" + al[1])
                        p.append(i-1)
                        p.append(j)
                    alins.append(p)
        return res
 
    def smithWaterman (self, seq1, seq2):
        if (seq1.tipo != seq2.tipo): return None
        self.S = MatrixNum(len(seq1)+1, len(seq2)+1)
        self.T = MatrixNum(len(seq1)+1, len(seq2)+1)
        self.seq1 = seq1
        self.seq2 = seq2
        maxscore = 0
        for j in range(1, len(seq2)+1):
            self.S[0, j] = 0
            self.T[0, j] = 0
        for i in range(1, len(seq1)+1):
            self.S[i, 0] = 0
            self.T[i, 0] = 0
        for i in range(0, len(seq1)):
            for j in range(len(seq2)):
                s1 = self.S[i,j] + self.scorePos(seq1[i], seq2[j]) 
                s2 = self.S[i,j+1] + self.g
                s3 = self.S[i+1,j] + self.g
                b = max(s1, s2, s3)
                if b <= 0:
                    self.S[i+1,j+1] = 0
                    self.T[i+1,j+1] = 0
                else:
                    self.S[i+1,j+1] = b
                    self.T[i+1,j+1]= max3t(s1, s2, s3)
                    if b > maxscore: 
                        maxscore = b
        return maxscore
    
    def smithWatermanTies (self, seq1, seq2, ties = False):
        if (seq1.tipo != seq2.tipo): return None
        self.S = [[0]]
        self.T = [[0]]
        self.seq1 = seq1
        self.seq2 = seq2
        maxscore = 0
        for j in range(1, len(seq2)+1):
            self.S[0].append(0)
            if ties: self.T[0].append([0])
            else: self.T[0].append(0)
        for i in range(1, len(seq1)+1):
            self.S.append([0])
            if ties: self.T.append([[0]])
            else: self.T.append([0])
        for i in range(0, len(seq1)):
            for j in range(len(seq2)):
                s1 = self.S[i][j] + self.scorePos(seq1[i], seq2[j]) 
                s2 = self.S[i][j+1] + self.g
                s3 = self.S[i+1][j] + self.g
                b = max(s1, s2, s3)
                if b <= 0:
                    self.S[i+1].append(0)
                    self.T[i+1].append(0)
                else:
                    self.S[i+1].append(b)
                    if ties:
                        self.T[i+1].append(max3t_with_ties(s1, s2, s3))
                    else:
                        self.T[i+1].append(max3t(s1, s2, s3))
                    if b > maxscore: 
                        maxscore = b
        return maxscore

    def recoverLocalAlignment (self):
        res = ["", ""]
        i, j = self.S.maxIndexes()
        while self.T[i,j]>0:
            if self.T[i,j]==1:
                res[0] = self.seq1[i-1] + res[0]
                res[1] = self.seq2[j-1] + res[1]
                i -= 1
                j -= 1
            elif self.T[i,j] == 3:
                res[0] = "-" + res[0];
                res[1] = self.seq2[j-1] + res[1]; 
                j -= 1
            elif self.T[i,j] == 2:
                res[0] = self.seq1[i-1] + res[0];
                res[1] = "-" + res[1]; 
                i -= 1
            else: break
        return MyAlign(res, self.seq1.tipo)


    def recoverAlignLocal_with_ties (self):
        maxval = self.S[0][0]
        maxtups = []
        for i in range(0,len(self.S)):
            for j in range(0, len(self.S[i])):
                if self.S[i][j] > maxval:
                    maxval = self.S[i][j]
                    maxtups = [(i,j)]
                elif self.S[i][j] == maxval:
                    maxtups.append((i,j))
        alins = []
        for (i,j) in maxtups:
            alins.append(["", "", i,j])        
        res = []
        while alins:
            al = alins.pop(0)
            i = al[2]
            j = al[3]   
            if (i==0 and j==0) or (0 in self.T[i][j]):
                res.append(al[:2])
            else:
                for t in self.T[i][j]:
                    p = []
                    if t==1:
                        p.append(self.seq1[i-1] + al[0])
                        p.append(self.seq2[j-1] + al[1])
                        p.append(i-1)
                        p.append(j-1)
                    elif t == 3:
                        p.append("-" + al[0])
                        p.append(self.seq2[j-1] + al[1])
                        p.append(i)
                        p.append(j-1)
                    else:
                        p.append(self.seq1[i-1] + al[0])
                        p.append("-" + al[1])
                        p.append(i-1)
                        p.append(j)
                    alins.append(p)
        return res


    def alignQuery (self, query, setSeqs):
        bestScore = -1
        bestAlin = None
        for seq in setSeqs:
            sc = self.smithWaterman(query, seq)
            if sc > bestScore:
                bestScore = sc
                bestAlin = self.recoverLocalAlignment()
        return bestAlin, bestScore
 
   
    def alignSeqWithSet(self, query, setSeqs):
        scores = []
        res = [[0,0]]*len(setSeqs)
        k = 0
        for s in setSeqs:
            scores[k] = self.smithWaterman(query, s)
            k = k + 1
        for i in range(len(setSeqs)):
            j = 0;
            while j < i and scores[ res[j][0] ] > scores[i]: 
                j = j + 1
            for k in range(i,j,-1):
                res[k][0] = res[k-1][0]
            res[j][0] = i
        for i in range(len(setSeqs)):
            res[i][1] = scores[ res[i][0] ]
        return res    
    
    def alignSets(self, c1, c2):
        res = []
        for i in range(len(c1)):
            maxscore = -1
            bestseq = None
            for j in range(len(c2)):
                score = self.smithWaterman(c1[i], c2[j])
                if score > maxscore:
                    maxscore = score
                    bestseq = c2[j]
            res.append(bestseq)
        return res
    
    def alignProtDNA(self, prot, dna):
        res = None
        orfs = dna.orfs();
        prots = []
        for i in range(len(orfs)): 
            prots.append( orfs[i].extraiTodasProts() )
        mx = 0
        for p in prots:
            score = self.smithWaterman(prot.seq, p.seq)
            if score > mx:
                mx = score
                res = self.getLocalAlignment()
        return res
   
    
def max3t (v1, v2, v3):
    if v1 > v2:
        if v1 > v3: return 1
        else: return 3
    else:
        if v2 > v3: return 2
        else: return 3

def max3t_with_ties(v1, v2, v3):
    if v1 > v2:
        if v1 > v3: 
            return [1]
        elif v1 == v3:
            return [1,3]
        else:
            return [3]
    elif v1 == v2:
        if v1 > v3: 
            return [1,2]
        elif v1 == v3:
            return [1,2,3]
        else:
            return [3]
    else:
        if v2 > v3: return [2]
        elif v2 == v3:
            return [2,3]
        else: return [3]

def printMat (mat):
    for i in range(0, len(mat)):
        print(mat[i])

def lcs (seq1, seq2):
    sm = SubstMatrix()
    sm.createFromMatchPars(1,0, seq1.alfabeto())
    aseq = AlignSeq(sm, 0)
    aseq.needlemanWunsch(seq1, seq2)
    alin = aseq.recoverAlignment ()
    sizeal = len(alin[0])
    lcs = ""
    for i in range(sizeal):
        if alin[0][i] != "-" and alin[0][i] == alin[1][i]:
            lcs += alin[0][i]
    return lcs        

    
def edit_distance(seq1, seq2):
    sm = SubstMatrix()
    sm.createFromMatchPars(0,-1,seq1.alfabeto())
    aseq = AlignSeq(sm, -1)
    sc = aseq.needlemanWunsch(seq1, seq2)
    return -sc
    
    
    
#### TESTS #####

def testScoresProt():
    s1 = MySeq("LGPS-GCASGIWTKSA", "protein")
    s2 = MySeq("TGPSGG--SRIWTKSG", "protein")
    alin = MyAlign([s1,s2], "protein")
    
    sm = SubstMatrix()
    sm.loadFromFile("blosum62.mat", "\t")
    
    aseq = AlignSeq(sm, -8) 
    score = aseq.scoreAlin(alin)
    
    print("Score", score)
 
# 29
    
def testScoresDNA():
    s1 = MySeq("ATGA-AGGT", "dna")
    s2 = MySeq("A-GAGAGGC", "dna")
    alin = MyAlign([s1,s2], "dna")
    
    sm = SubstMatrix()
    sm.createFromMatchPars(1,0, "ACGT")

    aseq = AlignSeq(sm, 0) 
    score = aseq.scoreAlin(alin)

    print("Score", score)

# 6

def test1():
    #seq1 = MySeq("NCRD","protein")
    #seq2 = MySeq("CRNC","protein")
    #seq1 = MySeq("DCNR","protein")
    #seq2 = MySeq("CNDR","protein")
    seq1 = MySeq("PHSWG","protein")
    seq2 = MySeq("HGWAG","protein")
    sm = SubstMatrix()
    sm.loadFromFile("blosum62.mat", "\t")
    #alin = AlignSeq(sm, -4)
    alin = AlignSeq(sm, -8)

    print("Alinhamento global")

    print("Score:", alin.needlemanWunsch(seq1, seq2))
    print("Matriz S:")
    alin.S.printmat()
    print("Matriz T:")
    alin.T.printmat()
    print("Alinhamento otimo:")
    print(alin.recoverAlignment())
    print()
    
    print("Alinhamento local")
    
    print("Score:", alin.smithWaterman(seq1, seq2))
    print("Matriz S:")
    alin.S.printmat()
    print("Matriz T:")
    alin.T.printmat()
    print("Alinhamento otimo:")
    print(alin.recoverLocalAlignment())

#Alinhamento global
#Score: 9
#Matriz S:
#[0, -8, -16, -24, -32, -40]
#[-8, -2, -10, -18, -25, -33]
#[-16, 0, -4, -12, -20, -27]
#[-24, -8, 0, -7, -11, -19]
#[-32, -16, -8, 11, 3, -5]
#[-40, -24, -10, 3, 11, 9]
#
#Matriz T:
#[0, 3, 3, 3, 3, 3]
#[2, 1, 3, 3, 1, 3]
#[2, 1, 1, 3, 3, 1]
#[2, 2, 1, 1, 1, 3]
#[2, 2, 2, 1, 3, 3]
#[2, 2, 1, 2, 1, 1]
#
#Alinhamento otimo:
#
#PHSW-G
#-HGWAG
#
#Alinhamento local
#Score: 19
#Matriz S:
#[0, 0, 0, 0, 0, 0]
#[0, 0, 0, 0, 0, 0]
#[0, 8, 0, 0, 0, 0]
#[0, 0, 8, 0, 1, 0]
#[0, 0, 0, 19, 11, 3]
#[0, 0, 6, 11, 19, 17]
#
#Matriz T:
#[0, 0, 0, 0, 0, 0]
#[0, 0, 0, 0, 0, 0]
#[0, 1, 0, 0, 0, 0]
#[0, 0, 1, 0, 1, 0]
#[0, 0, 0, 1, 3, 3]
#[0, 0, 1, 2, 1, 1]
#
#Alinhamento otimo:
#
#HSW
#HGW


# Result from test 1
# 7
# [0, -4, -8, -12, -16]
# [-4, -3, -4, -2, -6]
# [-8, 5, 1, -3, 7]
# [-12, 1, 10, 6, 3]
# [-16, -3, 6, 11, 7]
# 
# [0, 3, 3, 3, 3]
# [2, 1, 1, 1, 3]
# [2, 1, 3, 3, 1]
# [2, 2, 1, 3, 2]
# [2, 2, 2, 1, 3]
# 
# NCRD-
# -CRNC

def test2():
    seq1 = MySeq("ATGATATGATGATT")
    seq2 = MySeq("GATGAATAGATGTGT")
    sm = SubstMatrix()
    sm.createFromMatchPars(3, -1, "ACGT")
    alin = AlignSeq(sm, -3)
    print(alin.smithWaterman(seq1, seq2))
    printMat(alin.S)
    print(alin.recoverLocalAlignment())
    
# Results from test 2
# 25
# [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# [0, 0, 3, 0, 0, 3, 3, 0, 3, 0, 3, 0, 0, 0, 0, 0]
# [0, 0, 0, 6, 3, 0, 2, 6, 3, 2, 0, 6, 3, 3, 0, 3]
# [0, 3, 0, 3, 9, 6, 3, 3, 5, 6, 3, 3, 9, 6, 6, 3]
# [0, 0, 6, 3, 6, 12, 9, 6, 6, 4, 9, 6, 6, 8, 5, 5]
# [0, 0, 3, 9, 6, 9, 11, 12, 9, 6, 6, 12, 9, 9, 7, 8]
# [0, 0, 3, 6, 8, 9, 12, 10, 15, 12, 9, 9, 11, 8, 8, 6]
# [0, 0, 0, 6, 5, 7, 9, 15, 12, 14, 11, 12, 9, 14, 11, 11]
# [0, 3, 0, 3, 9, 6, 6, 12, 14, 15, 13, 10, 15, 12, 17, 14]
# [0, 0, 6, 3, 6, 12, 9, 9, 15, 13, 18, 15, 12, 14, 14, 16]
# [0, 0, 3, 9, 6, 9, 11, 12, 12, 14, 15, 21, 18, 15, 13, 17]
# [0, 3, 0, 6, 12, 9, 8, 10, 11, 15, 13, 18, 24, 21, 18, 15]
# [0, 0, 6, 3, 9, 15, 12, 9, 13, 12, 18, 15, 21, 23, 20, 17]
# [0, 0, 3, 9, 6, 12, 14, 15, 12, 12, 15, 21, 18, 24, 22, 23]
# [0, 0, 0, 6, 8, 9, 11, 17, 14, 11, 12, 18, 20, 21, 23, 25]
# 
# ATGATAT-GATGATT
# ATGA-ATAGATGTGT
   

def testTies():
    sm = SubstMatrix()
    sm.loadFromFile("blosum62.mat", "\t")
    alin = AlignSeq(sm, -1)
    
    seq1 = MySeq("GKYESVI")
    seq2 = MySeq("KYVSSWI")
    sc = alin.needlemanWunschTies(seq1, seq2, True)
    printMat(alin.S)
    printMat(alin.T)
    print("Melhor score do alinhamento otimo global:", sc)
    alins = alin.recoverAlignment_with_ties()
    for a in alins: print(a)
    
    alin.g = -3
    sc = alin.smithWaterman(seq1, seq2, True)
    print("Melhor score do alinhamento otimo local: " , sc)

    alinsL = alin.recoverAlignLocal_with_ties()
    for a in alinsL: print(a)


def testLCS():
    s1 = MySeq("ATTAGCT")
    s2 = MySeq("ATAAGCT")
    print(lcs(s1,s2))


def testED():
    s1 = MySeq("ATTAGCT")
    s2 = MySeq("ATTAAAGCT")
    print(edit_distance(s1,s2))

def testQ2():
    print("Introduza uma lista de sequencias de DNA terminada pela sequencia vazia")
    listaseqs = []
    seq = input("Introduza sequencia:").upper()
    while seq != "":
        listaseqs.append(MySeq(seq))
        seq = input("Introduza sequencia:")
    query = MySeq(input("Introduza query:"))
    print("Introduza parametros do alinhamento")
    g = int ( input("g = ") )
    match = int ( input("match = "))
    mismatch = int ( input("mismatch ="))
    sm = SubstMatrix()
    sm.createFromMatchPars(match, mismatch, "ACGT")
    alin = AlignSeq(sm, g)
    al, score = alin.alignQuery(query, listaseqs)
    print("Resultados:")
    print("Score: " , score)
    print("Alinhamento")
    print(al[0])
    print(al[1])

if __name__ == "__main__":   
#    print("Teste 1")
#    test1()
#    print("Teste 2")
#    test2()
#    print("Teste 3")
#    testLCS()
#    print("Teste 4")
#    testED()
#    print("Teste 5")
#    testTies()
    # print("Teste 6")
    # testQ2()

    testScoresProt()
    print()
    testScoresDNA()
    print()
    test1()