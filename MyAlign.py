class MyAlign:

    def __init__(self, lseqs, tipo="protein"):
        self.listseqs = lseqs
        self.tipo = tipo
    
    def __len__(self):# number of columns
        return len(self.listseqs[0])
    
    def __getitem__(self, n):
        if type(n) is tuple and len(n) ==2: 
            i, j = n
            return self.listseqs[i][j]
        elif type(n) is int: return self.listseqs[n]
        return None
    
    def __str__(self):
        res = ""
        for seq in self.listseqs:
            res += "\n" + seq 
        return res
    
    def numSeqs(self):
        return len(self.listseqs)
    
    def consensus (self):
        cons = ""
        for i in range(len(self)):
            cont = {}
            for k in range(len(self.listseqs)):
                c = self.listseqs[k][i]
                if c in cont:
                    cont[c] = cont[c] + 1
                else: 
                    cont[c] = 1
            maximum = 0
            cmax = None
            for ke in cont.keys():
                if ke != "-" and cont[ke] > maximum: 
                    maximum = cont[ke]
                    cmax = ke
            cons = cons + cmax
        return cons
    
    def column (self, indice):
        res = []
        for k in range(len(self.listseqs)):
            res.append(self.listseqs[k][indice])
        return res
    
    def consensus_column(self,indice_col):
        listaCar=self.column(indice_col)
        res=True
        i=1
        while res and i<len(listaCar):
            if listaCar[i]!= listaCar[0]: 
                res=False
            else: 
                i=i+1
        return res
    
    def consensus_columns(self):
        res=[]
        for c in range(len(self)):
            if self.consensus_column(c):
                res.append(c)  #como se não for não fazemos nada, não pomos else
        return res
    
    def column_has_kgaps(self,indice_col, k):
        listaCar=self.column(indice_col)
        numgaps=0
        for c in len(listaCar):
            if c=="-":
                numgaps+=1
        if numgaps>=k:
            return True
        else:
            return False
        
    def columns_with_kgaps(self,k):
        res=[]
        for c in range(len(self)):
            if self.column_has_kgaps(c,k):
                res.append(c)
        return res
        
    def add_gap(self,seq,pos):
         self.listseqs[seq]= self.listseqs[seq][:pos]+"-"+self.listseqs[seq][pos:]
         
         #MyAlign seqs como lista de strings

    def remove_gap(self,seq,pos):
        if self.listseqs[seq][col]=="-":
            self.listseqs[seq]= self.listseqs[seq][:pos]+self.listseqs[seq][pos+1:] 

    def add_sequence_heuristic(self,newseq):  #saiu no testes!!!!!
        newrow=" "
        i=0 #i do alinhamento
        j=0 # j da nova sequencia
        while i< len(self) and j< len(newseq):  #é melhor porque não sabemos qual acaba primeiro
            carscol=self.column[i]
            if newseq[j] in carscol:
                newrow += newseq[j]
                j+=1
                
            else:
                newrow+="-"
            i+=1
        while j< len(newseq):
            newrow+= newseq[j]
            for k in range(len(self.listseqs)):
                self.listseqs[k]+="-"
            j+=1
        if i< len(self):
            for c in range(i,len(self)):
                newrow +="-"
        #while i< len(self): # chama self.__len__()
        #newrow+="-"
        #i+=1
        self.listseqs.append(newrow)


if __name__ == "__main__":
    alig = MyAlign(["ATGA-A","AA-AT-"], "dna")
    print(alig)
    print(len(alig))
    print(alig.column(2))
    print(alig[1,1])
    print(alig[0,2])

# Results
# ATGA-A
# AA-AT-
# 6
# ['G', '-']
# A
# G
