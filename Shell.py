#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 21:39:03 2020

@author: fernandavieira
"""

from cmd import *
import re
from Blast import My_Blast
import functionsRE
from MySeq import MySeq
from ficheiro import Ficheiro
from SubstMatrix import SubstMatrix
from AlignSeq import AlignSeq
from MatrixNum import MatrixNum
from MyBlast import MyBlast
from upgma import UPGMA
from MultipleAlign import MultipleAlign
from BD import BD
from INOUT import INOUT
from AnaliseNCBI import AnaliseNCBI
import sys
import pickle
from Bio.ExPASy import ScanProsite
from Blast import My_Blast

try:
    pickle_in = open('dict.pickle', 'rb')
    dict = pickle.load(pickle_in)
except:
    dict = {}
    pickle_out = open('dict.pickle', 'wb')
    pickle.dump(dict, pickle_out)
    
class Shell(Cmd):
    intro= """ To start this program type 'start' """
    prompt="OPTION > "
    
    def do_start(self,arg):
        while True:
            print(''' \n
    
             Algoritmos para Análise de Seq Biológicas

    1. Insert a new sequence mannually
    2. Insert and export sequences from file txt or fasta   
    3. Search sequence by ID (visualize and change properties)
    4. Search from database the protein that is more similar to input
    5. Calculate the frequency of sub-sequencies (size k) in a sequence or DB
    6. Search ocorrences of patterns (regular expression)
    7. Filogenetic Tree 
    8. Multiple Alignment
    9.Search for possible functional domains.
    
    #Optional:

    10. Blast - (insert: Ficheiro.fasta Ficheiro.xml Database)
    11. sair
    ''')
            x=int(input('\n Choose the option you want:'))
            if x == 1:
                self.add_sequence()
            elif x == 2:
                self.ficheiro()
            elif x == 3:
                self.annotations()
            elif x == 4:
                self.search_best_align()
            elif x == 5:
                self.freq_symbol()
            elif x == 6:
                self.findpattern()
            elif x == 7:
                self.tree()
            elif x == 8:
                self.multiplealign1()
            elif x == 9:
                self.prosite()
            elif x==10:
                self.blast()
            elif x==11:
                self.sair()
                # sys.exit()
            else:
                print('\n Opção não válida, escolha outra vez')
                self.do_start()

        
    def add_sequence(self):
        # try:
        id=input('ID: ')
        seq1=input("Sequence: ")
        seq= MySeq((seq1.split()[0]),'protein')
        print(seq.valida())
        if seq.valida():
            dict[id] = {'seq': seq}
            dictionary = {'tamanho': seq.__len__(), 'maior proteína': seq.maiorProteina()}
            dict[id].update(dictionary)
        
    def fasta(filename):
        fh = open(filename)
        name = fh.readline()[1:-1]
        sequence = ''
        for line in fh:
            sequence += sub('\n', '', line)     # sub - substititui white spaces (\s) por '' em cada linha...penso eu
        #print('Nome: ' + name)
        #print(name)
        #print('Sequência: ' + sequence)
        #print(sequence)
        fh.close()
        return sequence
    
    def text(filename):
        file = open(filename)
        linhas = file.readlines()
        seq = ''
        if linhas[0][0] != '>':
            for li in linhas:
                seq += sub('\n', '', li)
        else:
            for li in range(1, len(linhas)):
                seq += sub('\n', '', linhas[li])
        file.close()
        #print(seq)
        return seq
        
    
    def ficheiro(self):
        fich=input("Qual é o ficheiro? ")
        # ficheiros=Ficheiro(fich)
        # bd=BD()
        BDS.ficheiros(fich) #coloca as seq do ficheiro na BD (base de dados) devolvendo um dicionario
        dict_seq = BDS.displayseqs()
        print(BDS.write_seqs())
        BDS.see_fileseqs()
        # seq=(bd.get_seq_id('P27748')) # vai buscar a sequencia pelo id
        # print(''.join(seq)) #Devolve a seq sem parênteses nem '', como string

    def annotations(self):
        id=input("ID:")
        op1=int(input("Do you want do insert(1), edit(2) or visualize(3)? "))
        if op1 == 1:
            featname=input("type of annotation: ")
            feature=input("insert: ")
            if id in dict.keys():
                v=dict[id]
                print(v)
                dictionary = {featname: feature}
                dict[id].update(dictionary)
                print(dict)
        elif op1 == 2:
            feature = input('Feature that you want to change: ')
            new = input("new value: ")
            if dict[id][feature] is not None:
                dict[id][feature] = new
                print("\nFeature editada!!\n", dict[id])
        elif op1==3:
            print(dict[id]['seq'].seq)
        else: 
            return None
        
    def search_best_align(self):
        op = int(input("from seq(1) or file(2): "))
        # seq_alin=input("Sequence for alignment: ")
        blast = MyBlast()
        for iD in dict.keys():
            seq = (dict[iD]['seq'])
            blast.addSequenceDB(seq)
        if op == 1:
            seq_alin=input("Sequence for alignment: ")
            res = blast.bestAlignment(seq_alin)
            seqdb = blast.db[res[4]]
            for ID in dict.keys():
                if dict[ID]['seq'] == seqdb:
                    print("best match: ", dict[ID]['seq'])
                    
        elif op == 2:
            file = input('file: ')
            seq = functionsRE.text(file)
            res = blast.bestAlignment(seq)
            seqdb = blast.db[res[4]]
            for ID in dict.keys():
              if dict[ID]['seq'] == seqdb:
                  print("best match: ", dict[ID]['seq'])
                  
    def freq_symbol(self):
        op = int(input('Specific protein (id)(1) or in BD(2) ?'))
        k = int(input('choose size: '))
        if op == 1:
            dic={}
            id=str(input('sequence id: '))
            if id in dict.keys():
                seq1 = dict[id]["seq"]
                for i in range((int(dict[id]['tamanho']))-k+1):
                    aa=seq1[i:i+k]
                    if aa in dic:
                        dic[aa] += 1
                    else:
                        dic[aa] = 1
            return print(dic)
        
        if op==2:
            dic={}
            lista=[]
            id=[]
            for i in dict.keys():
                print(i)
                sequencias=dict[i]['seq'].seq
                print(sequencias)
                id.append(i)
                lista.append(sequencias)
                
            for i in range(len(lista)):
                for j in range((dict[id[i]]['tamanho'])-k+1):
                    aa=lista[i][j:j+k]
                    if aa in dic:
                       dic[aa] += 1
                    else:
                      dic[aa] = 1
            return print(dic)
                
    def findpattern(self):
        op=int(input('Specific protein (id)(1) or in BD(2)? '))
        expreg=str((input('Find pattern: ')))
        if op == 1:
            id=str(input('sequence id: '))
            seq1 = dict[id]["seq"].seq
            lista=re.findall(expreg,seq1)
            return print(lista)
    
        elif op == 2:
            lista1=[]
            id=[]
            padroes=[]
            for i in dict.keys():
                print(i)
                sequencias=dict[i]['seq'].seq
                print(sequencias)
                id.append(i)
                lista1.append(sequencias)
            for i in range(len(lista1)):
                # seq1 = dict[id[i]]["seq"].seq
                # print(seq1)
                lista2=re.findall(expreg,(dict[id[i]]["seq"].seq))
                padroes.append(lista2)
            return print(padroes)
            
    def tree(self):
        bd=input("choose ids for tree: ")
        bd_list=bd.split(' ')
        bd_seq=[]
        for iD in bd_list:
            bd_seq.append (dict[iD]['seq'])
        sm = SubstMatrix()
        sm.createFromMatchPars (1,0,"ACDEFGHIKLMNPQRSTVWY")
        alin = AlignSeq (sm,-1)
        up = UPGMA (bd_seq,alin)
        up.matdist.printmat()
        arv = up.run()
        arv.printtree()
        
    def multiplealign1(self):
        bd=input("choose ids for multiple align: ")
        bd_list=bd.split(' ')
        bd_seq=[]
        for iD in bd_list:
            bd_seq.append (dict[iD]['seq'])
        sm = SubstMatrix()
        sm.createFromMatchPars (1,0,"ACDEFGHIKLMNPQRSTVWY")
        alin= AlignSeq (sm,-1)
        ma = MultipleAlign(bd_seq, alin)
        alinm = ma.alignConsensus()
        print("score:", ma.scoreSP(alinm))
        print(alinm)
        
    def prosite(self):
        seq=input("select sequences id: ")
        seq_select = seq.split(' ')
        res = []
        for iD in seq_select:
            handle = ScanProsite.scan(seq = seq_select[iD]['seq'].seq)
            res.append(ScanProsite.read(handle))
        return res
    
    def blast(self):
        arg1= input('File Fasta:')
        arg2=input('File export: ')
        arg3=input('Database: ')
        Blast=My_Blast(arg1,arg2,arg3)
        TRESH=input('Qual é o valor do e-value Tresh: ')
        Blast.blast(TRESH)

    def sair(self):
        pickle_out = open('dict.pickle', 'wb')
        pickle.dump(dict, pickle_out)
        pickle_out.close()
        print(dict)
        sys.exit()
    
if __name__ == "__main__":
    BDS=BD()
    a=AnaliseNCBI()
    janela=None
    sh=Shell()
    sh.cmdloop()
  

    
    
    
    


    