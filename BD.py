# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 22:30:38 2020

@author: mines
"""

import re 
import operator
from MySeq import MySeq

#dicionário com todas as nossas sequências de proteínas
class BD:
    def __init__(self):
        self.BD={}
        self.alfabeto = "ACDEFGHIKLMNPQRSTVWY"

    
    def __del__(self):
        pass
    
    def add_seq(self, id, seq):
        self.BD[str(id)] =(seq)
        
    def del_seq(self, id):
        del self.BD[str(id)]

    def displayseqs(self):
        return self.BD
    
    def get_seq_id(self, id):
        return self.BD[id]
    
    def valida(self,seq):
        alf=self.alfabeto
        pos=0
        while pos < len(seq):
            if seq[pos] not in alf:
                return False
            else:
                pos +=1
        return True
    
    # def createBD(self,dic):
    #     sequence=list(dic.values())
    #     id=list(dic.keys())
    #     for i in range(len(id)):
    #         if id[i] not in dic:
    #             self.BD[id[i]]=(sequence[i])
    #         else: 
    #             continue
    #     return self.BD
    
    def ficheiros(self,ficheiro):
        fasta={}
        with open(str(ficheiro)) as file_one:
            for line in file_one:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    active_sequence_name = line[1:]
                    if active_sequence_name not in fasta:
                        fasta[active_sequence_name] = []
                    continue
                sequence = line
                fasta[active_sequence_name].append(sequence)
    
        sequences=list(fasta.values())
        names=list(fasta.keys())
        for i in range(len(names)):
            if MySeq(''.join(sequences[i]),'protein'): # confirma se a sequencia é proteina
                sp=names[i].split('|')
                id=(sp[1])
                # a= sp[2].split(' ')
                # entry=a[0]
                # protein=re.search(' (.*?)OS',names[i]).group(1)
                # organism=re.search('OS=(.*?)GN',names[i])
                # if organism:
                #     organism=organism.group(1)
                # else:
                #     organism=re.search('OS=(.*?)PE',names[i]).group(1)
                # gene=re.search('GN=(.*?) ',names[i])
                # if gene:
                #     gene=gene.group(1)
                # else:
                #     gene=None
                
                # dic[id]=(entry,protein,organism,gene,sequences[i])
                self.add_seq(id,sequences[i])
                #if seq.valida():
                seq=MySeq(str(sequences[i]),'protein')
                dict[id] = {'seq': seq}
                dictionary = {'tamanho': seq.__len__(), 'maior proteína': seq.maiorProteina()}
                dict[id].update(dictionary)
        
            else: 
                return 'A sequência não é valida'
            #self.BD[id]=(sequences[i])
        return self.BD


    def write_seqs(self):
            
        try:
                # import re
                # file=open('sequencias.txt','a+') # abre o ficheiro sequencias.txt caso exista o ficheiro o apontador inicia no final caso não exista cria um novo em modo de leitura e escrita
                # seqs=file.readlines()  #lê as linhas do ficheiro
                # lista_seqs=[] # cria uma lista nova para colocar as seqs
                # for linha in seqs: # por cada linha no ficheiro 
                #     lista_seqs.append(re.split('[:\n]',linha)) #adicionamos à lista os seqs existentes (sem simbolos /n e :)
                #     print(lista_seqs) # NÃO CONSEGUE LER O QUE ESTÁ NO FICHEIRO ANTES!!
                # file.close() # fecha o ficheiro
                # ids=list(self.BD.keys()) # ids
                # sequencias=list(self.BD.values()) #sequencias
                # for i in range(len(ids)):
                #     if ids[i] not in lista_seqs:
                #         lista_seqs.append((ids[i],sequencias[i]))
                # tabela_seqs=[]
                # # print(lista_seqs)
                # new_file=open('sequencias.txt','a')
                # for i in range(len(lista_seqs)):
                #     # print(lista_seqs[i])
                #     idss=lista_seqs[i][0]
                #     # print(idss)
                #     seqss= ''.join(lista_seqs[i][1])
                #     print(idss+ ':' + seqss, file=new_file)
                # # print(res)
                # new_file.close()
            
        
            import re
            lista_ids=[]
            lista_seqss=[]
            file=open('sequencias.txt','r') # abre o ficheiro sequencias.txt caso exista o ficheiro o apontador inicia no final caso não exista cria um novo em modo de leitura e escrita
            seqs=file.readlines()  #lê as linhas do ficheiro
            lista_seqs=[] # cria uma lista nova para colocar as seqs
            for linha in seqs: # por cada linha no ficheiro 
                idss=re.findall('(.*?):',linha)
                seqss=re.findall(':(.*)',linha)
                lista_ids.append(idss[0])
                lista_seqss.append(seqss[0])
                lista_seqs.append(re.split('[:\n]',linha)) #adicionamos à lista os seqs existentes (sem simbolos /n e :)
            # print(lista_ids)
            # print(lista_seqss)    
            file.close()
           
            ids=list(self.BD.keys()) # ids
            sequencias=list(self.BD.values()) #sequencias
            print(ids)
            print(lista_ids)
            for i in range(len(ids)):
                if ids[i] in lista_ids: ##ISTO NÃO ESTÁ A FUNCIONAR
                    continue
                else:
                    lista_ids.append((ids[i]))
                    lista_seqss.append(''.join(sequencias[i]))
            print(lista_ids)
            print(lista_seqss)
            new_file=open('sequencias.txt','a')
            for i in range(len(lista_ids)):
                # print(lista_seqs[i])
                # idss=lista_seqs[i][0]
                # print(idss)
                # seqss= ''.join(lista_seqs[i][1])
                print(lista_ids[i]+ ':' + lista_seqss[i], file=new_file)
            # print(res)
            new_file.close()  
                    

                
        except:
            print('Problema no código de write_seqs')
    
    def see_fileseqs(self):
        ficheiros = open("sequencias.txt", "r") #abre o ficheiro em modo de leitura
        seqs_total = ficheiros.read()
        print('''      
                  
            **********  Sequências:  *************
                             
        ''')
        print(seqs_total) #lê todo o ficheiro
        ficheiros.close() #fecha o ficheiro
    


            
    
    
        