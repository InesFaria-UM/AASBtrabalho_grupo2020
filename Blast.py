#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 21:06:41 2020

@author: fernandavieira
"""

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

class My_Blast:
    def __init__(self,file_fasta,file_xml,db="nr"):
        """Recebe como parâmetros a proteína em ficheiro fasta, um ficheiro de output xml e a base de dados a utilizar """
        self.__file_fasta = file_fasta     #("__" antes, fica privado)
        self.__file_xml = file_xml
        self.__db = db
        
    def blast(self,TRESH):
        """faz o blast"""
        try:
            blast_records = SeqIO.parse(self.__file_fasta, "fasta")  
            save_file = open(self.__file_xml, "w")      #abre um ficheiro para posteriormente guardar os resultados em xml
            for blast_record in blast_records:
                result_handle = NCBIWWW.qblast("blastp", self.__db, blast_record.format("fasta"))
                save_file.write(result_handle.read())
            save_file.close()
            blast_records.close()
                                                    
        except:
            print("Erro na configuração do Blast")