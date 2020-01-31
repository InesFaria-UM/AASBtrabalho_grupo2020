# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 11:56:07 2019

@author: mines
"""

'''

Análise da sequência e das features presentes no NCBI

Deverá desenvolver scripts em BioPython que lhe permitam:
aceder ao NCBI e guardar os ficheiros correspondentes aos genes escolhidos, podendo
explorar possíveis variantes;
• verificar as anotações correspondentes aos genes de interesse;
• verificar e analisar a informação complementar fornecida pela lista de features e seus
qualifiers; pode usar os campos de referências externas para identificar identificadores
de outras bases de dados que permitam solidificar o conhecimento em relação a cada
gene.

'''
from Bio import Entrez
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation

class AnaliseNCBI:    
    def __init__(self):
        pass 
    
    def search(self,query):
        """Retira do NCBI a lista de ids por ordem de relevância"""
        Entrez.email='marinapt88@gmail.com'
        handle=Entrez.esearch(db='protein',sort='relevance',retmax='20',retmode='xml',term=query)
        results=Entrez.read(handle)
        return results
    
    def fetch_details(self,id_list):
        """Associa o id ao organismo correspondente"""
        ids=','.join(id_list)
        Entrez.email='marinapt88@gmail.com'
        handle=Entrez.esummary(db='protein',id=ids)
        record=Entrez.read(handle)
        return record
    
    def get_Chr(self,record,gene,pos):
        '''
         Retira ficheiro Genbank do gene id, pegando nas posições do cromossoma onde este se localiza 
        '''
        ChrAccVer=record['DocumentSummarySet']['DocumentSummary'][pos]['GenomicInfo'][0]['ChrAccVer']
        StartPos=record['DocumentSummarySet']['DocumentSummary'][pos]['GenomicInfo'][0]['ChrStart']
        EndPos=record['DocumentSummarySet']['DocumentSummary'][pos]['GenomicInfo'][0]['ChrStop']
        from Bio import Entrez
        Entrez.email = "marinapt88@gmail.com"  # Always tell NCBI who you are
        file_name='GB_'+gene+'_'+ str(ChrAccVer)+'_'+ '.txt'
        file=open(file_name,'w')
        handle = Entrez.efetch(db="nucleotide", id=str(ChrAccVer), rettype="gb", retmode="text",seq_start=StartPos, seq_stop=EndPos ) #Podemos colocar isto mais interativo e podermos escolher o formato final ex fasta ou genbank,etc
        a=handle.read()
        file.write(a)
        file.close()
        return file_name
        
    
    def get_annotations(self,file_name):
        records=SeqIO.parse(file_name,"gb")
        for seq_record in records:
            print(seq_record.id,seq_record.description[:100], "...")
            print("Sequence length:", len(seq_record))
            print(len(seq_record.features), "features")
            print("from: ", seq_record.annotations["source"])
            print('''
            
            *** Features ****
            
            ''')
            featcds=[]
            protein_id=[]
            for i in range(len(seq_record.features)):
                if seq_record.features[i].type=='CDS':
                    featcds.append(i)
            for k in featcds:
                print('\n')
                print(seq_record.features[k])
                print(seq_record.features[k].extract(seq_record.seq))
            for k in featcds:
                print(seq_record.features[k].location)
                protein_id.append(seq_record.features[k].qualifiers['protein_id'])
            proteinas=[]
            for i in protein_id:
                a=str(i).strip('[]')
                proteinas.append(a)
            print(proteinas)
            return proteinas #posso colocar a lista do __init__ e guardar lá depois
        
            # for k in featcds:
            #     print(seq_record.features[k].extract(seq_record.seq))
    
    
            
    def get_fasta(self,proteins,gene,protein_id): #lista com os protein_id de cada CDS
        Entrez.email ="marinapt88@gmail.com"
        file_name='fasta_'+gene+'_'+protein_id+'_'+ '.txt'
        file=open(file_name,'w')
        handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id=protein_id)
        # seq_record = SeqIO.read(handle, "fasta")
        a=handle.read()
        print(a)
        # print('handle:', handle)
        # print('seq_record: ', seq_record)
        # print (seq_record.id, " com ", len(seq_record.features), " features ")
        file.write(a)
        file.close() 
        handle.close() 
        return file_name              
        
  
    

    
    
    

# def seq_web(db,i):
#     '''
#      Retira ficheiro Genbank do gene id
#     '''
#     # Entrez.email='marinapt88@gmail.com'
#     file_name='GB_'+ gene+'_'+i+ '.txt'
#     file=open(file_name,'w')
    # handle=Entrez.efetch(db=db,rettype='gb',retmode='text', id=i)
    # # seq_record=SeqIO.read(handle,'gb')
    # file.write(handle)
    # file.close()
    # return handle
    # return seq_record
        
    # for seq_record in SeqIO.parse(file_name,"gb"):
    #     print(seq_record.id,seq_record.description[:100], "...")
    #     print("Sequence length:", len(seq_record))
    #     print(len(seq_record.features), "features")
    #     print("from: ", seq_record.annotations["source"])
    # handle.close()
    


#if __name__=="__main__":
#    gene=str(input("O que deseja procurar?"))
#    results=search(gene)
#    id_list=results['IdList']
#    record=fetch_details(id_list)
#    print(record)
#    for i, rec in enumerate(range(len(id_list))):
#         n_id = record['DocumentSummarySet']['DocumentSummary'][i].attributes['uid']
#         sc_name = record['DocumentSummarySet']['DocumentSummary'][i]['Organism']['ScientificName']
#         # ChrAccVer=record['DocumentSummarySet']['DocumentSummary'][0]['GenomicInfo'][0]['ChrAccVer']
#         '''Fazer com que a posição escolhida 1 seja o humano'''
#         # StartPosition=record['DocumentSummarySet']['DocumentSummary'][i]['GenomicInfo'[2]]
#         # EndPosition=record['DocumentSummarySet']['DocumentSummary'][i]['GenomicInfo'[3]]
#         print( '%d) %s - %s ' % (i+1, n_id, sc_name))
#    dic={1:0,2:1,3:2,4:3,5:4,6:5,7:6,8:7,9:8,10:9,11:10,12:11,13:12,14:13,15:14,16:15,17:16,18:17,19:18,20:19}
#    l=int(input('Em que posição se encontra da lista?'))
#    pos=dic[l]
#    file_name=get_Chr(record,gene,pos)
#    proteins=get_annotations(file_name)
#    get_fasta(proteins,gene)
#        


                
            
    
    
    
    
    
    
    
    
    
    
    # for id in id_list:
    #     record= fetch_details(id)
    #     print(id, record['DocumentSummarySet']['DocumentSummary'][0]['Organism']['ScientificName'])
        
  











# from Bio import SeqIO
# def get_record(file, format):
#     record=SeqIO.read(file,format)
#     print( "ID: ", record.id)
#     print("Name: ", record.name)
#     print("Sequence: ", record.seq)
#     print("Description: ", record.description) 
#     print("dbxrefs: ", record.dbxrefs)
#     print("Annotations: ", record.annotations)
#     print("Letter Annotations: ", record.letter_annotations)
#     print("Features: ", record.features) 
    
    


# if __name__=='__main__':
    
#     file=float(input("Insert the name of your file: "))
#     format= float(input("Format of file,'fasta', 'genbank': "))
    