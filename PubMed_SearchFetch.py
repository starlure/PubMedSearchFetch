# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 22:39:19 2018

@author: Liang Wang
"""

from Bio import Entrez
from Bio import Medline
import csv
import re
import sys
import argparse

class Article():
    
    def __init__(self, record):
        self.TI = record.get('TI',None)
        self.FAU = record.get('FAU',None)
        self.AB = record.get('AB',None)
        self.JT = record.get('JT',None)
        self.DP = record.get('DP',None)
        self.AD = record.get('AD',None)
        
    def get_LastAuthorEmail(self):
        match = re.findall(r'[\w\.-]+@[\w\.-]+\.\w+',self.AD)
        if not match:
            return ''
        else:
            return match[-1]
    
    def get_FirstLastAffi(self, LastAuthorEmail):
        FirstAuthorAffi = self.AD.partition('.')[0]
        if LastAuthorEmail == '':
            LastAuthorAffi = self.AD.split('.')[-2]
        else:
            tmp = self.AD.replace(LastAuthorEmail,'')
            tmp_split = tmp.split('.')
            if len(tmp_split) >= 3:
                LastAuthorAffi = tmp_split[-3]
            else:
                LastAuthorAffi = tmp_split[-2]
        return (FirstAuthorAffi, LastAuthorAffi)
    
    def print_to_csv(self, writer, field_names):
        if self.AD == None:
            (LastAuthorEmail,FirstAuthorAffi,LastAuthorAffi) = [None] * 3
        else:
            LastAuthorEmail = self.get_LastAuthorEmail()
            (FirstAuthorAffi,LastAuthorAffi) = self.get_FirstLastAffi(LastAuthorEmail)
#           print (LastAuthorEmail,FirstAuthorAffi,LastAuthorAffi)
        values = [self.TI, self.FAU, FirstAuthorAffi, LastAuthorAffi, LastAuthorEmail, 
                  self.AB, self.JT, self.DP]
#        print (dict(zip(field_names, values)))
        writer.writerow(dict(zip(field_names, values)))
    
def search_term(keyword):
    # if certain key terms needs to be searched in database
    # output article PMID list in "idlist"
    search_handle = Entrez.esearch(db="pubmed", term = keyword) # default retmax=20
    record = Entrez.read(search_handle)
#    print (record["Count"]) # number of available results
    search_handle.close()
    idlist = record["IdList"]
#    print (idlist)
    return idlist

def search_journal(journal_name="JAMA[ta]"):
    # if certain key terms needs to be searched in database
    # output article PMID list in "idlist"
    search_handle = Entrez.esearch(db="pubmed", term = journal_name, mindate='2017/01/01', 
                                   maxdate='2017/12/31', datetype='pdat', retmax=100000) # default retmax=20
    record = Entrez.read(search_handle)
#    print (record["Count"]) # number of available results
    search_handle.close()
    idlist = record["IdList"]
#    print (idlist)
    return idlist
    
def fetch_PMID(idlist):
    fetch_handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
    records = list(Medline.parse(fetch_handle))
    return records

def main():
    Entrez.email = "jasonwang.info@gmail.com"
    
#    parser = argparse.ArgumentParser(description='Get data of papers from Pubmed database')
#    parser.add_argument('-i', dest='inputfile', type=str, nargs=1,
#                    help='output data with given input PMID file')
#    parser.add_argument('-s', dest='keyword', type=str, nargs='+', default ='orchid',
#                    help='search papers of the keywords (default: output the first 20 papers)')
#    args = parser.parse_args()
##    print (sys.argv)
##    print (args)
#    if sys.argv[1] == '-s':
#        keyword = ' '.join(args.keyword)
#        idlist = search_term(keyword)
#    elif sys.argv[1] == '-i':
#        with open(args.inputfile[0]) as fobj:
#            idlist = [s.rstrip() for s in fobj.readlines()]
#    else:
#        print ('test a random paper')
#        idlist = ['29769726']
#    print (idlist)
        
    idlist = search_journal("JAMA[ta]")
    if not idlist:
        print ('no results for the keyword search, try another one.')
    else:
        records = fetch_PMID(idlist)
        field_names = ['Title','Authors','FirstAuthorAffiliation','LastAuthorAffiliation',
                       'LastAuthorEmail','Abstract','Journal','PublicationDate']
        with open('outputs.csv','w',newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=field_names)
            writer.writeheader()
            for record in records:
                article = Article(record)
    #            print (record)
                article.print_to_csv(writer, field_names)
    
if __name__ == "__main__":
    main()