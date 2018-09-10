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
from Get_Hindex import Get_Hindex
from Get_Affi import Get_Affi
from ArticleRetrieval_Elsevier import Get_Abstract
import requests


class Article():
    
    def __init__(self, record, MY_API_KEY, authtoken):
        self.PMID = record.get('PMID',None)
        self.TI = record.get('TI',None) # Title
        self.FAU = record.get('FAU',None) # Full Author Name
        self.AB = record.get('AB',None) # Abstract
        self.JT = record.get('JT',None) # Journal Name
        self.DP = record.get('DP',None) # Publication Date
        self.AD = record.get('AD',None) # Affiliations
        self.GR = record.get('GR',None) # Grant
        self.MY_API_KEY = MY_API_KEY
        self.authtoken = authtoken
        
    def get_LastAuthorEmail(self):
        match = re.findall(r'[\w\.-]+@[\w\.-]+\.\w+',self.AD)
        if not match:
            return ''
        else:
            return match[-1]
    
    def get_FirstLastAuthor(self):
        if self.FAU == None:
            (FirstAuthor, LastAuthor, FirstAuthorAffi, LastAuthorAffi) = [None] * 4
        else:
            FirstAuthor, LastAuthor = self.FAU[0], self.FAU[-1]
            last_name, first_name = [x.strip() for x in FirstAuthor.split(',')]
            FirstAuthorAffi = Get_Affi(first_name, last_name, self.MY_API_KEY, self.authtoken)
            last_name, first_name = [x.strip() for x in LastAuthor.split(',')]
            LastAuthorAffi = Get_Affi(first_name, last_name, self.MY_API_KEY, self.authtoken)
        
        return (FirstAuthor, LastAuthor, FirstAuthorAffi, LastAuthorAffi)
    
    def get_LastAuthor_Hindex(self):
        if self.FAU == None:
            return None
        else:
            last_name, first_name = [x.strip() for x in self.FAU[-1].split(',')]
            Hindex = Get_Hindex(first_name, last_name, self.MY_API_KEY, self.authtoken)
            return Hindex
    
    def get_FirstLastAffi(self):
        author_len = len(self.FAU)
        if author_len == 0:
            FirstAuthorAffi, LastAuthorAffi = None, None
        elif author_len == 1:
            FirstAuthorAffi, LastAuthorAffi = self.AD, self.AD
        else:
            FirstAuthorAffi, LastAuthorAffi = self.AD.partition('.')[0], self.AD.partition('.')[-1]
            
#        FirstAuthorAffi = self.AD.partition('.')[0]
#        if LastAuthorEmail == '':
#            LastAuthorAffi = self.AD.split('.')[-2]
#        else:
#            tmp = self.AD.replace(LastAuthorEmail,'')
#            tmp_split = tmp.split('.')
#            if len(tmp_split) >= 3:
#                LastAuthorAffi = tmp_split[-3]
#            else:
#                LastAuthorAffi = tmp_split[-2]
        return (FirstAuthorAffi, LastAuthorAffi)
    
    def print_to_csv(self, writer, field_names):
        if self.AD == None:
            LastAuthorEmail = None
        else:
            LastAuthorEmail = self.get_LastAuthorEmail()
        (FirstAuthor,LastAuthor, FirstAuthorAffi, LastAuthorAffi) = self.get_FirstLastAuthor()
        LastAuthorHindex = self.get_LastAuthor_Hindex()
        Abstract = Get_Abstract(self.PMID, self.MY_API_KEY, self.authtoken)
#        print (LastAuthorEmail,FirstAuthor,FirstAuthorAffi,LastAuthor,LastAuthorAffi)
        values = [self.TI, self.FAU, FirstAuthor, FirstAuthorAffi, LastAuthor, LastAuthorAffi, 
                  LastAuthorEmail, LastAuthorHindex, Abstract, self.JT, self.DP, self.GR]
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
    search_handle = Entrez.esearch(db="pubmed", term = journal_name, mindate='2017/12/01', 
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
    MY_API_KEY = 'e526b5cb576ab5b4c1841c22e1939479'
    resp_auth = requests.get("http://api.elsevier.com/authenticate?platform=SCOPUS&choice=40772",
                        headers={'Accept':'application/json',
                                 'X-ELS-APIKey': MY_API_KEY})
    authtoken = resp_auth.json().get("authenticate-response").get("authtoken")
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
        field_names = ['Title','Authors','FirstAuthor','FirstAuthorAffiliation','LastAuthor','LastAuthorAffiliation',
                       'LastAuthorEmail','LastAuthorHindex','Abstract','Journal','PublicationDate','Grant']
        with open('outputs.csv','w',newline='') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=field_names)
            writer.writeheader()
            for record in records:
                article = Article(record, MY_API_KEY, authtoken)
    #            print (record)
                article.print_to_csv(writer, field_names)
    
if __name__ == "__main__":
    main()