# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 11:25:17 2018
The function works, and needs to connect with "PubMed_SearchFetch.py"

@author: Liang Wang
"""

import requests
import json
from elsapy.elsclient import ElsClient
from elsapy.elssearch import ElsSearch


def Get_Abstract(pubmed_id, MY_API_KEY, authtoken):
    
    # Convert author_name to author_id in scopus
    ## Initialize author search object and execute search
#    pubmed_id = '1372328'
    resp = requests.get("http://api.elsevier.com/content/abstract/pubmed_id/" + pubmed_id,
                        headers={'Accept':'application/json',
                                 'X-ELS-APIKey': MY_API_KEY,
                                 'X-ELS-Authtoken': authtoken})
    print (resp.json())
    if resp.json()==None or resp.json().get("abstracts-retrieval-response").get("coredata")==None :
        return None
    else:
        abstract = resp.json().get("abstracts-retrieval-response").get("coredata").get("dc:description")
        return abstract
    
#if __name__ == "__main__":
#    Get_Hindex()