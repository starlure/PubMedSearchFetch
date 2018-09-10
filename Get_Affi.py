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


def Get_Affi(first_name, last_name, MY_API_KEY, authtoken):

    # Convert author_name to author_id in scopus
    ## Initialize author search object and execute search
    resp = requests.get("http://api.elsevier.com/content/search/author?query=authfirst(" + first_name + ") and authlastname(" + last_name + ")",
                        headers={'Accept':'application/json',
                                 'X-ELS-APIKey': MY_API_KEY,
                                 'X-ELS-Authtoken': authtoken})
    print (resp.json())
    if resp.json()==None or resp.json().get("search-results").get("entry")[0].get('affiliation-current')==None:
        return None
    else:
        affi = resp.json().get("search-results").get("entry")[0].get('affiliation-current').get('affiliation-name')
        return affi
    
#if __name__ == "__main__":
#    Get_Hindex()