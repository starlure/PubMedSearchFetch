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


def Get_Hindex(first_name, last_name, MY_API_KEY, authtoken):

    # Convert author_name to author_id in scopus
    ## Initialize author search object and execute search
    resp = requests.get("http://api.elsevier.com/content/search/author?query=authfirst(" + first_name + ") and authlastname(" + last_name + ")",
                        headers={'Accept':'application/json',
                                 'X-ELS-APIKey': MY_API_KEY,
                                 'X-ELS-Authtoken': authtoken})
    auth_eid = resp.json().get("search-results").get("entry")[0].get("eid")
    
    # query author metrics using author_id in scopus database
    if auth_eid == None:
        return None
    else:
        resp = requests.get("http://api.elsevier.com/content/author?eid="+auth_eid+"&view=metrics",
                        headers={'Accept':'application/json',
                                 'X-ELS-APIKey': MY_API_KEY,
                                 'X-ELS-Authtoken': authtoken})
#    print (resp.json().get("author-retrieval-response")[0])
    
    #print (json.dumps(resp.json(),
    #                 sort_keys=True,
    #                 indent=4, separators=(',', ': ')))
        h_index = resp.json().get("author-retrieval-response")[0].get("h-index")
    #    print (h_index)
        return h_index
    
#if __name__ == "__main__":
#    Get_Hindex()