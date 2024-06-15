import geopandas as gpd
import pandas as pd
import requests
import xml.etree.ElementTree as ET

def flatten_list(lst):
    return [item for sublist in lst for item in sublist]

def get_EGID(geo_pack_path, layer_):
    geo_pack = gpd.read_file(geo_pack_path, layer=layer_)
    if 'RegBL_EGID' in geo_pack.columns:
        return(geo_pack['RegBL_EGID'])
    else:
        egrid_full_list = []  # define egrid_full_list before if statement
        if 'EGRIS_EGRID_2' in geo_pack.columns:
            print("Column exists")
            list_EGRID = geo_pack['EGRIS_EGRID']
            list_EGRID2 = geo_pack['EGRIS_EGRID_2']
            egrid_full_list = list_EGRID.fillna(list_EGRID2)  # assign value to egrid_full_list
        else:
            egrid_full_list = geo_pack['EGRIS_EGRID']  # assign value to egrid_full_list

        egid_list = get_egid_list_from_egrids(egrid_full_list)

    return egid_list


#def get_egid_from_response(response_xml):
    #root = ET.fromstring(response_xml)
    #ns = {'eCH-0206': 'http://www.ech.ch/xmlns/eCH-0206/2'}
    #egid = root.find('.//eCH-0206:EGID', ns).text
    #return egid

def get_egid_list_from_response(response_xml):
    root = ET.fromstring(response_xml)
    ns = {'eCH-0206': 'http://www.ech.ch/xmlns/eCH-0206/2'}
    egid_list = []
    for egid in root.findall('.//eCH-0206:EGID', ns):
        egid_list.append(egid.text)
    return egid_list

def get_egid_list_from_egrids(egrids):
    url = 'https://madd.bfs.admin.ch/eCH-0206'
    headers = {'Content-Type': 'text/xml'}
    egid_list = []
    for egrid in egrids:
        madd_request_xml = f'''<?xml version="1.0" encoding="UTF-8"?>
<eCH-0206:maddRequest xmlns:eCH-0058="http://www.ech.ch/xmlns/eCH-0058/5" xmlns:eCH-0206="http://www.ech.ch/xmlns/eCH-0206/2" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.ech.ch/xmlns/eCH-0206/2 eCH-0206-2-0.xsd">
     <eCH-0206:requestHeader>
   	   	<eCH-0206:messageId>1710276222297</eCH-0206:messageId>
   	   	<eCH-0206:businessReferenceId>1710276222297</eCH-0206:businessReferenceId>
   	   	<eCH-0206:requestingApplication>
   	   	   	<eCH-0058:manufacturer>FSO</eCH-0058:manufacturer>
   	   	   	<eCH-0058:product>MADDAssist</eCH-0058:product>
   	   	   	<eCH-0058:productVersion>1.0.0</eCH-0058:productVersion>
   	   	</eCH-0206:requestingApplication>
   	   	<eCH-0206:requestDate>2024-03-12T21:43:42</eCH-0206:requestDate>
   	</eCH-0206:requestHeader>
   	<eCH-0206:requestContext>building</eCH-0206:requestContext>
   	<eCH-0206:requestQuery>
   	   	<eCH-0206:condition>
            	<eCH-0206:attributePath>/eCH-0206:maddResponse/eCH-0206:buildingList/eCH-0206:buildingItem/eCH-0206:realestateIdentificationList/eCH-0206:realestateIdentificationItem/eCH-0206:EGRID</eCH-0206:attributePath>
            	<eCH-0206:operator>equalTo</eCH-0206:operator>
            	<eCH-0206:attributeValue>{egrid}</eCH-0206:attributeValue>
        	</eCH-0206:condition>
   	</eCH-0206:requestQuery>
</eCH-0206:maddRequest>
'''
        response = requests.post(url, headers=headers, data=madd_request_xml)
        if response.status_code == 200:
            egid = get_egid_list_from_response(response.content)
            egid_list.append(egid)
            
        else:
            print(f'MADD Request failed for EGRID {egrid} with status code {response.status_code}')
    egid_list = flatten_list(egid_list)
    return egid_list

