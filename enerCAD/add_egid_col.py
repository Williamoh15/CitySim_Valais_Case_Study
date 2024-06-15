#import necessary packages

import geopandas as gpd
import pandas as pd
import requests
import xml.etree.ElementTree as ET
from shapely.geometry import Polygon, Point

#Breaks the sublists to have a clean list
def flatten_list(lst):
    return [item for sublist in lst for item in sublist]


#from the response of the madd request gets the egid and the coordinates of the RegBL (reference -> swiss grid)
def get_egid_list_from_response(response_xml):
    root = ET.fromstring(response_xml)
    ns = {'eCH-0206': 'http://www.ech.ch/xmlns/eCH-0206/2'}
    egid_list = []
    coord_list = []
    for building_item in root.findall('.//eCH-0206:buildingItem', ns):
        egid = building_item.find('.//eCH-0206:EGID', ns).text
        east = building_item.find('.//eCH-0206:east', ns).text
        north = building_item.find('.//eCH-0206:north', ns).text
        egid_list.append(egid)
        coord_list.append((east, north))
    return egid_list, coord_list


# Takes the geometry of the MO as argument and the corresponding egrid
def get_egid_list_from_egrids(egrids, geom):
    url = 'https://madd.bfs.admin.ch/eCH-0206'
    headers = {'Content-Type': 'text/xml'}
    egid_list = []
    coord_list = []
    for egrid in egrids:
        #Makes the Madd request for the egrid list
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
        #recuperates the response from the MADD request
        response = requests.post(url, headers=headers, data=madd_request_xml)
        if response.status_code == 200:
            #Get the list of each egid from the xml response file and the coordinates of each egid
            egids, coords = get_egid_list_from_response(response.content)
            if len(egids) > 1:
                # check which EGID has coordinates located in the geometrical figure given by the geopackage
                for egid, coord in zip(egids, coords):
                    if geom.contains(Point(coord)):
                        egid_list.append(egid)
                        coord_list.append(coord)
                        break
            else:
                # check that an egid corresponded to the egrid
                if len(egids):
                    egid_list.append(egids[0])
                    coord_list.append(coords[0])
                else:
                    print(f'No EGIDs found for EGRID {egrid}')
        else:
            print(f'MADD Request failed for EGRID {egrid} with status code {response.status_code}')
    return egid_list, coord_list

#Takes a geo_pack & layer to add an egid column for each building
def get_EGID(geo_pack_path, layer_):
    geo_pack = gpd.read_file(geo_pack_path, layer=layer_)
    if 'RegBL_EGID' in geo_pack.columns:
        return(geo_pack['RegBL_EGID'])
    else:
        egrid_full_list = []
        if 'EGRIS_EGRID_2' in geo_pack.columns:
            print("Column exists")
            list_EGRID = geo_pack['EGRIS_EGRID']
            list_EGRID2 = geo_pack['EGRIS_EGRID_2']
            egrid_full_list = list_EGRID.fillna(list_EGRID2)
        else:
            egrid_full_list = geo_pack['EGRIS_EGRID']

        egid_list = []
        for idx, row in geo_pack.iterrows():
            egrid = row['EGRIS_EGRID']
            geom = row['geometry']
            egids, _ = get_egid_list_from_egrids([egrid], geom)
            egid_list.append(egids[0] if egids else None)

        geo_pack['RegBL_EGID'] = egid_list
        return geo_pack

# Function that add the column to the geopackage   
def add_egid_to_geo (input_file_path,layer_name, output_file_path):

    geo_pack = get_EGID(input_file_path, layer_name)
    geo_pack.to_file(output_file_path, layer=layer_name, driver='GPKG')

    return None

def prepare_geo (input_file_path, layer):
    add_egid_to_geo (input_file_path, layer, input_file_path)
    return(None)



