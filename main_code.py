# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 17:16:28 2023

@author: Olivier Chavanne
"""

import geopandas as gpd
import pandas as pd
from shapely import box
import os
import matplotlib.pyplot as plt
import numpy as np
from shapely import wkt

# Local libraries
from enerCAD.building import generate_envelope
from enerCAD.building import generate_buildings
import enerCAD.xml as xml
import enerCAD.result as result
import enerCAD.network as network
import enerCAD.production as prod
import enerCAD.KPI as KPI
from enerCAD.add_egid_col import prepare_geo

# URL for RegBL API request
GEOADMIN_BASE_URL = "https://api.geo.admin.ch/rest/services/ech/MapServer/ch.bfs.gebaeude_wohnungs_register/"
    
##################################################
# 
#                  Functions
#
##################################################

def process_column_name(col):
    parts = col.split("(")
    if len(parts) > 1:
        # Join all parts except the first one, using '(' as separator
        return "(".join(parts[1:])
    else:
        return col



def calculate_azimuth(polygon):
    # Extract the coordinates of the polygon
    coordinates = np.array([list(point) for point in polygon.exterior.coords])

    # Calculate the surface normal of the polygon
    # (Assumes that the polygon is flat and that the vertices are in order)
    collinear = True
   

    v1 = coordinates[0] - coordinates[1]
    v2 = coordinates[1] - coordinates[2]

    normal = np.cross(v1, v2)

    azimuth = np.arctan2(normal[1], normal[0])

    angle_degrees = np.rad2deg(azimuth)

    azimuth_deg = 90-angle_degrees
    if azimuth_deg<0 :
        azimuth_deg = azimuth_deg + 360
        
    return azimuth_deg 

def create_xml_root(xml_file_to_copy, climate_file, horizon_file):
    '''
    Parameters                                                          
    ----------
    xml_file_to_copy : TYPE
        DESCRIPTION.
    climate_file : TYPE
        DESCRIPTION.
    horizon_file : TYPE
        DESCRIPTION.

    Returns
    -------
    root : TYPE
        DESCRIPTION.
    district : TYPE
        DESCRIPTION.
    '''
    
    # Write XML file for CitySim :
    print("Writing XML file...")    
    # Add Root 
    root = xml.add_root()
    # Add Simulation days
    xml.add_simulation_days(root)
    # Add Climate
    xml.add_climate(root, climate_file)
    # Add District
    district = xml.add_district(root)
    
    # Horizon
    # read in the tab-separated file as a dataframe
    horizon_df = pd.read_csv(horizon_file, sep='\t', header=None)
    # assign column names to the dataframe
    horizon_df.columns = ['phi', 'theta']
    # Add Far field obstructions
    xml.add_far_field_obstructions(district, horizon_df)
    
    # Add all the composites and profiles, taken from a source XML
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'Composite')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'OccupancyDayProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'OccupancyYearProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DeviceType')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'ActivityType')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DHWDayProfile')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DHWYearProfile')
    
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'Building')
    xml.add_child_from_xml_to_district(district, xml_file_to_copy, 'DistrictEnergyCenter')
    
    print("Xml source copied")
    
    return root, district 

def Module_1(gpkg_filepath, GEOADMIN_BASE_URL,
             directory_path, xml_name,
             xml_base_file, climate_file, horizon_file,
             create_geometry_3D=False, calculate_volume_3D=False,
             EGID_column='RegBL_EGID', renov_window = False, renov_wall = False):
    '''
    Parameters
    ----------
    gpkg_filepath : TYPE
        DESCRIPTION.
    GEOADMIN_BASE_URL : TYPE
        DESCRIPTION.
    directory_path : TYPE
        DESCRIPTION.
    xml_file_to_create : TYPE
        DESCRIPTION.
    xml_base_file : TYPE
        DESCRIPTION.
    climate_file : TYPE
        DESCRIPTION.
    horizon_file : TYPE
        DESCRIPTION.
    create_geometry_3D : TYPE, optional
        DESCRIPTION. The default is False.
    calculate_volume_3D : TYPE, optional
        DESCRIPTION. The default is False.
    EGID_column : TYPE, optional
        DESCRIPTION. The default is 'RegBL_EGID'.

    Returns
    -------
    None.
    '''
    
    ### Exctract geopackage ###
    
    print("Exctracting geopackage layers...")
    
    # MO Cadaster
    MO_all = gpd.read_file(gpkg_filepath, layer = "zone_tout")

    #check if the layer zone_tout contains the column regBL_EGID
    if EGID_column not in MO_all.columns:
        print('Adding EGID column to geopackage...')
        prepare_geo(gpkg_filepath, 'zone_tout')
        print('EGID column added to geopackage')
        MO_all = gpd.read_file(gpkg_filepath, layer = "zone_tout")
    
    MO_dhn = gpd.read_file(gpkg_filepath, layer = "zone_cad")

    #check if the layer zone_cad contains the column regBL_EGID
    if EGID_column not in MO_dhn.columns:
        print('Adding EGID column to geopackage...')
        prepare_geo(gpkg_filepath, 'zone_cad')
        print('EGID column added to geopackage')
        MO_dhn = gpd.read_file(gpkg_filepath, layer = "zone_cad")
        
    centrale = gpd.read_file(gpkg_filepath, layer = "centrale")

    EGID_column = 'RegBL_EGID'
    
    # Split Multipolygons into Polygons
    zone_all = MO_all.explode(index_parts=False)
    zone_dhn = MO_dhn.explode(index_parts=False)
    
    # List containing EGID of buildings to simulate
    EGID_list = MO_dhn[EGID_column].tolist()
    
    # Save EGID list of buildings connected to CAD
    df_EGID = pd.DataFrame(EGID_list)
    df_EGID.columns = ['EGID']
    EGID_path = os.path.join(directory_path, 'EGID.csv')     
    df_EGID.to_csv(EGID_path, index=False)
    print("EGID.csv created")
    
    # Swissbuildings3D
    print("Swissbuildings3D processing...")
    try:
        floor_data = gpd.read_file(gpkg_filepath, layer = "floor")
        roof_data = gpd.read_file(gpkg_filepath, layer = "roof")
        wall_data = gpd.read_file(gpkg_filepath, layer = "wall")
        
        # Filter on the zone with 10m buffer around surrounding square box 
        zone_bounds = MO_all.geometry.buffer(10).values.total_bounds
        zone_box = box(zone_bounds[0], zone_bounds[1], zone_bounds[2], zone_bounds[3])
        
        # Cut swissbuildings3D to zone of concern
        floor_data_intersection = floor_data[floor_data.geometry.intersects(zone_box)]
        roof_data_intersection = roof_data[roof_data.geometry.intersects(zone_box)]
        wall_data_intersection = wall_data[wall_data.geometry.intersects(zone_box)]
    
        # Split Multipolygons into Polygons
        zone_floor = floor_data_intersection.explode(index_parts=True).reset_index()
        zone_roof = roof_data_intersection.explode(index_parts=True).reset_index()
        zone_wall = wall_data_intersection.explode(index_parts=True).reset_index()


        print('Swissbuildings3D cut to zone of interest \n')
    
    except: print('Error : Swissbuildings3D not provided')

    ### Envelope processing ###
    
    try:
        # Get z coordinates of 1st vertex from 1st surface of 1st building's floor polygon as altitude by default for MO footprints
        altitude_default = zone_floor.loc[0].geometry.exterior.coords[0][2]
    except:
        altitude_default = 0
    
    # Create DataFrames containing all necessary information for each building
    print("Creating Buildings GeoDataFrame...")
    footprints, buildings = generate_buildings(zone_all, EGID_list, GEOADMIN_BASE_URL, altitude_default,
                                               create_geometry_3D, calculate_volume_3D, zone_floor, zone_roof, zone_wall,renov_window = renov_window, renov_wall = renov_wall)
    print("Buildings GeoDataFrame created \n") 
    
    # Generate the envelope surfaces
    print("Generating Buildings envelope...")
    envelope, buildings_volume_3D, center_coordinates = generate_envelope(footprints, buildings, calculate_volume_3D)
    print("Envelope created \n")
    
    # Merge "volume_3D" and "n_occupants" to main buildings geodataframe according to 'bid'
    merged_buildings = buildings.merge(buildings_volume_3D, left_on='bid', right_on='bid', how='left')    
    if not merged_buildings.empty:
        columns_to_add = ['volume_3D', 'n_occupants']
        for column in columns_to_add:
            buildings[column] = merged_buildings[column]
        print("Buildings 3D volume calculated and merged \n")
    
    ### Buildings XML processing ###
        
    root, district = create_xml_root(xml_base_file, climate_file, horizon_file)
    
    print("Adding buildings...")


    envelope["azimuth"] = envelope['geometry'].apply(calculate_azimuth)
    envelope['is_north_facing'] = np.where((envelope['azimuth'] > 45) & (envelope['azimuth'] < 315), True, False)
    
    # Add the buildings
    xml.add_all_buildings(district, buildings, envelope, center_coordinates)
    # write envelope to csv
    envelope_path = os.path.join(directory_path, 'Envelope.csv')
    envelope.to_csv(envelope_path, index=False)
        
    # Write XML file
    xml_path = os.path.join(directory_path, xml_name+".xml")     
    xml.write_xml_file(root, xml_path)
    print(f"{xml_name}.xml file created \n")

    # create geopackage in output folder 
    zone_all.to_file(output_gpkg_path, layer='zone_all', driver='GPKG')
    zone_dhn.to_file(output_gpkg_path, layer='zone_cad', driver='GPKG')
    centrale.to_file(output_gpkg_path, layer='centrale', driver='GPKG')

    
    return buildings, zone_dhn, centrale
   
def simulate_citysim(directory_path, xml_file, citysim_filepath):
    '''
    Parameters
    ----------
    xml_file : TYPE
        DESCRIPTION.

    Returns
    -------
    None.
    '''
    
    import subprocess
    import time
    start = time.time()
    print('Process started')
    print(f'Simulation of {xml_file}.xml...')

    #run CitySim.exe with xml file
    xml_path = os.path.join(directory_path, xml_file+".xml")
    result = subprocess.run([citysim_filepath, '-q', f"{xml_path}"])
    
    end = time.time()
    duration = end - start
    m, s = divmod(duration, 60)
    print('Simulation ended. Time :', "%.0f" %m,'min', "%.0f" %s,'s \n')



#------------------Part 2 iterating----------------------------------------------------------

def Module_2(directory_path, xml_name, zone_dhn, centrale,
             xml_DHN, climate_file, horizon_file,
             scenarios_list):   

    # Open the 1st iteration results file (without DHN)
    results_filepath = os.path.join(directory_path, xml_name+"_TH.out")
    results = pd.read_csv(results_filepath, delimiter="\t")
    results = results.set_index("#timeStep")
    print(f'{xml_name}_TH.out opened')
    
    # Get P_max for every building in the network
    power_EGID = result.get_Pmax_per_EGID(results)
    power_EGID_path = os.path.join(directory_path, 'Power_EGID.csv')
    power_EGID.to_csv(power_EGID_path, index=False)
        
    # Create and size network
    graph, lines_gdf, nodes_gdf, points, pipes, substations = network.get_trace(zone_dhn, centrale, power_EGID)
    
    # Save to csv file
    points_path = os.path.join(directory_path, 'Points.csv')
    pipes_path = os.path.join(directory_path, 'Pipes.csv')
    substations_path = os.path.join(directory_path, 'Substations.csv')
    points.to_csv(points_path, index=False)
    pipes.to_csv(pipes_path, index=False)
    substations.to_csv(substations_path, index=False)
    print("csv files saved")
    
    # Get Load duration curve of the network
    load_curve = result.get_thermal_load_curve(results)
    load_curve_path = os.path.join(directory_path, 'Load_curve.csv')
    load_curve.to_csv(load_curve_path, index=False)
     
    # Compute production scenarios
    scenarios = prod.get_scenarios(load_curve)
    
    # Calculate production water storage
    volume_storage, capacity_storage = prod.get_storage(load_curve)
    
    for i in range(len(scenarios)):
        scenario = scenarios[i]
        sc_id = scenario[0]
        if sc_id in scenarios_list:

            ### Scenarios XML processing ###
            
            xml_to_copy_path = os.path.join(directory_path, xml_name+'.xml' )
            root, district = create_xml_root(xml_to_copy_path, climate_file, horizon_file)
               
            # Add District heating center
            district_heating_center = xml.add_district_heating_center(district)
            
            # Add Network (nodes and pipes)
            xml.add_network(district_heating_center, points=points.copy(), pipes=pipes.copy())

            # Change Boilers into Substations
            xml.change_boiler_to_substation(district, substations=substations.copy(), points=points.copy())
            
            # Add Thermal station
            ts_node_id = points.loc[points['Type']=='start heating station']['npid'].iloc[0]
            network_p_max = pipes['power_line[W]'].max()
            production_stages = [scenario[1],scenario[2]]
            
            # Efficiency data
            technology_parameters = pd.read_csv("KPI.csv", delimiter=",")
            eff_columns = ['T_eff','efficiency']
            efficiency_parameters = technology_parameters[eff_columns]
            
            xml.add_thermal_station(district_heating_center, ts_node_id, p_max=network_p_max, 
                                    c_storage=capacity_storage, efficiencies=efficiency_parameters, stages=production_stages)

            # Write XML file
            scenario_path = os.path.join(directory_path,f"Scenario_{sc_id}")
            os.makedirs(scenario_path, exist_ok=True)
            xml_to_create_path = os.path.join(scenario_path, xml_DHN+f"_sc_{sc_id}"+".xml")
            xml.write_xml_file(root, xml_to_create_path)
            print(f'{xml_DHN}_sc_{sc_id}.xml file created \n')

    return graph, lines_gdf, nodes_gdf, results, scenarios, volume_storage

#--------------------- Part 3

def Module_results_network(scenario_path, sc_id, xml_DHN, zone_dhn, centrale, graph, lines_gdf, nodes_gdf):
    
    # Open the 2nd iteration results file (with DHN)
    results_filepath = os.path.join(scenario_path, xml_DHN+f"_sc_{sc_id}"+"_TH.out")
    results = pd.read_csv(results_filepath, delimiter="\t")
    results_final = results.set_index("#timeStep")
        
    # Get Data for every node and pipe in the network
    Pipes_mass_flow, Nodes_supply_temp, index_max = result.get_network_data_max(results_final)
    
    # Save to csv file
    Pipes_mass_flow_filepath = os.path.join(scenario_path, f'Pipes_mass_flow_sc_{sc_id}.csv')
    Nodes_supply_temp_filepath = os.path.join(scenario_path, f'Nodes_supply_temp_sc_{sc_id}.csv')
    Pipes_mass_flow.to_csv(Pipes_mass_flow_filepath, index=False) 
    Nodes_supply_temp.to_csv(Nodes_supply_temp_filepath, index=False)

    return results_final, Pipes_mass_flow, Nodes_supply_temp, index_max

#--------------------- KPI calculation

def Module_KPI(results_production, volume_storage, 
               scenarios, sc_id, scenario_path, do_plot):

    # Get production consumption in kWh
    pump_cons, fuel_cons, elec_cons, th_prod, df_th_prod, df_elec = result.get_energy_data(results_production, sc_id, scenarios)

    # Save thermal results to csv file
    th_production_results_filepath = os.path.join(scenario_path, f'Thermal_sc_{sc_id}.csv')
    df_th_prod.to_csv(th_production_results_filepath, index=True) 

    if do_plot == True:
    # Plot energy production data
        print('Energy production plot...')
        result.plot_energy_data(results_production, sc_id, scenarios, scenario_path)

    # Calculate KPI (key performance indicators)
    technology_parameters = pd.read_csv("KPI.csv", delimiter=",")

    df_KPI = KPI.calculate_KPI(sc_id, scenarios, volume_storage, technology_parameters,
                                                 pump_cons, fuel_cons, elec_cons, th_prod)
    print('KPI calculated')
    
    # Save electrical results to csv file
    electricity_results_filepath = os.path.join(scenario_path, f'Electricity_sc_{sc_id}.csv')
    df_elec.to_csv(electricity_results_filepath, index=True) 

    # Save KPI results to csv file
    KPI_results_filepath = os.path.join(scenario_path, f'KPI_results_sc_{sc_id}.csv')
    df_KPI.to_csv(KPI_results_filepath, index=False) 
    
    return df_KPI


##################################################
# 
#         Information to provide
#
##################################################

# Geopackage filepath
gpkg_filepath = r"C:\Users\Administrateur\Documents\Semester_Project\CAD-O-main\essai_new_mo\Martigny_clean.gpkg"       #TODO

# Create geometry with swissbuildings3D
create_geometry_3D = True                             #TODO

# Calculate volume from swissbuildings3D
calculate_volume_3D = True                               #TODO

# CitySim.exe filepath
citysim_filepath = r"C:\Program Files (x86)\CitySimSolver\CitySim.exe" #TODO

# XML name to export
directory_path = "Martigny_current"                                   #TODO

os.makedirs(directory_path, exist_ok=True)
                                      
xml_name = directory_path                                       
xml_DHN = "DHN_"+xml_name

# XML source files
xml_base_file = "xml_base.xml"                                
climate_file = r"C:\Users\Administrateur\Documents\Semester_Project\CAD-O-main\CH_Martigny_Verbier\CH_Martigny_CONTEMP.cli"                    #TODO
horizon_file = r"CH_Martigny_Verbier\Martigny.hor"                                 #TODO


#renovation_scenario
renov_window = False                                        #TODO
renov_wall = False                                        #TODO


# Scenarios to simulate
scenarios_list = [1,2,3,4,5,6]                                  #TODO

do_plot = True

# define output gpkg file
output_gpkg_path = os.path.join(directory_path, directory_path + ".gpkg")

def main(): 
    
    # Generate individual buildings XML
    print('***Module 1*** \n')
    buildings, zone_dhn, centrale = Module_1(gpkg_filepath, GEOADMIN_BASE_URL, 
                                             directory_path, xml_name,
                                             xml_base_file, climate_file, horizon_file,
                                             create_geometry_3D, calculate_volume_3D,
                                             EGID_column='RegBL_EGID', renov_window = renov_window, renov_wall = renov_wall)
    
    buildings.to_csv(directory_path  + "buildings.csv", index=False)


  

 
    # 1st CitySim simulation
    simulate_citysim(directory_path, xml_name, citysim_filepath)

    print ("writing back in geopackage")

    # Set the price and coefficient of performance (COP) of HP
    price = 0.23
    cop = 3

    # Rename and convert columns for merging
    buildings_copy = buildings.rename(columns={'egid': 'RegBL_EGID'})
    buildings_copy['RegBL_EGID'] = buildings_copy['RegBL_EGID'].astype(int)
    zone_dhn['RegBL_EGID'] = zone_dhn['RegBL_EGID'].astype(int)

    # Merge the buildings data with the zone data
    gdf_building_merge = zone_dhn.merge(buildings_copy[["wall_type", "habitable_surface", "RegBL_EGID"]], on='RegBL_EGID', how='inner')
    

    # Define the column names for the yearly results
    col_names = ['RegBL_EGID', 'heatingNeeds', 'coolingNeeds', 'solarPVProduction', 'solarThermalProduction']
    yearly_result_building = directory_path + '_YearlyResultsPerBuilding.out'
    output_TH = directory_path + '_TH.out'

    # Read and process hourly PV production data
    hourly_PV_prod = pd.read_csv(directory_path + '\\' + output_TH, sep='\t', header=0)
    hourly_PV_prod.columns = [process_column_name(col) for col in hourly_PV_prod.columns]
    delimiter = '_'

    # Filter columns for heating and solar PV production
    filtered_columns = hourly_PV_prod.filter(regex=r'SolarPVProduction\(Wh\)$|Heating\(Wh\)$')
    filtered_columns.columns = filtered_columns.columns.astype(str)
    filtered_columns.columns = filtered_columns.columns.str.replace('(', '').str.replace(')', '')
    filtered_columns.columns = [col[0:-2] for col in filtered_columns.columns]
    filtered_columns.columns = [delimiter.join([col.split(':')[0], col.split(':')[-1], "Wh"]) for col in filtered_columns.columns]

    # Extract unique building IDs
    column_names = [col.split('_')[0] for col in filtered_columns.columns]
    column_names = list(dict.fromkeys(column_names))

    # Insert additional columns for energy calculations
    for name in column_names:
        # Calculate the percentage of heating needs covered by PV
        filtered_columns.insert(filtered_columns.columns.get_loc(f'{name}_SolarPVProduction_Wh') + 1, f'{name}_Need_cov_pv(%)', 
                                ((100 - ((filtered_columns[f'{name}_Heating_Wh'] / cop - filtered_columns[f'{name}_SolarPVProduction_Wh']) * 100 / filtered_columns[f'{name}_Heating_Wh']) / cop).clip(upper=100)).fillna(0))
        # Calculate the electricity import needed
        filtered_columns.insert(filtered_columns.columns.get_loc(f'{name}_SolarPVProduction_Wh') + 2, f'{name}_Elec_import_pv_Wh', 
                                (filtered_columns[f'{name}_Heating_Wh'] / cop - filtered_columns[f'{name}_SolarPVProduction_Wh']))
        # Calculate the cost of heating without PV
        filtered_columns.insert(filtered_columns.columns.get_loc(f'{name}_SolarPVProduction_Wh') + 3, f'{name}_Cost_heat_no_pv_CHF', 
                                (filtered_columns[f'{name}_Heating_Wh'] / (1000 * cop)) * price)
        # Calculate the cost of electricity import
        filtered_columns.insert(filtered_columns.columns.get_loc(f'{name}_SolarPVProduction_Wh') + 4, f'{name}_Cost_elec_imp_CHF',
                                np.where(filtered_columns[f'{name}_Need_cov_pv(%)'] == 100,
                                        ((filtered_columns[f'{name}_Heating_Wh'] / cop) - filtered_columns[f'{name}_SolarPVProduction_Wh']) * 0.1690 / 1000,
                                        filtered_columns[f'{name}_Elec_import_pv_Wh'] * 0.2960 / 1000))

    # Sum the values for cost calculation
    sum_values = filtered_columns.sum()
    cost_w_PV = sum_values[sum_values.index.str.endswith('imp_CHF')]
    cost_w_PV.index = cost_w_PV.index.str.replace('_Cost_elec_imp_CHF', '')
    cost_w_PV = pd.DataFrame(cost_w_PV.reset_index())
    cost_w_PV.columns = ['RegBL_EGID', 'Cost_w_PV_CHF']
    cost_w_PV['RegBL_EGID'] = cost_w_PV['RegBL_EGID'].astype(int)

    # Read and process area data for PV panels
    output_Area = directory_path + '_Area.out'
    df_area = pd.read_csv(directory_path + '\\' + output_Area, sep='\t')
    df_PV_area = df_area[df_area['type'] == 'roofPVArea(m2)']
    df_PV_area = df_PV_area.melt(id_vars=['type'], var_name='RegBL_EGID', value_name='PV_panels_area_m2')

    # Remove the 'type' column as it's not needed anymore in the vertical format
    df_PV_area.drop(columns=['type'], inplace=True)
    df_PV_area['installation_cost_CHF'] = df_PV_area['PV_panels_area_m2'] * 300
    df_PV_area['RegBL_EGID'] = [process_column_name(col).replace(')', '') for col in df_PV_area['RegBL_EGID']]
    df_PV_area = df_PV_area[df_PV_area['RegBL_EGID'].str.isnumeric()]
    df_PV_area['RegBL_EGID'] = df_PV_area['RegBL_EGID'].astype(int)

    # Read and process yearly result data
    df_year_result = pd.read_csv(directory_path + '\\' + yearly_result_building, sep='\t', skiprows=1, names=col_names)
    df_year_result['RegBL_EGID'] = df_year_result['RegBL_EGID'].str.extract(r'\((\d+)\)', expand=False)
    df_year_result['heatingNeeds'] = pd.to_numeric(df_year_result['heatingNeeds'], errors='coerce')
    df_year_result.rename(columns={'heatingNeeds': 'heatingNeeds_kWh'}, inplace=True)
    df_year_result['coolingNeeds'] = pd.to_numeric(df_year_result['coolingNeeds'], errors='coerce')
    df_year_result.rename(columns={'coolingNeeds': 'coolingNeeds_kWh'}, inplace=True)
    df_year_result['solarPVProduction'] = pd.to_numeric(df_year_result['solarPVProduction'], errors='coerce')
    df_year_result.rename(columns={'solarPVProduction': 'solarPVProduction_kWh'}, inplace=True)
    df_year_result['RegBL_EGID'] = df_year_result['RegBL_EGID'].astype(int)

    # Merge yearly results with building data
    gdf_merge_2 = gdf_building_merge.merge(df_year_result, how='inner', on='RegBL_EGID')
    gdf_merge_2['heating_kWh_per_m2'] = gdf_merge_2['heatingNeeds_kWh'] / gdf_merge_2['habitable_surface']
    gdf_merge_2['habitable_surface'] = round(gdf_merge_2['habitable_surface'])
    gdf_merge_2.rename(columns={'habitable_surface': 'habitable_surface_m2'}, inplace=True)

    # Merge with cost and area data
    gdf_merge_2 = gdf_merge_2.merge(cost_w_PV, how='inner')
    gdf_merge_2 = gdf_merge_2.merge(df_PV_area, how='inner')

    # Calculate additional columns for cost and savings
    gdf_merge_2['heating_cost_no_PV_CHF'] = (gdf_merge_2['heatingNeeds_kWh'] / cop) * price
    gdf_merge_2['Energy_Savings_CHF'] = gdf_merge_2['heating_cost_no_PV_CHF'] - gdf_merge_2['Cost_w_PV_CHF']
    gdf_merge_2['PV_payback_time_years'] = gdf_merge_2['installation_cost_CHF'] / gdf_merge_2['Energy_Savings_CHF']

    # Write the final GeoDataFrame to a GeoPackage file
    gdf_merge_2.to_file(output_gpkg_path, layer='zone_cad_added_infos', driver='GPKG')



    #I commented out the module 2 as I didn't need it


    
    # Generate DHN XML for each scenario
    #print('***Module 2*** \n')
    #graph, lines_gdf, nodes_gdf, results, scenarios, volume_storage = Module_2(directory_path, xml_name, 
                                                                       # zone_dhn, centrale, 
                                                                       # xml_DHN, climate_file, horizon_file,
                                                                       # scenarios_list)
    # Intermediate results
    #network.plot_network_sizing(directory_path, graph, zone_dhn, lines_gdf, centrale)     
    #result.plot_load_curve(results, directory_path)
    
    #KPI_result_list = []
        
    # DHN simulation for each scenario
    #for i in range(len(scenarios_list)):
        #sc_id = scenarios_list[i]
        #print(f'***Scenario {sc_id}*** \n')
        
        # CitySim simulation        
        #scenario_path = os.path.join(directory_path,f"Scenario_{sc_id}")
       # simulate_citysim(scenario_path, f'{xml_DHN}_sc_{sc_id}', citysim_filepath)
                
        # Analysing final results  
        #print('Results processing...')
        #results_prod_path = os.path.join(scenario_path,  f'{xml_DHN}_sc_{sc_id}_TH.out')
       # results_production = pd.read_csv(results_prod_path, delimiter="\t")
            
        # KPI calculation
        #df_KPI = Module_KPI(results_production, volume_storage, 
                            #scenarios, sc_id, scenario_path, do_plot)
        
        #KPI_result_list.append(df_KPI)

        # Network final plots
       # results_final, Pipes_mass_flow, Nodes_supply_temp, index_max = Module_results_network(scenario_path, sc_id, xml_DHN, 
                                                                               # zone_dhn, centrale, 
                                                                               # graph, lines_gdf, nodes_gdf)   
       # if do_plot == True:
           # network.plot_network_data(scenario_path, sc_id, graph, zone_dhn, 
                                     # lines_gdf, nodes_gdf, centrale, index_max)
        
        #plt.close("all")
        #print(f"Scenario {sc_id} processed \n")
    
  #  if len(KPI_result_list)>1:
       # KPI.plot_KPI(KPI_result_list, scenarios_list, directory_path)
    
    #print("***Overall processing finished***")
    #print(f"Find all results and graphs in directory : {directory_path}")

if __name__ == "__main__":
    plt.close("all")
    
    import subprocess
    import time
    start_overall = time.time()
    print('Main code started')
    print('-----------------')
    
    main()    

    print('-----------------')
    print('Main code ended')
    print('-----------------')
    end_overall = time.time()
    duration_overall = end_overall - start_overall
    m, s = divmod(duration_overall, 60)
    print('Overall run time :', "%.0f" %m,'min', "%.0f" %s,'s \n')














