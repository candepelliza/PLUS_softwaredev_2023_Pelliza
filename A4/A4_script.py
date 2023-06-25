
"""Geocoding and Isochrones dunctions

This script contains the definition of the functions used in the A4 notebook.
It can be imported as a module and requires the dependencies from the ox environment to be installed
Contained functions:

    * get_address - gives a pair of locations por a pair of coordinates corresponding to the intersections of a street
    * get_isochrones - gives a pair of polygons corresponding to the areas reachable in a given walking distance for each of both street intersections
    * polygon_disolve - disolves both polygons for each segment, obtaining an isochrone corresponding to the total area reachable in a given walkable distance from the segment
"""


#Import libraries
from geopy.geocoders import Nominatim
import geopandas as gpd
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon, MultiLineString
import shapely.geometry as sg
from shapely.ops import unary_union, split, linemerge
import networkx as nx
import osmnx as ox
from plotnine import ggplot, aes, geom_line, geom_map
import fiona
import pyproj
from pyproj import CRS
import warnings


def get_address(dataset):
    """This function gets the addresses from a dataset containing a pair of coordinates defininf the intersections of a street

    Args:
        dataset (dataFrame): dataframe following the predefined format

    Returns:
        DataFrame: existing dtaframe with two new columns containing the address for intersection 1 and intersection 2
    """
    # Create an empty list to store the addresses
    address1 = []
    address2 = []

    for index, i in dataset.iterrows():
        lat1 = i['intrsc_1_lat']
        long1 = i['intrsc_1_long']
        lat2 = i['intrsc_2_lat']
        long2 = i['intrsc_2_long']
        # Geolocate
        geolocator = Nominatim(user_agent="geocode_A4")
        coords1 = "{}, {}".format(lat1, long1)
        coords2 = "{}, {}".format(lat2, long2)
        location1 = geolocator.reverse(coords1)
        location2 = geolocator.reverse(coords2)
        # Append the addresses to the list
        address1.append(location1.address)
        address2.append(location2.address)
        
    # Assign the addresses list as a new column in the "data" DataFrame
    dataset['address1'] = address1
    dataset['address2'] = address2

    return dataset

#------------------------------------------------------------------------------------------------

def get_isochrones(dataset, meters):
    """Based on functions provided by the OSMNx library, the function obtains one polygon (isochrone) for each of the 2 segment intersections or endpoints. This polygons represent the reachable area within a given walking distance starting from the analysis intersect point.
    Args:
        dataset (DataFrame): dataframe following the predefined format
        meters (int): Isochrone walking distance in meters

    Returns:
        DataFrame: dataframe containing the id of the analyzed segment ("sgmntid") and the geometry of the isochrone polygons for both intersections ("isochrone1" and "isochrone2")
    """

    # Create empty lists
    polygon1 = []
    polygon2 = []
    sgmntids = []

    for index, row in dataset.iterrows():

        try:

            #ISOCRHONE INTERSECTION 1
            lat1 = row['intrsc_1_lat']
            lon1 = row['intrsc_1_long']
            loc1 = (lat1, lon1)
            #retrieves a walking network graph using the osmnx library.
            G1 = ox.graph_from_point(loc1, simplify=True, network_type='walk', dist = meters)
            #projects the graph onto the EPSG:4326
            G1 = ox.project_graph(G1, to_crs="4326")        
            #extracts the coordinates of each node in the subgraph and creates a Point object for each one using a list comprehension.
            node_points1 = [Point(data['x'], data['y']) for node, data in G1.nodes(data=True)]
            #creates a GeoSeries object from the list of Point objects and computes the convex hull of the GeoSeries object 
            polys1 = gpd.GeoSeries(node_points1).unary_union.convex_hull
            polygon1.append(polys1)

        except:
            #if it doesnt find any node in the given distance, and cannot create the polygon, it returns 0
            polygon1.append("0")

        try:
            # ISOCHRONE INTERSECTION 2
            lat2 = row['intrsc_2_lat']
            lon2 = row['intrsc_2_long']
            loc2 = (lat2, lon2)
            #retrieves a walking network graph using the osmnx library.
            G2 = ox.graph_from_point(loc2, simplify=True, network_type='walk', dist = meters)
            #projects the graph onto the EPSG:4326
            G2 = ox.project_graph(G2, to_crs="4326") # Use this line if the coordinates sistem returned from polys is changed from the original (check which crs you are using)
            #extracts the coordinates of each node in the subgraph and creates a Point object for each one using a list comprehension.
            node_points2 = [Point(data['x'], data['y']) for node, data in G2.nodes(data=True)]
            #creates a GeoSeries object from the list of Point objects and computes the convex hull of the GeoSeries object 
            polys2 = gpd.GeoSeries(node_points2).unary_union.convex_hull
            polygon2.append(polys2) 

        except:
            #if it doesnt find any node in the given distance, and cannot create the polygon, it returns 0
            polygon2.append("0") 

        sgmntid = row['sgmntid']
        sgmntids.append(sgmntid)

    # Add columns to dataset
    dataset['isochrone1'] = polygon1
    dataset['isochrone2'] = polygon2
    
    return dataset


#---------------------------------------------------------------------------------------------------------

# POLYGON_DISOLVE: 

def polygon_disolve(dataset):
    """Merges intersection 1 and intersection 2 polygons obtained with the get_isochrones() function, to get one single polygon representing the whole area reachable in the given distance from any point of the analyzed segment

    Args:
        dataset (DataFrame): dataframe obtained as output from the get_isochrones() function execution

    Returns:
        GeoDataFrame: geo dataframe containing the isochrone for the given distance for each analyzed segment
    """

    # Iterate over each row in the DataFrame
    for index, row in dataset.iterrows():
        # Get the values of 'isochrone1' and 'isochrone2' for the current row
        isochrone1 = row['isochrone1']
        isochrone2 = row['isochrone2']

        # Check the types of 'isochrone1' and 'isochrone2'
        if isinstance(isochrone1, sg.Polygon) and isinstance(isochrone2, sg.Polygon):
            # Both 'isochrone1' and 'isochrone2' are polygons, no action needed
            pass
        elif isinstance(isochrone1, sg.Polygon) and not isinstance(isochrone2, sg.Polygon):
            # 'isochrone1' is a polygon, 'isochrone2' is not, so copy 'isochrone1' to 'isochrone2'
            row['isochrone2'] = row['isochrone1']
        elif isinstance(isochrone2, sg.Polygon) and not isinstance(isochrone1, sg.Polygon):
            # 'isochrone2' is a polygon, 'isochrone1' is not, so copy 'isochrone2' to 'isochrone1'
            row['isochrone1'] = row['isochrone2']
        else:
            # Neither 'isochrone1' nor 'isochrone2' is a polygon, no action needed
            pass

    # Filter out rows where 'isochrone1' or 'isochrone2' is '0'
    isochrones_df = dataset[dataset['isochrone1'] != '0']
    isochrones_df = dataset[dataset['isochrone2'] != '0']

    # Perform unary union operation on 'isochrone1' and 'isochrone2' for each row
    for index, row in dataset.iterrows():
        dataset['geometry'] = dataset.apply(lambda row: unary_union([row['isochrone1'], row['isochrone2']]), axis=1)
        dataset['geometry'] = dataset['geometry'].apply(unary_union)

    # Create a GeoDataFrame with the 'geometry' column from 'dataset'
    isochrones_sgmnt_df = gpd.GeoDataFrame(geometry=dataset['geometry'])

    # Update the 'geometry' column in 'dataset' with the GeoDataFrame
    dataset['geometry'] = isochrones_sgmnt_df

    # Filter the desired columns from 'dataset'
    isochrones_dis = dataset.filter(items=['sgmntid', 'address1', 'address2', 'geometry'])

    # Create a GeoDataFrame with the filtered data and set the CRS
    isochrones_dis = gpd.GeoDataFrame(isochrones_dis, geometry='geometry', crs="EPSG:4326")

    # Return the final GeoDataFrame
    return isochrones_dis
