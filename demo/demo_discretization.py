import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
import numpy as np

import laspated as spated
from shapely.geometry import Polygon, MultiPolygon, Point
import sys

if len(sys.argv) != 3:
    raise ValueError(f"demo usage: python demo.py <region> <geo_disc>\n\twhere <region> in ['rectangle', 'convex', 'custom'] and <geo_disc> in ['rectangle','hexagon','custom']")    
    
region_type = sys.argv[1]
disc_type = sys.argv[2]

if not region_type in ['rectangle', 'convex', 'custom']:
    raise ValueError(f"<region> must be in ['rectangle', 'convex', 'custom']")

if not disc_type in ['rectangle','hexagon','custom']:
    raise ValueError(f"<disc_type> must be in ['rectangle','hexagon','custom']")

app = spated.DataAggregator(crs="epsg:4326")
events = pd.read_csv(r"../Data/sorted_events.csv", encoding="ISO-8859-1", sep=",")
app.add_events_data(
    events,
    datetime_col="data_hora",
    lat_col="lat",
    lon_col="long",
    feature_cols=["prioridade"],
    datetime_format="%m/%d/%y %H:%M:%S",
)  # %m/%d/%y %H:%M:%S

if region_type == "custom":
    max_borders = gpd.read_file(
        r"../Data/rj/rj.shp"
    )  # Load the geometry of region of interest
    app.add_max_borders(max_borders)
else:
    app.add_max_borders(method=region_type)

# Time discretizations
app.add_time_discretization("m", 30, 60 * 24, column_name="hhs")

app.add_time_discretization("D", 1, 7, column_name="dow")

time_disc_df = pd.DataFrame([
    ["2016-01-01", "2016-01-01", 1, "yearly"],
    ["2016-02-06", "2016-02-11", 2, None],
    ["2017-02-24", "2017-03-06", 2, None],
], columns=["start", "end", "t", "repetition"])

app.add_time_discretization(time_disc_df)

#Geo discretization
if disc_type == "rectangle":
    app.add_geo_discretization(
        discr_type="R", rect_discr_param_x=10, rect_discr_param_y=10
    )
elif disc_type == "hexagon":
    app.add_geo_discretization(
        discr_type='H',
        hex_discr_param=7
    ) 
else:
    custom_map = gpd.read_file(r'../Data/rio_de_janeiro_neighborhoods/rio_neighborhoods.shp')
    app.add_geo_discretization('C', custom_data=custom_map)

# Plot discretization
app.plot_discretization()

# Geo Features
population = gpd.read_file(r"../Data/regressores/populacao/")
population = population[["population", "geometry"]].copy()
app.add_geo_variable(population)

print("Writing output files...")
# Writing output
app.write_arrivals("disc_data/arrivals.dat")
app.write_regions("disc_data/neighbors.dat")
app.write_info(obs_index_column="dow", path="disc_data/info.dat")