import pandas as pd
import geopandas as gpd
import numpy as np

import matplotlib.pyplot as plt

from shapely.ops import unary_union
# from scipy.spatial import Voronoi
from shapely.geometry import Polygon, MultiPolygon, Point

from geovoronoi import voronoi_regions_from_coords, points_to_coords



def get_voronoi_regions(voro_points: gpd.GeoDataFrame, borders: gpd.GeoDataFrame):
    bases = voro_points[["geometry"]].copy()
    max_borders = borders[["geometry"]].copy()

    bases["id"] = list(range(len(bases)))
    if bases.crs is None:
        bases = bases.set_crs("epsg:4326")
    else:
        bases = bases.to_crs("epsg:4326")
        
    if max_borders.crs is None:
        max_borders = max_borders.set_crs(epsg=3395)
    else:
        max_borders = max_borders.to_crs(epsg=3395)

    bases_proj = bases.to_crs(max_borders.crs)


    boundary_shape = unary_union(max_borders.geometry)
    coords = points_to_coords(bases_proj.geometry)

    # gc = voronoi_polygons(bases["geometry"], extend_to=max_borders["geometry"])
    poly_shapes, pts, unassigned = voronoi_regions_from_coords(coords, boundary_shape, return_unassigned_points=True,per_geom=False)

    items = poly_shapes.items()
    bases_poly_map = {}
    for key,poly in poly_shapes.items():
        bases_poly_map[key] = []
        for i,row in bases_proj.iterrows():
            if poly.contains(row["geometry"]):
                bases_poly_map[key].append(i)

    voros = pd.DataFrame()
    voros["id"] = [x[0] for x in items]
    voros["geometry"] = [x[1] for x in items]
    voros["id"].replace(bases_poly_map,inplace=True)
    voros = voros.sort_values(by="id")


    voros = gpd.GeoDataFrame(voros, geometry="geometry")
    bases = bases.merge(voros, left_on="id", right_on="id")

    bases = bases.rename(columns={"geometry_y": "geometry"})
    bases = gpd.GeoDataFrame(bases[["geometry"]].copy(), geometry="geometry")
    bases = bases.set_crs(epsg=3395)
    bases = bases.to_crs("epsg:4326")

    return bases


def distance(p1, p2):
    R = 6370
    pi = np.pi
    lat1 = p1[0]
    long1 = p1[1]
    lat2 = p2[0]
    long2 = p2[1]
    cos = np.cos
    sin = np.sin
    asin = np.arcsin
    point1=(R*cos((pi/180)*lat1)*cos((pi/180)*long1), R*cos((pi/180)*lat1)*sin((pi/180)*long1), R*sin((pi/180)*lat1))
    point2=(R*cos((pi/180)*lat2)*cos((pi/180)*long2), R*cos((pi/180)*lat2)*sin((pi/180)*long2), R*sin((pi/180)*lat2))
    d=np.linalg.norm(np.array(point1)-np.array(point2))
    dearth = 2*R*asin(d/(2*R))


    return dearth