import pandas as pd
import geopandas as gpd
import numpy as np

import matplotlib.pyplot as plt

from shapely.ops import split, unary_union
from shapely.geometry import LineString, MultiPolygon, Point, Polygon


def points_to_coords(points):
    coords = []
    for point in points:
        if point is None or point.is_empty:
            raise ValueError("Voronoi points cannot be empty.")
        if point.geom_type != "Point":
            raise TypeError("Voronoi discretization expects Point geometries.")
        coords.append((point.x, point.y))
    return np.asarray(coords, dtype=float)


def _line_split_scale(boundary_shape):
    minx, miny, maxx, maxy = boundary_shape.bounds
    diagonal = np.hypot(maxx - minx, maxy - miny)
    return max(diagonal * 4, 1.0)


def _point_belongs_to_site(point_coords, site_coords, other_coords):
    return np.linalg.norm(point_coords - site_coords) <= np.linalg.norm(
        point_coords - other_coords
    )


def _clip_with_bisector(cell, site_coords, other_coords, split_scale):
    if np.allclose(site_coords, other_coords):
        return cell

    midpoint = (site_coords + other_coords) / 2
    normal = other_coords - site_coords
    direction = np.array([-normal[1], normal[0]], dtype=float)
    direction_norm = np.linalg.norm(direction)
    if direction_norm == 0:
        return cell

    direction = direction / direction_norm * split_scale
    bisector = LineString([midpoint - direction, midpoint + direction])
    parts = list(split(cell, bisector).geoms)
    if len(parts) <= 1:
        return cell

    kept_parts = []
    for part in parts:
        representative = part.representative_point()
        representative_coords = np.array([representative.x, representative.y])
        if _point_belongs_to_site(representative_coords, site_coords, other_coords):
            kept_parts.append(part)

    if not kept_parts:
        return cell

    return unary_union(kept_parts)


def voronoi_regions_from_coords(
    coords,
    boundary_shape,
    return_unassigned_points=False,
    per_geom=False,
):
    del per_geom

    if len(coords) == 0:
        result = ({}, [], [])
        if return_unassigned_points:
            return result
        return result[0], result[1]

    unique_coords, inverse = np.unique(coords, axis=0, return_inverse=True)
    split_scale = _line_split_scale(boundary_shape)

    region_polygons = {}
    for unique_idx, site_coords in enumerate(unique_coords):
        cell = boundary_shape
        for other_idx, other_coords in enumerate(unique_coords):
            if unique_idx == other_idx:
                continue
            cell = _clip_with_bisector(cell, site_coords, other_coords, split_scale)
            if cell.is_empty:
                break
        region_polygons[unique_idx] = cell

    point_regions = {point_idx: region_polygons[unique_idx] for point_idx, unique_idx in enumerate(inverse)}
    unassigned_points = []
    if return_unassigned_points:
        return point_regions, coords, unassigned_points
    return point_regions, coords


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

    poly_shapes, _, _ = voronoi_regions_from_coords(
        coords,
        boundary_shape,
        return_unassigned_points=True,
        per_geom=False,
    )

    voros = pd.DataFrame(
        {
            "id": list(poly_shapes.keys()),
            "geometry": list(poly_shapes.values()),
        }
    ).sort_values(by="id")

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
