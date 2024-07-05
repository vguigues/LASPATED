import h3
from geopandas import GeoDataFrame
from shapely.geometry import Polygon, Point


# helper function: gets geojson like (H3 expects as input) from shapely polygon
def polygon_to_geojson(polygon: Polygon):
    """
    Auxiliar function: Gets geojson like (H3 expect as input) from
    Shapely Polygon object
    """
    temp_coord_list = list(polygon.exterior.coords)
    # we have to manually convert from a list of list of tuples
    # to a list of list of lists (and also invert lat, long ordering)
    coords = [[lat, long] for (long, lat) in temp_coord_list]
    # we dont want the loop around here
    coords.pop()
    # generate geojson
    geoJson = {"type": "Polygon", "coordinates": [coords]}

    return geoJson


def generate_H3_discretization(gdf: GeoDataFrame, resolution: int = 7):
    """
    Generate a hexagonal discretization of the area using Uber's H3 library
    and returns a new GeoDataFrame with it.

    Parameters
    ----------
    param gdf: GeoDataFrame
        Original geodataframe whose boundaries will be considered to generate the discretization
    param resolution: int
        Desired resolution level passed to the H3 library. A bigger number means smaller hexagons.
        Valid range : [0,15]

    Returns
    -------
    return: GeoDataFrame
        New GeoDataFrame containing the H3 hexagons that approximately cover the original region
    """

    if resolution < 0 or resolution > 15:
        raise ValueError(
            "Hex resolution must be in range [0, 15]. Got {}".format(resolution)
        )

    # following code assumes crs epsg = 4326 to interface with H3
    h3_gdf = gdf.to_crs(epsg=4326)

    # loop through all observation in the geoseries
    # For every polygon there, compute using h3 the a list of h3 indexes
    # and add it to the hex_indexes structure, used later to construct a new goodataframe

    # this procedure "flattens" the original geoseries, ie, multipolygon's are
    # treated as a sequence of polygon without differentiating them in any way

    hex_indexes = set()
    for observation in h3_gdf.geometry:
        if observation.geom_type == "MultiPolygon":
            for polygon in observation.geoms:
                geoJson = polygon_to_geojson(polygon)
                # h3.polyfill is the important method in the H3 library
                # that does the heavy work of finding a good hex-cover
                hex_indexes.update(h3.polyfill(geoJson, resolution))
        elif observation.geom_type == "Polygon":
            geoJson = polygon_to_geojson(observation)
            hex_indexes.update(h3.polyfill(geoJson, resolution))
        else:
            raise ValueError(
                "GeoDataFrame's geometry is limited to either polygon of multi-polygon. Got {}".format(
                    observation.geom_type
                )
            )

    # we have the H3 indexes in hex_indexes
    # Now we just need to transform them into an geodataframe with any relevant info we may need
    hex_indexes = list(hex_indexes)
    polygons = []
    pol_areas = []
    neighbors = []
    center_points = []
    c_neighbors = [[], [], [], [], [], []]
    for hex in hex_indexes:
        # send a warning if there is a pentagon in the study region!
        if h3.h3_is_pentagon(hex):
            print(
                "A H3 cell in the study region is a pentagon. See H3's documentation for further details."
            )

        # once again: h3 uses 2-tuples for points, but shapely uses 2-lists
        hex_coords_h3 = h3.h3_set_to_multi_polygon([hex], geo_json=False)
        hex_coords_sh = []
        for coords in hex_coords_h3[0]:
            for coord in coords:
                # and in a different (lat,long) order
                hex_coords_sh.append([coord[1], coord[0]])

        pol = Polygon(hex_coords_sh)
        polygons.append(pol)
        # cell area in default km2
        pol_areas.append(h3.cell_area(hex))

        # collects h3 center points. There might be a very slight difference between these points
        # and the ones calculated by shapely based on the polygons!
        lat_long_center = h3.h3_to_geo(hex)
        center_points.append(Point(lat_long_center[1], lat_long_center[0]))

        hex_ring = h3.hex_ring(hex, k=1)
        neighbors.append(
            [
                hex_indexes.index(neighbor)
                for neighbor in hex_ring
                if neighbor in hex_indexes
            ]
        )
        i_n = 0
        i_filled = 0
        for neighbor in hex_ring:
            if neighbor in hex_indexes:
                c_neighbors[i_filled].append(neighbor)
                i_filled += 1
            i_n += 1
        assert (
            i_n == 6
        ), "Cell found that did not have 6 neighbors. Does the study area contain H3 pentagons? Check H3 official documentation for details"
        while i_filled < i_n:
            c_neighbors[i_filled].append(None)
            i_filled += 1

    # is there an other relevant info that could be calculated here?
    temp_dict = {"geometry": polygons, "neighbors": neighbors}

    hex_gdf = (
        GeoDataFrame(temp_dict, crs="EPSG:4326")
        .to_crs(gdf.crs)
        .reset_index()
        .rename({"index": "id"}, axis=1)
    )

    return hex_gdf
