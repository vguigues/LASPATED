import h3
from geopandas import GeoDataFrame
from shapely.geometry import Polygon


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

    hex_indexes = set()
    for observation in h3_gdf.geometry:
        if observation.geom_type not in {"Polygon", "MultiPolygon"}:
            raise ValueError(
                "GeoDataFrame's geometry is limited to either polygon of multi-polygon. Got {}".format(
                    observation.geom_type
                )
            )
        hex_indexes.update(h3.geo_to_cells(observation.__geo_interface__, resolution))

    # we have the H3 indexes in hex_indexes
    # Now we just need to transform them into an geodataframe with any relevant info we may need
    hex_indexes = list(hex_indexes)
    polygons = []
    neighbors = []
    for hex in hex_indexes:
        # send a warning if there is a pentagon in the study region!
        if h3.is_pentagon(hex):
            print(
                "A H3 cell in the study region is a pentagon. See H3's documentation for further details."
            )

        boundary = h3.cell_to_boundary(hex)
        pol = Polygon([(lng, lat) for lat, lng in boundary])
        polygons.append(pol)
        hex_ring = h3.grid_ring(hex, k=1)
        neighbors.append(
            [
                hex_indexes.index(neighbor)
                for neighbor in hex_ring
                if neighbor in hex_indexes
            ]
        )

    # is there an other relevant info that could be calculated here?
    temp_dict = {"geometry": polygons, "neighbors": neighbors}

    hex_gdf = (
        GeoDataFrame(temp_dict, crs="EPSG:4326")
        .to_crs(gdf.crs)
        .reset_index()
        .rename({"index": "id"}, axis=1)
    )

    return hex_gdf
