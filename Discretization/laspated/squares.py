import geopandas as gpd
from shapely.geometry import Polygon


def rectangle_discretization(
    gdf: gpd.GeoDataFrame,
    nx: int,
    ny: int,
    neighborhood: int = 8
):
    '''
    Construct a square discretized gdf from the original region. The final
    discretization is intended to be a grid of nx per ny rectangles.
    Also construct neighbors data for each rectangle.

    Please be aware that each cell in the returned GeoDataFrame is not
    necessarily a rectangle: edge regions or 'islands' will fit to their
    original formats.  A data structural consequence of this is that the
    shapes forming the regions are not all Polygons.
    They may be Multi-polygons, or Points / Multipoint depending on the
    original topology.

    Arguments:
        gdf : GeoDataFrame
            The GeoDataFrame describing the original study region
        nx, ny
            Define how many subdivision to use on the discretization,
            on the x and y axis respectively
    '''
    # check for negative discretization sizes
    if nx <= 0 or ny <= 0:
        raise ValueError("Squares discretization sizes must be positive")
    # get bounds
    minx, miny, maxx, maxy = gdf.geometry.total_bounds
    # calculate square sizes
    deltay = (maxy - miny) / ny
    deltax = (maxx - minx) / nx

    # construct 'squares' GeoDataFrame
    # Lots of rectangles covering the entire study region
    squares = []
    indices = []
    neighbors = []
    # valid_index = lambda ix, iy : True if ix >= 0 and iy >= 0 and ix < nx
    # and iy < ny else False
    # index = lambda ix, iy : iy * ny + ix if valid_index(ix, iy) else None

    def get_index(ix, iy):
        '''
        Auxiliary function that gets square list index
        based on position (ix, iy)
        '''
        if ix >= 0 and iy >= 0 and ix < nx and iy < ny:
            return iy * ny + ix
        else:
            return None

    for iy in range(ny):
        for ix in range(nx):
            cx, cy = minx + ix * deltax, miny + iy * deltay
            pol = Polygon([
                (cx, cy),
                (cx + deltax, cy),
                (cx + deltax, cy + deltay),
                (cx, cy + deltay)
            ])
            squares.append(pol)
            indices.append(get_index(ix, iy))

            if neighborhood == 4:
                my_neighbors = [
                    get_index(ix + 1, iy),
                    get_index(ix - 1, iy),
                    get_index(ix, iy + 1),
                    get_index(ix, iy - 1)
                ]
            elif neighborhood == 8:
                my_neighbors = [
                    get_index(ix + 1, iy),
                    get_index(ix - 1, iy),
                    get_index(ix, iy + 1),
                    get_index(ix, iy - 1),
                    get_index(ix + 1, iy + 1),
                    get_index(ix + 1, iy - 1),
                    get_index(ix - 1, iy + 1),
                    get_index(ix - 1, iy - 1)
                ]
            else:
                raise ValueError("Invalid neighborhood type " + neighborhood)
            my_neighbors = [n for n in my_neighbors if n is not None]
            neighbors.append(my_neighbors)

    squares = gpd.GeoSeries(squares)
    squares_gdf = gpd.GeoDataFrame({
        'geometry': squares,
        'myindex': indices,
        'neighbors': neighbors
        }, crs=gdf.crs)

    # calculate the intersection between this squares gdf and
    # the original region
    res_intersection = gpd.overlay(squares_gdf, gdf, how='intersection')

    # When calculating the intersection, some rectangles may have been
    # dropped entirely. Thus, we need to adjust the 'neighbor' column
    # to correctly reflect the new indices. That's why we keep a auxiliary
    # 'myindex' column. We just need to adjust 'myindex' to actual indexes.
    new_neighbors = []
    for i, row in res_intersection.iterrows():
        new_list = [
            res_intersection.index[res_intersection['myindex'] == n]
            for n in row['neighbors']
        ]
        new_list = [item for sublist in new_list for item in sublist]
        new_neighbors.append(new_list)
    res_intersection['neighbors'] = new_neighbors

    res_intersection.drop(['myindex'], axis=1, inplace=True)
    res_intersection.reset_index(inplace=True)
    res_intersection.rename(columns={'index': 'id'}, inplace=True)
    return res_intersection[['id', 'neighbors', 'geometry']]
