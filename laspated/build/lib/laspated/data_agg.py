import numpy as np
import pandas as pd
import geopandas as gpd

from typing import List
from shapely.ops import unary_union
from shapely.geometry import Polygon, MultiPolygon, Point

from geovoronoi.plotting import subplot_for_map, plot_voronoi_polys_with_points_in_area
from geovoronoi import voronoi_regions_from_coords, points_to_coords

from .time_discretization_utils import calculate_seasonality, apply_custom_time_events
from .geo_discretization_utils import get_voronoi_regions, distance
from .squares import rectangle_discretization
from .h3_utils import generate_H3_discretization
from .add_regressors import addRegressorUniformDistribution, addRegressorWeightedAverage


class DataAggregator:

    def __init__(self, crs: str = "epsg:4326"):
        """
        Class initialization.

        Parameters
        ----------
        crs : The value can be any accepted by pyproj.CRS.from_user_input()
            Coordinate Reference System to be assured and applied.
        """
        self.crs = crs
        self.time_indexes = []
        self.geo_index = "gdiscr"
        self.events_features = []
        self.geo_features = []

        self.max_borders = None
        self.events_data = None
        self.geo_discretization = None

    def add_max_borders(self, data: gpd.GeoDataFrame = None, method: str = None):
        """
        Adds limit borders to internal map. Will be used on multiple
        methods to limit events and discretization reach.

        Parameters
        ----------
        data : geopandas.GeoDataFrame
            GeoDataFrame containing Polygon representing limit borders

        method : str
            Method of estimating borders using events data:

                rectangle => Calculates smallest rectangle containing all event points
                convex => Calculates aproximate of smallest convex polygon containing all event points
        """

        if data is not None:
            # If no CRS set, will assume class set CRS
            if data.crs is None:
                print(
                    f"""Limit borders data is not set with CRS information. Will assume "{self.crs}"."""
                )
                self.max_borders = data

            # If CRS identified, set new one if needed
            elif data.crs != self.crs:
                self.max_borders = data.to_crs(self.crs)

            else:
                self.max_borders = data

        elif method == "rectangle":

            # if no events data was passed, raise error
            if not (hasattr(self, "events_data")):
                raise AttributeError(
                    "Please inform events data using `add_events_data` method."
                )

            # get limits and compute rectangle
            minx, miny, maxx, maxy = self.events_data.geometry.total_bounds
            pol_max_borders = Polygon(
                [(minx, miny), (minx, maxy), (maxx, maxy), (maxx, miny)]
            )
            # replace self.max_borders
            self.max_borders = gpd.GeoDataFrame()
            self.max_borders["geometry"] = [pol_max_borders]
            self.max_borders = self.max_borders.set_crs(self.crs)

        elif method == "convex":

            # if no events data was passed, raise error
            if not (hasattr(self, "events_data")):
                raise AttributeError(
                    "Please inform events data using `add_events_data` method."
                )

            x_std = self.events_data.geometry.x.std() * 0.1
            y_std = self.events_data.geometry.y.std() * 0.1
            # concatenate list of mini event-polygons
            d_polys = self.events_data.drop(
                self.events_data.loc[self.events_data.geometry.x.isna()].index
            ).geometry.apply(
                lambda pt: (
                    Polygon(
                        [
                            (pt.x - x_std, pt.y - y_std),
                            (pt.x - x_std, pt.y + y_std),
                            (pt.x + x_std, pt.y + y_std),
                            (pt.x + x_std, pt.y - y_std),
                        ]
                    )
                    if pt.x >= -np.inf
                    else pt
                )
            )
            # replace self.max_borders
            self.max_borders = gpd.GeoDataFrame()
            self.max_borders["geometry"] = [
                MultiPolygon(list(d_polys.values)).convex_hull
            ]
            self.max_borders = self.max_borders.set_crs(self.crs)

        else:
            raise ValueError("Please inform `data` or `method` parameter.")

    def add_events_data(
        self,
        events_data: pd.DataFrame,
        datetime_col: str = None,
        lat_col: str = "lat",
        lon_col: str = "lon",
        feature_cols: List[str] = [],
        datetime_format: str = None,
    ):
        """
        Processes and save events registers.

        Parameters
        ----------
        events_data : pandas.DataFrame, geopandas.GeoDataFrame
            DataFrame containing latitude longitude information
            -> If pandas dataframe, should contain two
            columns (lat_col, lon_col) with latitude and longitude
            coordinates
            -> If geopandas geodataframe, should contain "geometry" column

        datetime_col : str
            Column containing time information from events
            Can be not informed if not applicable

        lat_col : str
            Column containing latitude from events

        lon_col : str
            Column containing longitude from events

        feature_cols : List[str]
            Column(s) containing adicional information from events
            that will be considered in discretization process

        datetime_format : str
            The strftime to parse time, e.g. "%d/%m/%Y".
            See https://docs.python.org/3/library/datetime.html#strftime-and-strptime-behavior
            for reference
        """
        # configure geodataframe
        # if pandas dataframe type, will create geometry columns
        if type(events_data) == pd.DataFrame:
            events_data["geometry"] = gpd.points_from_xy(
                events_data[lon_col], events_data[lat_col]
            )
            self.events_data = (
                gpd.GeoDataFrame(events_data).set_crs("epsg:4326").to_crs(self.crs)
            )
        # if geopandas object, will assure correct CRS
        else:
            if events_data.crs is None:
                self.events_data = events_data.set_crs(self.crs)
            else:
                self.events_data = events_data.to_crs(self.crs)

        # transform datetime_col data to datetime
        if datetime_col is not None:
            ts = self.events_data[datetime_col]
            if ts.dtype == "O":
                self.events_data["ts"] = pd.to_datetime(ts, format=datetime_format)
            else:
                self.events_data["ts"] = ts.copy()
            self.events_data.sort_values("ts", ascending=True, inplace=True)

        # filter only necessary cols
        self.events_features += feature_cols
        special_cols = ["geometry", "ts"] if datetime_col is not None else ["geometry"]
        self.events_data = self.events_data[feature_cols + special_cols].copy()

    def add_time_discretization(
        self,
        seasonality_type: str,
        window: int = 1,
        frequency: int = 1,
        column_name: str = None,
    ):
        """
        Applies seasonality index as time discretization new column.

        Parameters
        ----------
        seasonality_type: str or DataFrame
            If string, represents seasonality frequency type
                Y -> year
                M -> month
                W -> week
                D -> day
                H -> hour
                m -> minute
                S -> second

            If dataframe, represents all necessary information for each event.

        window: (List[int], int)
            Seasonality window. Can be unique, passed as an int, or variable,
            passed as a list of int.

        frequency: int
            Seasonality frequency.

            eg: If seasonality is between 12 months in a year:
                seasonality_type = 'Y'
                window = 1
                frequency = 12

                If time discretization must consider the first 12 hours in a day as first index,
                next 8 hours as second and last 4 hours as third:
                seasonality_type = 'D'
                window = [12, 8, 4]
                frequency = 24
        """
        # throw error if events_data is None
        if self.events_data is None:
            raise AttributeError(
                "Please inform events data using `add_events_data` method."
            )

        # get new column name
        if column_name is None:
            i_col = 1
            while True:
                t_col = f"tdiscr_{i_col}"
                if t_col in self.events_data.columns:
                    i_col += 1
                else:
                    break
        else:
            t_col = column_name

        # if dataframe passed, apply custom events time discretization
        if type(seasonality_type) == pd.DataFrame:
            self.time_indexes.append(t_col)
            self.events_data[t_col] = apply_custom_time_events(
                ts=self.events_data["ts"], time_disc_df=seasonality_type, nan_idx=0
            )

        else:
            # throw error if frequency is not compatible with window
            if (type(window) is int) and (frequency % window != 0):
                raise ValueError(
                    "Parameter `frequency` must be a multiple of `window`."
                )
            if type(window) is list and sum(window) != frequency != 0:
                raise ValueError("Parameter `frequency` must be equal sum of `window`.")

            # add new column
            self.time_indexes.append(t_col)
            self.events_data[t_col] = calculate_seasonality(
                ts=self.events_data["ts"],
                seasonality_type=seasonality_type,
                window=window,
                frequency=frequency,
            )
        self.events_data = self.events_data.dropna()

    def add_geo_discretization(
        self,
        discr_type: str,
        hex_discr_param: int = None,
        rect_discr_param_x: int = None,
        rect_discr_param_y: int = None,
        custom_data: gpd.GeoDataFrame = None,
    ) -> gpd.GeoDataFrame:
        """
        Calculates geographical discretization and apply index
        to events data.

        Parameters
        ----------
        discr_type : str
            Indicator of method of geographical discretization.

                R => Rectangular discretization. Parameters `rect_discr_param_x`
                and `rect_discr_param_y` are required.

                H => Hexagonal discretization using H3 package. Parameter
                `hex_discr_param` is required.

                C => Custom discretization using user GeoDataFrame. Parameter
                `custom_data` is required.

                G => Custom discretization using graph data. Geolocated nodes
                must be informed in dataframe format `custom_data`.

                V => Voronoi discretization using GeoDataFrame with points.
                `custom_data` must contain the points that will define the voronoi regions.

        hex_discr_param : int
            Granularity level for hexagonal discretization.

        rect_discr_param_x : int
            Granularity level for rectangular horizontal discretization.

        rect_discr_param_y : int
            Granularity level for rectangular vertical discretization.
        """

        # if there is no max borders set throw error
        if self.max_borders is None:
            raise ValueError(
                "Please inform geographical limits using `add_max_borders` method."
            )

        # square discretization
        if discr_type == "R":
            self.geo_discretization = rectangle_discretization(
                gdf=self.max_borders, nx=rect_discr_param_x, ny=rect_discr_param_y
            )

        # hexagonal discretization
        elif discr_type == "H":
            full_area = self.max_borders.copy()  # keep original
            full_area.geometry = (
                full_area.geometry.convex_hull
            )  # apply convex hull so no area is cutted out
            self.geo_discretization = generate_H3_discretization(
                self.max_borders, hex_discr_param
            )

            # cut hexagons borders that dont belong to original shape
            limited_discretization = (
                gpd.overlay(
                    self.geo_discretization[["id", "geometry"]],
                    self.max_borders[["geometry"]],
                    how="intersection",
                )
                .dissolve(by="id")
                .reset_index()
            )
            self.geo_discretization = pd.merge(
                self.geo_discretization.drop("geometry", axis=1),
                limited_discretization,
                on="id",
                how="right",
            )
            # return type to geoDataFrame
            self.geo_discretization = gpd.GeoDataFrame(self.geo_discretization)

            # correct polygons index, ordering in neighbors list
            corr_index_dict = (
                self.geo_discretization.reset_index().set_index("id")["index"].to_dict()
            )
            neighbors_list = list(self.geo_discretization["neighbors"])
            for ind, n_list in enumerate(neighbors_list):
                new_list = []
                for n in n_list:
                    if corr_index_dict.get(n):
                        new_list.append(corr_index_dict[n])
                neighbors_list[ind] = new_list
            self.geo_discretization["neighbors"] = neighbors_list
            self.geo_discretization.drop("id", axis=1, inplace=True)
            self.geo_discretization.reset_index(inplace=True)
            self.geo_discretization.rename({"index": "id"}, axis=1, inplace=True)

        elif discr_type == "C":
            # check for custom_data
            if custom_data is None:
                raise ValueError(
                    "Please inform `custom_data` to be used as geographical discretization."
                )
            else:
                self.geo_discretization = custom_data[["geometry"]].copy()
                self.geo_discretization["id"] = list(range(len(custom_data)))

                # If no CRS set, will assume class set CRS
                if self.geo_discretization.crs is None:
                    print(
                        f'Custom data is not set with CRS information. Will assume "{self.crs}".'
                    )
                    self.geo_discretization = self.geo_discretization.set_crs(self.crs)

                # If CRS indetified, set new one if needed
                elif self.geo_discretization.crs != self.crs:
                    self.geo_discretization = self.geo_discretization.to_crs(self.crs)

            # compute neighbors
            neighbors = []
            for _, row in self.geo_discretization.iterrows():
                row_neighbors = list(
                    self.geo_discretization[
                        ~self.geo_discretization.geometry.disjoint(row.geometry)
                    ]["id"]
                )
                row_neighbors = [n for n in row_neighbors if n != row["id"]]
                neighbors.append(row_neighbors)
            self.geo_discretization["neighbors"] = neighbors

        elif discr_type == "G":
            # check for custom_data
            if custom_data is None:
                raise ValueError(
                    "Please inform `custom_data` to be used as grpah discretization."
                )
            else:
                self.geo_discretization = custom_data[["geometry"]].copy()
                self.geo_discretization["id"] = list(range(len(custom_data)))

                # If no CRS set, will assume class set CRS
                if self.geo_discretization.crs is None:
                    print(
                        f'Custom data is not set with CRS information. Will assume "{self.crs}".'
                    )
                    self.geo_discretization = self.geo_discretization.set_crs(self.crs)

                # If CRS indetified, set new one if needed
                elif self.geo_discretization.crs != self.crs:
                    self.geo_discretization = self.geo_discretization.to_crs(self.crs)

            # join events with nodes using nearest node
            self.events_data = (
                gpd.sjoin_nearest(
                    self.events_data.to_crs("epsg:29193"),
                    self.geo_discretization.to_crs("epsg:29193")[["id", "geometry"]],
                    how="left",
                )
                .to_crs(self.crs)
                .drop("index_right", axis=1)
                .rename(columns={"id": "node_id"})
            )

            # create map from node_id to gdiscr index
            node_idx_dict = dict(
                enumerate(self.events_data["node_id"].dropna().unique())
            )
            node_idx_dict = {v: k for k, v in node_idx_dict.items()}
            self.events_data["gdiscr"] = self.events_data["node_id"].map(node_idx_dict)
            self.events_data.drop("node_id", axis=1, inplace=True)

            # drop duplicated events
            self.events_data = (
                self.events_data.reset_index()
                .drop_duplicates(subset=["index"], keep="first")
                .set_index("index")
            )

            # not consider events outside max borders
            self.events_data.loc[
                ~self.events_data["geometry"]
                .dropna()
                .within(self.max_borders["geometry"].values[0]),
                "gdiscr",
            ] = pd.NA

            # column gdiscr is alread copmuted so function must be terminated
            return None
        elif discr_type == "V":
            if custom_data is None:
                raise ValueError(
                    "Please inform `custom_data` to be used as geographical discretization."
                )
            else:
                disc = get_voronoi_regions(custom_data, self.max_borders)
                self.geo_discretization = disc[["geometry"]].copy()
                self.geo_discretization["id"] = list(range(len(custom_data)))
                # self.geo_discretization = custom_data[['geometry']].copy()
                # self.geo_discretization['id'] = list(range(len(custom_data)))

                # If no CRS set, will assume class set CRS
                if self.geo_discretization.crs is None:
                    print(
                        f'Custom data is not set with CRS information. Will assume "{self.crs}".'
                    )
                    self.geo_discretization = self.geo_discretization.set_crs(self.crs)

                # If CRS identified, set new one if needed
                elif self.geo_discretization.crs != self.crs:
                    self.geo_discretization = self.geo_discretization.to_crs(self.crs)

            # compute neighbors
            neighbors = []
            for _, row in self.geo_discretization.iterrows():
                row_neighbors = list(
                    self.geo_discretization[
                        ~self.geo_discretization.geometry.disjoint(row.geometry)
                    ]["id"]
                )
                row_neighbors = [n for n in row_neighbors if n != row["id"]]
                neighbors.append(row_neighbors)
            self.geo_discretization["neighbors"] = neighbors
        else:
            raise ValueError(f"Invalid `discr_type` value {discr_type}.")

        # fills center_lat and center_lon (using projected crs for distance precision)
        centroids = self.geo_discretization.to_crs("epsg:4088").geometry.centroid
        self.geo_discretization["center_lat"] = centroids.x
        self.geo_discretization["center_lon"] = centroids.y

        # apply index to events data
        self.events_data = (
            gpd.sjoin(
                self.events_data.drop("gdiscr", axis=1, errors="ignore"),
                self.geo_discretization[["geometry", "id"]],
                how="left",
                predicate="within",
            )
            .drop("index_right", axis=1)
            .rename({"id": "gdiscr"}, axis=1)
        )

        self.events_data = self.events_data.dropna()

    def add_geo_variable(
        self, data: gpd.GeoDataFrame, type_geo_variable: str = "feature"
    ):
        """
        Merge information from external geographical data with computed
        geographical discretization.

        Parameters
        ----------
        data: geopandas.GeoDataFrame
            Geodataframe with desired variables.
        """

        if type_geo_variable != "feature" and type_geo_variable != "area":
            raise ValueError(
                'Parameter `type_geo_variable` must be "feature" or "area".'
            )

        # If no CRS set, will assume class set CRS
        if data.crs is None:
            print(
                f'Regressor data is not set with CRS information. Will assume "{self.crs}".'
            )
            regr_data = data.set_crs(self.crs)

        # If CRS indetified, set new one if needed
        else:
            regr_data = data.to_crs(self.crs)
        # add features to geo_features list
        self.geo_features += list(data.drop(["geometry"], axis=1).columns)

        # Calculate intersection between regressor data and geo discretization.
        # For math precision, if current crs is not projected,
        # will transform into proper (projected) CRS before using area metric
        # to calcultate intersections
        geo_disc = self.geo_discretization[["id", "geometry"]].copy()
        if not (self.geo_discretization.crs.is_projected):
            geo_disc = geo_disc.to_crs("epsg:4088")
            regr_data = regr_data.to_crs("epsg:4088")

        regr_intersection = addRegressorUniformDistribution(
            geo_disc, regr_data, discr_id_col="id", type_geo_variable=type_geo_variable
        )
        # regr_intersection = addRegressorWeightedAverage(
        #     geo_disc,
        #     regr_data,
        #     discr_id_col='id'
        # )
        self.geo_discretization = pd.merge(
            self.geo_discretization,
            regr_intersection.drop("geometry", axis=1),
            on="id",
            how="left",
        )

    def write_files(self, geo_discretization_path: str, events_data_path: str):
        """Write discretization DataFrames into files provided by geo_discretization_path and events_data_path.
            Files are saved as csv.
        Parameters
        ----------
        geo_discretization_path: str
            path containing location where to save the geo_discretization DataFrame.
        events_data_path: str
            path containing location where to save the events DataFrame.
        """
        self.geo_discretization.to_csv(geo_discretization_path, index=False)
        self.events_data.to_csv(events_data_path, index=False)

    def get_intersection(self, df_geo1: gpd.GeoDataFrame, df_geo2: gpd.GeoDataFrame):
        """
            Calculates the intersection between GeoDataFrames df_geo1 and df_geo2
        Parameters
        ----------
        df_geo1: gpd.GeoDataFrame
        df_geo2: gpd.GeoDataFrame
        """
        crs1 = df_geo1.crs
        crs2 = df_geo2.crs
        df_geo1 = df_geo1.assign(row_number_1=range(len(df_geo1)))
        df_geo2 = df_geo2.assign(row_number_2=range(len(df_geo2)))
        df_geo1 = df_geo1.to_crs("epsg:32723")
        df_geo2 = df_geo2.to_crs("epsg:32723")
        intersec = df_geo1.overlay(df_geo2, how="intersection")

        A = np.zeros((len(df_geo1), len(df_geo2)))
        intersec = intersec.to_crs("epsg:32723")
        for k, row in intersec.iterrows():
            i = int(row["row_number_1"])
            j = int(row["row_number_2"])
            A[i, j] += row["geometry"].area / 10**6

        df_geo1.drop(columns=["row_number_1"])
        df_geo2.drop(columns=["row_number_2"])
        df_geo1 = df_geo1.to_crs(crs1)
        df_geo2 = df_geo1.to_crs(crs2)
        return A

    def get_events_aggregated(self):
        """Returns np.array of lists with number of events given by time indexes, geodiscretization index and features indexes for each observation."""
        limits = []
        for time_index in self.time_indexes:
            limits.append(np.max(self.events_data[time_index]) + 1)
        limits.append(int(np.max(self.events_data["gdiscr"]) + 1))
        for feature in self.events_features:
            limits.append(np.max(self.events_data[feature]) + 1)

        samples = np.empty(tuple(limits), dtype=object)
        for index, val in np.ndenumerate(samples):
            samples[index] = []

        i = 0
        self.events_data.dropna()
        while i < len(self.events_data):
            row = self.events_data.iloc[i]
            times_i = [row[time_index] for time_index in self.time_indexes]
            index_i = (
                times_i
                + [int(row["gdiscr"])]
                + [row[feature] for feature in self.events_features]
            )
            count = 1
            for j in range(i + 1, len(self.events_data)):
                row_j = self.events_data.iloc[j]
                times_j = [row_j[time_index] for time_index in self.time_indexes]
                index_j = (
                    times_j
                    + [int(row_j["gdiscr"])]
                    + [row_j[feature] for feature in self.events_features]
                )

                if index_i == index_j:
                    count += 1
                elif times_i != times_j:
                    break
            samples[tuple(index_i)].append(count)
            i += count

        return samples

    def write_info(self, obs_index_column, path="info.dat"):
        """Write file with informations about given discretization.

        Args:
            obs_index_column (str): time index that defines the number of observations.
            path (str, optional): path of the written file. Defaults to "info.dat".
        """
        info_file = open(path, "w")
        for time_index in self.time_indexes:
            info_file.write(f"{np.max(self.events_data[time_index]) + 1} ")
        info_file.write(f"{int(np.max(self.events_data['gdiscr']) + 1)} ")
        for feature in self.events_features:
            info_file.write(f"{np.max(self.events_data[feature]) + 1} ")
        info_file.write(f"{len(self.geo_features)} ")
        nb_holidays = 0
        info_file.write(f"{nb_holidays}\n")
        events = self.events_data
        nb_obs_column_indexes = np.max(events[obs_index_column]) + 1
        obs_count = [0] * nb_obs_column_indexes
        for d in range(nb_obs_column_indexes):
            filtered_events = events.loc[events[obs_index_column] == d].copy()
            filtered_events["date_only"] = pd.to_datetime(filtered_events["ts"]).copy()
            date_only = filtered_events["date_only"].dt.date
            obs_count[d] = len(date_only.unique())
            info_file.write(f"{obs_count[d]} ")
        # print(f"obs_count = {obs_count}")
        # info_file.write(obs_count)
        info_file.write("\nEND")
        info_file.close()

    def write_arrivals(self, path="arrivals.dat"):
        """Write events data, aggregated by time indexes, geo index and features.

        Args:
            path (str, optional): path of the written file. Defaults to "arrivals.dat".
        """
        samples = self.get_events_aggregated()
        arrivals_file = open(path, "w")
        # arrivals_file.write("%s\n" % (" ".join([str(x) for x in samples.shape])))
        for index, sample in np.ndenumerate(samples):
            for i, val in enumerate(sample):
                line = "%s %d %d\n" % (" ".join([str(x) for x in index]), i, val)
                arrivals_file.write(line)
        arrivals_file.write("END")
        arrivals_file.close()

    def write_regions(self, path="neighbors.dat"):
        """Write file with description of discretization subregions.
        Description contains for each subregion:
            id,
            centroid coordinates,
            geo features,
            pairs of neighbor subregion ids and distances (km).

        Args:
            path (str, optional): path of the written file. Defaults to "neighbors.dat".
        """
        centers = self.geo_discretization.geometry.to_crs("epsg:29193").centroid.to_crs(
            self.geo_discretization.crs
        )
        coords = centers.apply(lambda x: x.representative_point().coords[:][0])
        neighbors_file = open(path, "w")
        # neighbors_file.write("%d %d\n" % (len(centers), len(self.geo_features)))
        # print(f"Writing geo_features {self.geo_features}")
        for i, row in self.geo_discretization.iterrows():
            id_region = row["id"]
            lat = coords[id_region][1]
            lon = coords[id_region][0]
            neighbors = row["neighbors"]
            features = [row[feature] for feature in self.geo_features]
            neighbors_file.write(
                "%d %.6f %.6f 0 %s "
                % (id_region, lat, lon, " ".join([str(x) for x in features]))
            )
            for neighbor in neighbors:
                row_neighbor = self.geo_discretization[
                    self.geo_discretization["id"] == neighbor
                ].iloc[0]
                id_neighbor = row_neighbor["id"]
                lat2, lon2 = coords[id_neighbor][1], coords[id_neighbor][0]
                d = distance((lat, lon), (lat2, lon2))
                neighbors_file.write(" %d %.3f" % (row_neighbor["id"], d))
            neighbors_file.write("\n")
        neighbors_file.write("END")
        neighbors_file.close()

    def plot_discretization(self, to_file=None):
        """Plot the current discretization. The plot contains: events as points, max_borders and subregion borders and ids.

        Args:
            to_file (string, optional): path indicating where to save the discretization figure. 

        Returns:
            _type_: _description_
        """
        import matplotlib.pyplot as plt

        self.geo_discretization["center"] = self.geo_discretization.geometry.to_crs(
            "epsg:29193"
        ).centroid.to_crs(self.geo_discretization.crs)
        self.geo_discretization["coords"] = self.geo_discretization["center"].apply(
            lambda x: x.representative_point().coords[:][0]
        )

        fig, ax = plt.subplots(figsize=(24, 16))
        self.geo_discretization.boundary.plot(ax=ax, aspect=1)
        self.events_data.dropna(subset=["gdiscr"]).plot(
            markersize=5, color="red", ax=ax
        )
        for _, row in self.geo_discretization.iterrows():
            plt.annotate(
                text=row["id"],
                xy=row["coords"],
                horizontalalignment="center",
                verticalalignment="center",
                color="black",
                fontsize=11,
            )
        if to_file != None:
            plt.savefig(to_file, bbox_inches="tight")
            return f"saved discretization plot at {to_file}"
        else:
            plt.show()
            return None
