import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
import numpy as np

import laspated as spated
from shapely.geometry import Polygon, MultiPolygon, Point

# def read_emergencies():
#     df =  pd.read_csv("../Data/emergency_calls_rio_de_janeiro_old.csv")
#     df["data_hora"] = df["data"] + " " + df["hora"]
#     print(df)
#     df.to_csv("../Data/emergency_calls_rio_de_janeiro.csv", index=False)


def move_events_to_ny(max_borders, events):
    num_events_lat = events["lat"].count()
    num_events_long = events["long"].count()
    if num_events_lat != num_events_long:
        print("Weird lat and long count is different", num_events_lat, num_events_long)
        input()
    print(num_events_lat, num_events_long)
    new_coords = max_borders.geometry.sample_points(num_events_lat).explode(index_parts=True).to_crs("epsg:4326")
    coord_list = [(x,y) for x,y in zip(new_coords.x , new_coords.y)]
    print(coord_list[:10])
    events["long"] = pd.Series([x[0] for x in coord_list])
    events["lat"] = pd.Series([x[1] for x in coord_list])
    
    events.to_csv("new_york_emergencies.csv",index=False)
    input("File saved")


def read_calls():
    T = 48
    G = 7
    R = 76
    P = 3
    calls_shape = (T,G,R,P)
    nb_observations = np.zeros(calls_shape)
    nb_calls = np.zeros(calls_shape)
    sample = np.empty(shape=calls_shape, dtype=list)
    for t in range(T):
        for g in range(G):
            for r in range(R):
                for p in range(P):
                    sample[t,g,r,p] = []

    calls_file = open("calls.dat", "r")
    for line in calls_file.readlines():
        if line == "END":
            break
        t,g,r,p,j,val,h = [int(x) for x in line.split()]
        nb_observations[t,g,r,p] += 1
        nb_calls[t,g,r,p] += val
        if sample[t,g,r,p] is list:
            sample[t,g,r,p].append(val)
        else:
            sample[t,g,r,p] = [val]

    calls_file.close()
    

    for t in range(T):
        for g in range(G):
            for r in range(R):
                for p in range(P):
                    if len(sample[t,g,r,p]) > 0:
                        print(t,g,r,p, sample[t,g,r,p])
                        input()


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

def write_files(app):
    num_time_windows = 48
    num_days = 7
    print(pd.unique(app.events_data["prioridade"]), app.events_data.count())
    num_feature_types = 4
    print(app.geo_discretization.sample(10))
    num_regions = np.max(app.geo_discretization["id"]) + 1
    num_regressors = 1
    num_holidays = 0
    info_file = open("info.dat", "w")
    info_file.write("%d %d %d %d %d %d\n" % (num_time_windows, num_days, num_regions, num_feature_types, num_regressors, num_holidays))
    for i in range(7):
        info_file.write("%d " % (len(app.events_data[ app.events_data["dow"] == i])))
    info_file.write("END")
    info_file.close()

    centers = app.geo_discretization.geometry.to_crs("epsg:29193").centroid.to_crs(app.geo_discretization.crs)
    coords = centers.apply(lambda x: x.representative_point().coords[:][0])
    print(coords)
    neighbors_file = open("neighbors.dat", "w")
    for i,row in app.geo_discretization.iterrows():
        land_type = 1
        id_region = row["id"]
        lat = coords[id_region][1]
        lon = coords[id_region][0]
        neighbors = row["neighbors"]
        regs = [row["populacao_"], row["sub_group_0"], row["sub_group_1"], row["sub_group_2"], row["sub_group_3"]]
        neighbors_file.write("%d %.6f %.6f 0 " % (id_region, lat, lon))
        for reg in regs:
            neighbors_file.write("%.4f " % reg)
        for neighbor in neighbors:
            row_neighbor = app.geo_discretization[app.geo_discretization["id"] == neighbor].iloc[0]
            id_neighbor = row_neighbor["id"]
            lat2,lon2 = coords[id_neighbor][1], coords[id_neighbor][0]
            d = distance((lat,lon), (lat2,lon2))
            neighbors_file.write("%d %.3f " % (row_neighbor["id"], d))
        neighbors_file.write("\n")
    neighbors_file.write("END")
    neighbors_file.close()

    calls_file = open("arrivals.dat", "w")
    print(app.events_data.sample(10))
    nb_observations = np.zeros((num_time_windows, num_days, num_regions, num_feature_types))

    for i,row in app.events_data.iterrows():
        try:
            t = int(row["hhs"])
            g = int(row["dow"])
            r = int(row["gdiscr"])
            p = int(row["prioridade"])
        except ValueError:
            continue
        
        nb_obs = nb_observations[t,g,r,p]
        calls_file.write("%d %d %d %d %d 1 -1\n" % 
            (t, g, r, p, nb_obs))
        nb_observations[t,g,r,p] += 1
    calls_file.write("END")
    calls_file.close()


def write_files_ny(app):
    num_time_windows = 48
    num_days = 7
    num_feature_types = 3
    num_regions = np.max(app.geo_discretization["id"]) + 1
    num_regressors = 0
    num_holidays = 0
    info_file = open("info_ny.dat", "w")
    info_file.write("%d %d %d %d %d %d\n" % (num_time_windows, num_days, num_regions, num_feature_types, num_regressors, num_holidays))
    for i in range(7):
        info_file.write("%d " % (len(app.events_data[ app.events_data["dow"] == i])))
    info_file.write("END")
    info_file.close()

    centers = app.geo_discretization.geometry.to_crs("epsg:29193").centroid.to_crs(app.geo_discretization.crs)
    coords = centers.apply(lambda x: x.representative_point().coords[:][0])
    neighbors_file = open("neighbors_ny.dat", "w")
    for i,row in app.geo_discretization.iterrows():
        land_type = 1
        id_region = row["id"]
        lat = coords[id_region][1]
        lon = coords[id_region][0]
        neighbors = row["neighbors"]
        regs = [0,0,0,0,0]
        neighbors_file.write("%d %.6f %.6f 0 " % (id_region, lat, lon))
        for reg in regs:
            neighbors_file.write("%.4f " % reg)
        for neighbor in neighbors:
            row_neighbor = app.geo_discretization[app.geo_discretization["id"] == neighbor].iloc[0]
            id_neighbor = row_neighbor["id"]
            lat2,lon2 = coords[id_neighbor][1], coords[id_neighbor][0]
            d = distance((lat,lon), (lat2,lon2))
            neighbors_file.write("%d %.3f " % (row_neighbor["id"], d))
        neighbors_file.write("\n")
    neighbors_file.write("END")
    neighbors_file.close()

    calls_file = open("calls_ny.dat", "w")
    nb_observations = np.zeros((num_time_windows, num_days, num_regions, num_feature_types))
    for i,row in app.events_data.iterrows():
        try:
            t = int(row["hhs"])
            g = int(row["dow"])
            r = int(row["gdiscr"])
            p = int(row["prioridade"])
        except ValueError:
            continue
        nb_obs = nb_observations[t,g,r,p]
        calls_file.write("%d %d %d %d %d 1 -1\n" % 
            (t, g, r, p, nb_obs))
        nb_observations[t,g,r,p] += 1
    calls_file.write("END")
    calls_file.close()

def example_ny():
    app = spated.DataAggregator(crs="epsg:4326")
    max_borders = gpd.read_file(r'../Data/ny/')
    app.add_max_borders(max_borders)
    app.max_borders.plot()
    plt.savefig("ny.png")
    events = pd.read_csv(r'../Data/emergency_calls_new_york.csv', encoding = "ISO-8859-1", sep=",")
    events = events.drop(events[events["prioridade"] > 2].index)
    app.add_events_data(events, datetime_col='data_hora', lat_col='lat', lon_col="long", feature_cols=['prioridade'])

    app.add_time_discretization('D', 1, 7, column_name="dow")
    app.add_time_discretization('m', 30, 60*24, column_name="hhs")

    app.add_geo_discretization(
        discr_type='R',
        rect_discr_param_x=10,
        rect_discr_param_y=10
    )

    write_files_ny(app)


def random_points_in_shp(shp, n):
    within = False
    samples = []
    while len(samples) < n:
        x = np.random.uniform(shp.bounds[0], shp.bounds[2])
        y = np.random.uniform(shp.bounds[1], shp.bounds[3])
        within = shp.contains(Point(x, y))
        if within:
            samples.append((x,y))
    
    return samples

def example_rj():
    app = spated.DataAggregator(crs="epsg:4326") # initializes data aggregator
    max_borders = gpd.read_file(r'../Data/rj/rj.shp') # Load the geometry of region of interest
    app.add_max_borders(max_borders) # Add the border of region 
    events = pd.read_csv(r'sorted_events.csv', encoding = "ISO-8859-1", sep=",")
    # [46737, 102558, 41309, 85621]
    app.add_events_data(events, datetime_col='data_hora', lat_col='lat', lon_col="long", feature_cols=['prioridade'],
                        datetime_format="%m/%d/%y %H:%M:%S") # %m/%d/%y %H:%M:%S

    app.add_time_discretization('D', 1, 7, column_name="dow")
    app.add_time_discretization('m', 30, 60*24, column_name="hhs")

    bases = gpd.read_file("bases/bases.shp")
    app.add_geo_discretization('V', custom_data=bases)

    # app.plot_discretization()

    # population = gpd.read_file(r'../Data/regressores/populacao/')
    # population = population[['populacao_','geometry']].copy()
    # app.add_geo_variable(population)

    # land_use = gpd.read_file(r'../Data/regressores/uso_do_solo/')
    # land_use = land_use[['subgroup_0', 'subgroup_1', 'subgroup_2', 'subgroup_3','geometry']].copy()
    # app.add_geo_variable(land_use,type_geo_variable="area")
    # print(app.geo_discretization[["id", "subgroup_0","subgroup_1","subgroup_2","subgroup_3", "populacao_"]])
    # print(app.geo_discretization.isna().any())

    # A = app.get_events_aggregated()
    # print(A.shape)
    # app.write_arrivals("test.dat")
    # app.write_regions("testr.dat")

    print(app.geo_discretization)
    samples = []
    arq_sample = open(r"../Data/Bases_voronoi/samples.dat", "w")
    arq_csv = open(r"samples.csv", "w")
    arq_csv.write("region_id,sample_id,lat,long\n")
    for i,row in app.geo_discretization.iterrows():
        sample_i = random_points_in_shp(row["geometry"], 100)
        for j,coords in enumerate(sample_i):
            lon,lat = coords
            arq_sample.write(f"{i} {j} {lat} {lon}\n")
            arq_csv.write(f"{i},{j},{lat},{lon}\n")
    arq_sample.close()
    arq_csv.close()

def main():
    # read_calls()
    example_rj()
    # example_ny()

if __name__ == "__main__":
    main()



    # centers = app.geo_discretization.geometry.to_crs("epsg:29193").centroid.to_crs(app.geo_discretization.crs)
    # coords = centers.apply(lambda x: x.representative_point().coords[:][0])

    # fig, ax = plt.subplots()
    # app.geo_discretization.boundary.plot(ax=ax,aspect=1)
    # app.events_data.head(10).dropna(subset=['gdiscr']).plot(markersize=15, color='red', ax=ax)
    # for ind, row in app.geo_discretization.iterrows():
    #     # print(coords.loc[ind])
    #     plt.annotate(row["id"],coords.loc[ind], 
    #                  horizontalalignment="center", 
    #                  verticalalignment="center", 
    #                  color='black', 
    #                  fontsize=9)

    # increments_dic = {0: (0.015,0), 1: (0.01,0.025), 2: (0,0.014), 
    #                   3: (0.015,0.017), 4: (0.01,-0.02), 5: (0.01,0),
    #                   6: (-0.01,-0.02), 7: (-0.02,-0.02), 8: (-0.01,-0.02),
    #                   9: (-0.01,0.02)}


    # for ind,row in app.events_data.iterrows():
    #     if ind >= 10:
    #         break
    #     event_row = events.loc[ind]
    #     # print((event_row["long"], event_row["lat"]))
    #     inc_x, inc_y = increments_dic[ind]
    #     val = row["hhs"] if ind != 8 else 6
    #     plt.annotate(
    #         # ind, row["hhs"]
    #         val, (event_row["long"] + inc_x, event_row["lat"] + inc_y),
    #         horizontalalignment="center",
    #         verticalalignment="center",
    #         color='green', fontsize=13
    #     )
    # plt.savefig("Video/first_10_calls_with_discretization.png")
    # plt.show()
    # print("Saved new fig")
    # quit()

    # print(app.events_data)
    # input()