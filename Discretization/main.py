import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
import numpy as np
import geodatasets as gds

import laspated as spated

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
    num_feature_types = 3
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

def example_rj():
    app = spated.DataAggregator(crs="epsg:4326") # initializes data aggregator
    max_borders = gpd.read_file(r'../Data/rj/rj.shp') # Load the geometry of region of interest
    # app.add_max_borders(max_borders) # Add the border of region 
    events = pd.read_csv(r'../Data/emergency_calls_rio_de_janeiro.csv', encoding = "ISO-8859-1", sep=",")
    events = events.drop(events[events["prioridade"] > 2].index)
    app.add_events_data(events.sample(20), datetime_col='data_hora', lat_col='lat', lon_col="long", feature_cols=['prioridade']) # %m/%d/%y %H:%M:%S
    # print(app.events_data.sample(20))
    # app.add_max_borders(method="convex")
    # # # app.max_borders.plot()
    # fig, ax = plt.subplots()
    # app.max_borders.plot(ax=ax)
    # app.events_data.plot(markersize=10, color='red', ax=ax)
    # # plt.show()
    # plt.savefig("convex_rj.png")


    app.add_time_discretization('D', 1, 7, column_name="dow")
    app.add_time_discretization('m', 30, 60*24, column_name="hhs")
    

    app.add_geo_discretization(
        discr_type='R',
        rect_discr_param_x=10,
        rect_discr_param_y=10
    )

    # print(app.events_data.head(30))
    # input()

    # centers = app.geo_discretization.geometry.to_crs("epsg:29193").centroid.to_crs(app.geo_discretization.crs)
    # coords = centers.apply(lambda x: x.representative_point().coords[:][0])

    population = gpd.read_file(r'../Data/regressores/populacao/')
    population = population[['populacao_','geometry']].copy()
    app.add_geo_variable(population)
    land_use = gpd.read_file(r'../Data/regressores/uso_do_solo/')
    print(land_use.columns)
    print(np.unique(land_use["usoagregad"]))
    groups = [x for x in pd.unique(land_use["grupo"])]
    sub_groups = [['Afloramentos rochosos e depósitos sedimentares', 'Cobertura arbórea e arbustiva', 
                'Cobertura gramíneo lenhosa', 'Corpos hídricos', 'Áreas de exploração mineral',
                'Áreas não edificadas'], ['Favela', 'Áreas sujeitas à inundação', 'Áreas agrícolas'], 
                ['Áreas de educação e saúde', 'Áreas de lazer', 'Áreas de transporte', 'Áreas institucionais e de infraestrutura pública'],
                ['Áreas residenciais']]
    for i,sub_group in enumerate(sub_groups):
        land_use["sub_group_%d" % (i)] = land_use["usoagregad"]
    
    print(land_use.columns)
    for i,row in land_use.iterrows():
        uso = row["usoagregad"]
        for j,sub_group in enumerate(sub_groups):
            row["sub_group_%d" % (j)] = 1 if uso in sub_group else 0
        land_use.iloc[i] = row

    print(land_use.sample(10))
    # land_use["grupo_id"] = np.nan
    # land_use["grupo_id"] = land_use["grupo"]
    # land_use["grupo_id"] = np.where(land_use["grupo_id"] == "Áreas urbanizadas", 1, 0)
    land_use = land_use[['sub_group_0', 'sub_group_1', 'sub_group_2', 'sub_group_3','geometry']].copy()
    app.add_geo_variable(land_use)
    print(app.geo_discretization.sample(10))
    write_files(app)

def main():
    # read_calls()
    # example_rj()
    example_ny()

if __name__ == "__main__":
    main()

