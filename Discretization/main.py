import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
import numpy as np

import laspated as spated

# def read_emergencies():
#     df =  pd.read_csv("../Data/emergency_calls_rio_de_janeiro_old.csv")
#     df["data_hora"] = df["data"] + " " + df["hora"]
#     print(df)
#     df.to_csv("../Data/emergency_calls_rio_de_janeiro.csv", index=False)


# def move_events_to_ny(max_borders, events):
#     num_events_lat = events["lat"].count()
#     num_events_long = events["long"].count()
#     if num_events_lat != num_events_long:
#         print("Weird lat and long count is different", num_events_lat, num_events_long)
#         input()
#     print(num_events_lat, num_events_long)
#     print(max_borders.geometry)
#     new_coords = max_borders.geometry.sample_points(num_events_lat).explode(index_parts=True)
#     coord_list = [(x,y) for x,y in zip(new_coords.x , new_coords.y)]
#     print(coord_list[:10])
#     events["long"] = pd.Series([x[0] for x in coord_list])
#     events["lat"] = pd.Series([x[1] for x in coord_list])
    
#     events.to_csv("new_york_emergencies.csv",index=False)
#     input("File saved")


def distance(p1, p2):
    pass

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
    info_file.close()

    print(app.geo_discretization.head(10))
    centers = app.geo_discretization.geometry.to_crs("epsg:29193").centroid.to_crs(app.geo_discretization.crs)
    coords = centers.apply(lambda x: x.representative_point().coords[:][0])
    
    neighbors_file = open("neighbors.dat", "w")
    for row in app.geo_discretization.iterrows():
        land_type = 1
        neighbors_file.write("%d %.6f %.6f %d " % (row["id"], centers["lat"], centers["long"], land_type))
        

    input()


def example_ny():
    app = spated.DataAggregator(crs="epsg:4326")
    max_borders = gpd.read_file(r'../Data/ny/')
    app.add_max_borders(max_borders)
    app.max_borders.plot()
    plt.savefig("ny.png")
    events = pd.read_csv(r'../Data/emergency_calls_new_york.csv', encoding = "ISO-8859-1", sep=",")
    app.add_events_data(events, datetime_col='data_hora', lat_col='lat', lon_col="long", feature_cols=['prioridade'])

    app.add_time_discretization('D', 1, 7, column_name="dow")
    app.add_time_discretization('m', 30, 60*24, column_name="hhs")

    app.add_geo_discretization(
        discr_type='R',
        rect_discr_param_x=10,
        rect_discr_param_y=10
    )

    centers = app.geo_discretization.geometry.to_crs("epsg:29193").centroid.to_crs(app.geo_discretization.crs)
    coords = centers.apply(lambda x: x.representative_point().coords[:][0])

def example_rj():
    app = spated.DataAggregator(crs="epsg:4326")
    max_borders = gpd.read_file(r'../Data/rj/')
    app.add_max_borders(max_borders)
    app.max_borders.plot()
    plt.savefig("rj.png")
    events = pd.read_csv(r'../Data/emergency_calls_rio_de_janeiro.csv', encoding = "ISO-8859-1", sep=",")
    events = events.drop(events[events["prioridade"] > 2].index)
    app.add_events_data(events, datetime_col='data_hora', lat_col='lat', lon_col="long", feature_cols=['prioridade'])


    app.add_time_discretization('D', 1, 7, column_name="dow")
    app.add_time_discretization('m', 30, 60*24, column_name="hhs")

    app.add_geo_discretization(
        discr_type='R',
        rect_discr_param_x=10,
        rect_discr_param_y=10
    )

    centers = app.geo_discretization.geometry.to_crs("epsg:29193").centroid.to_crs(app.geo_discretization.crs)
    coords = centers.apply(lambda x: x.representative_point().coords[:][0])

    population = gpd.read_file(r'../Data/regressores/populacao/')
    population = population[['populacao_','populaca_1','populaca_2','populaca_3','geometry']].copy()
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
        print(land_use.groupby("sub_group_%d" % (i))["id"].apply(list))
        land_use["sub_group_%d" % (i)] = np.where(land_use["sub_group_%d" % (i)] in sub_group, 1, 0)

    print(land_use.sample(10))
    # land_use["grupo_id"] = np.nan
    # land_use["grupo_id"] = land_use["grupo"]
    # land_use["grupo_id"] = np.where(land_use["grupo_id"] == "Áreas urbanizadas", 1, 0)
    # land_use = land_use[['grupo_id','geometry']].copy()
    app.add_geo_variable(land_use)
    print(app.geo_discretization.sample(10))
    write_files(app)

def main():
    example_rj()

if __name__ == "__main__":
    main()

