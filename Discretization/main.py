import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
import numpy as np

import laspated as spated


def read_calls(calls_path):
    arq = open(calls_path, "r")

    arq.close()


# def read_emergencies():
#     df =  pd.read_csv("../Data/emergency_calls_rio_de_janeiro_old.csv")
#     df["data_hora"] = df["data"] + " " + df["hora"]
#     print(df)
#     df.to_csv("../Data/emergency_calls_rio_de_janeiro.csv", index=False)



def main():
    max_borders = gpd.read_file(r'../Data/rj/')
    app = spated.DataAggregator(crs="epsg:4326")
    app.add_max_borders(max_borders)
    # app.max_borders.plot()
    # plt.show()
    events = pd.read_csv(r'../Data/emergency_calls_rio_de_janeiro.csv', encoding = "ISO-8859-1", sep=",")
    # print(events.head(10))
    app.add_events_data(events, datetime_col='data_hora', lat_col='lat', lon_col="long", feature_cols=['prioridade'])
    app.add_time_discretization('D', 1, 7, column_name="dow")
    app.add_time_discretization('m', 30, 60*24, column_name="hhs")
    print(app.events_data.head(100))

    app.add_geo_discretization(
        discr_type='R',
        rect_discr_param_x=10,
        rect_discr_param_y=10
    )
    centers = app.geo_discretization.geometry.to_crs("epsg:29193").centroid.to_crs(app.geo_discretization.crs)
    coords = centers.apply(lambda x: x.representative_point().coords[:][0])

    population = gpd.read_file(r'../Data/regressores/populacao/')
    # print(population.columns)
    # print(population.head(10))
    population = population[['populacao_','populaca_1','populaca_2','populaca_3','geometry']].copy()
    app.add_geo_variable(population)
    # print(app.events_data.head(50))
    land_use = gpd.read_file(r'../Data/regressores/uso_do_solo/')
    print(land_use.columns)
    print(land_use.sample(100))
    print(land_use["ruleid"].unique())
    print(land_use["usoagregad"].unique())
    print(land_use["grupo"].unique())
    land_use = land_use[['shape_Area','geometry']].copy()
    app.add_geo_variable(land_use)
    print(app.geo_discretization)
    

if __name__ == "__main__":
    main()

