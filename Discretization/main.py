import pandas as pd
import geopandas as gpd

import matplotlib.pyplot as plt
import numpy as np

import laspated as spated


# def remove_empty(df):
#     df["lat"].replace('', np.nan, inplace=True)
#     df["long"].replace('', np.nan, inplace=True)
#     df.dropna(subset=["lat", "long"], inplace=True)


def main():
    max_borders = gpd.read_file(r'../Data/rj/')
    app = spated.DataAggregator(crs="epsg:4326")
    app.add_max_borders(max_borders)
    # app.max_borders.plot()
    # plt.show()
    events = pd.read_csv(r'../Data/emergency_calls_rio_de_janeiro.csv', encoding = "ISO-8859-1", sep=",")
    app.add_events_data(events.sample(1000), datetime_col='data', lat_col='lat', lon_col="long", feature_cols=['prioridade'],
        datetime_format="%m/%d/%y")
    
    print(events.columns)
    app.add_time_discretization('D', 1, 7, column_name="dow")
    print(app.events_data.sample(20))



if __name__ == "__main__":
    main()

