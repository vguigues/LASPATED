# spated

Spatio-Temporal Discretizator for spatio-Temporal data wrangling and analysis. If you have data from events with geographical and temporal information and want to create analytical datasets aggregating this information with other sources, `spated` can help you.

## Instalation

Just run `python -m pip install spated`.

## Usage

Use DataAggregator object to merge information from your multiple datasources.

```python
import spated

app = spated.DataAggregator()

# add your events dataset
# for example, a dataset with ambulance calls registers
app.add_events(ambulance_calls_df)

# add a base geographical limit as a geopandas.GeoDataFrame to locate the analysis
app.add_max_borders(base_map)
# you can also estimate a base map from our implemented methods
# app.add_max_borders(method='convex')

# compute time indexes from events
# for example, a 4-hour window daily sazonality
app.add_time_discretization(sazonality_type='H', window=4, frequency=24)

# compute geographical discretization
# for example, split your map into rectangles with 15 horizontal splits and 20 vertical splits
app.add_geo_discretization(discr_type='R', rect_discr_param_x=15, rect_discr_param_y=20)

# and finally, add more information from other geo-located data sources
# for example, the population from neighborhoods of original map
app.add_geo_features(neighborhoods_population)
```

For more details and options, check `usage_example.ipynb` notebook.
