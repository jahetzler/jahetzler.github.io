# Fiskeridirektoratet fangstdata

[Link to example map](arc_layers.html) (file with all points is too big (200mb) to upload to github).

I have used pydeck to plot in which municipality fishing boats have unloaded their catch and from which area the area the fish have been caught in.

datasets used in this is the [fangstdata](https://www.fiskeridir.no/Tall-og-analyse/AApne-data/Fangstdata-seddel-koblet-med-fartoeydata) from fiskeridirektoratet, and [robhop](https://github.com/robhop/) municipality geojson data [link](https://github.com/robhop/fylker-og-kommuner/blob/main/Kommuner-L.geojson).

## Libraries 

```
import pydeck as pdk
import pandas as pd
import geopandas as gpd
import json

#Libraries
import ast
import pandas as pd
import geopandas as gpd
import numpy as np
import json
import requests
import folium
import pydeck as pdk

import plotly.express as px
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib import animation
from matplotlib.widgets import CheckButtons, Slider
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

```
## Import and format data from fiskeridirektoratet
```
df = pd.read_csv('fangstdata_2024.csv.zip', sep=';')

dtypes = {
    "Landingskommune": "category",
    "Lengdegruppe": "category",
    "Art FAO": "category",
    "Anvendelse": "category"}

df = df.astype(dtypes)

df['Bruttovekt'] = df['Bruttovekt'].str.replace(',','.')
df['Bruttovekt'] = pd.to_numeric(df['Bruttovekt'])
df['Fangstår'] = pd.to_numeric(df['Fangstår'])
df['Landingskommune'] = df['Landingskommune'].str.capitalize()

df["Produksjonskommune (kode)"] = df["Produksjonskommune (kode)"].astype('str')
df['Produksjonskommune (kode)'] = df['Produksjonskommune (kode)'].str.replace('.0','')

df["Lon (hovedområde)"] = df["Lon (hovedområde)"].astype('str')
df["Lon (hovedområde)"] = df["Lon (hovedområde)"].str.replace(',','.')
df["Lat (hovedområde)"] = df["Lat (hovedområde)"].astype('str')
df["Lat (hovedområde)"] = df["Lat (hovedområde)"].str.replace(',','.')

df["Lon (lokasjon)"] = df["Lon (lokasjon)"].astype('str')
df["Lon (lokasjon)"] = df["Lon (lokasjon)"].str.replace(',','.')
df["Lat (lokasjon)"] = df["Lat (lokasjon)"].astype('str')
df["Lat (lokasjon)"] = df["Lat (lokasjon)"].str.replace(',','.')

```
## Import municipality map 

```

with open('Kommuner-L.geojson', encoding='utf-8') as fh:
    GEO = json.load(fh)
# GEO["features"][0]


gdf = gpd.GeoDataFrame.from_features(GEO,  crs='epsg:4326')

gdf["centroid"] = gdf.centroid
gdf["x"] = gdf.centroid.x
gdf["y"] = gdf.centroid.y

```
## Merge datasets

```
df = pd.merge(df, gdf[["id", "x", "y"]], how="right", left_on="Produksjonskommune (kode)", right_on = "id")
```

## Match up landing point with centroid for the municipality from the municipality data
```
df_pydeck = df[["Lon (lokasjon)", "Lat (lokasjon)", "x", "y", "Landingskommune", "Hovedområde"]]


df_pydeck["lat"] = df_pydeck["Lat (lokasjon)"]
df_pydeck["lon"] = df_pydeck["Lon (lokasjon)"]

df_pydeck = df_pydeck[["lon", "lat", "x", "y", "Landingskommune", "Hovedområde"]]

df_pydeck['x'] = df_pydeck['x'].apply(lambda x: round(x, 5))
df_pydeck['y'] = df_pydeck['y'].apply(lambda x: round(x, 5))


df_pydeck["lon"] = df_pydeck["lon"].astype('float64')
df_pydeck["lat"] = df_pydeck["lat"].astype('float64')
df_pydeck.dtypes
# df_pydeck.head(200)

df_pydeck["vekt"] = df["Bruttovekt"]

df_pydeck = df_pydeck.dropna()

```

## plot in pydeck

This code only runs a subset of the dataset (Bodø and Ørland) to upload an example map, the full html file for the dataset is around 200mb. 
```
# Specify a deck.gl ArcLayer
Bodoe = pdk.Layer(
    "ArcLayer",
    # "GreatCircleLayer",
    data=df_pydeck[df_pydeck["Landingskommune"] == "Bodø"],
    get_width="S000 * 2",
    get_source_position=["lon", "lat"],
    get_target_position=["x", "y"],
    # get_tilt=15,
    get_source_color=[0, 255, 0, 40],
    get_target_color=[0, 255, 0, 80],
    pickable=True,
    auto_highlight=True,
)

Ørland = pdk.Layer(
    "ArcLayer",
    # "GreatCircleLayer",
    data=df_pydeck[df_pydeck["Landingskommune"] == "Ørland"],
    get_width="S000 * 2",
    get_source_position=["lon", "lat"],
    get_target_position=["x", "y"],
    # get_tilt=15,
    get_source_color=[240, 100, 0, 40],
    get_target_color=[240, 100, 0, 80],
    pickable=True,
    auto_highlight=True,
)

view_state = pdk.ViewState(
    latitude=64,
    longitude=12,
    zoom=6,
    min_zoom=2,
    max_zoom=18,
    pitch=40.5,
    bearing=-27.36
)

TOOLTIP_TEXT = {'html': '<b>Brutto Vekt:</b> {vekt} <br> <b>Kommune:</b> {Landingskommune}',
        'style': {
            'color': 'white'}
}

r = pdk.Deck(initial_view_state=view_state, layers = [Bodoe, Ørland],  tooltip=TOOLTIP_TEXT)
r.to_html("arc_layers.html")

```
