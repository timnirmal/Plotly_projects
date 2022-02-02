import numpy as np
from netCDF4 import Dataset
from numpy import pi, sin, cos


# Import latitde, longitude and topography data from netCDF file
def Etopo(lat_area, lon_area, resolution):
    """":cvar
    lat_area: latitude of the area
    lon_area: longitude of the area
    Combined lat and lan: Region of the map you want [90,110]
    resolution: resolution of the topography for both latitude and longitude (deg) - (Original resolution is 0.0167 deg)

    Output:
    Mesh type, Longitude and Latitude of the topography
    """
    # Open netCDF file
    data = Dataset("ETOPO1_Bed_g_gdal.grd", "r")

    # Read Data
    lat_range = data.variables['y_range'][:]
    lon_range = data.variables['x_range'][:]
    topo = data.variables['z_range'][:]
    spacing = data.variables['spacing'][:]
    dimension = data.variables['dimension'][:]
    z = data.variables['z'][:]
    lon_num = dimension[0]
    lat_num = dimension[1]

    # Close netCDF file
    data.close()

    # Prepare array
    lon_input = np.zeros(lon_num)
    lat_input = np.zeros(lat_num)

    for i in range(lon_num):
        lon_input[i] = lon_range[0] + i * spacing[0]
    for i in range(lat_num):
        lat_input[i] = lat_range[0] + i * spacing[1]

    # 2D mesh-grid (array)
    lon, lat = np.meshgrid(lon_input, lat_input)

    # 2D to 1D for z
    topo = np.reshape(z, (lat_num, lon_num))

    # Skip the data for resolution
    if ((resolution < spacing[0]) | (resolution < spacing[1])):
        print('Set the highest resolution')
    else:
        skip = int(resolution / spacing[0])
        lon = lon[::skip, ::skip]
        lat = lat[::skip, ::skip]
        topo = topo[::skip, ::skip]

    topo = topo[::-1]

    # Select the range of map
    range1 = np.where((lon >= lon_area[0]) & (lon <= lon_area[1]))
    lon = lon[range1]
    lat = lat[range1]
    topo = topo[range1]
    range2 = np.where((lat >= lat_area[0]) & (lat <= lat_area[1]))
    lon = lon[range2]
    lat = lat[range2]
    topo = topo[range2]

    # Convert 2D again
    lon_num = len(np.unique(lon))
    lat_num = len(np.unique(lat))
    lon = np.reshape(lon, (lat_num, lon_num))
    lat = np.reshape(lat, (lat_num, lon_num))
    topo = np.reshape(topo, (lat_num, lon_num))

    return lon, lat, topo


# Functions from https://chart-studio.plotly.com/~empet/14813/heatmap-plot-on-a-spherical-map/#/

def degree2radians(degree):
    # convert degrees to radians
    return degree * pi / 180


def mapping_map_to_sphere(lon, lat, radius=1):
    # this function maps the points of coords (lon, lat) to points onto the  sphere of radius radius
    lon = np.array(lon, dtype=np.float64)
    lat = np.array(lat, dtype=np.float64)
    lon = degree2radians(lon)
    lat = degree2radians(lat)
    xs = radius * cos(lon) * cos(lat)
    ys = radius * sin(lon) * cos(lat)
    zs = radius * sin(lat)
    return xs, ys, zs


# Selecting whole globe with small resolution (or there will be lots of data)
lon_topo, lat_topo, topo = Etopo([-90, 90], [-180, 180], 0.8)

# Convert to spherical coordinates
xs, ys, zs = mapping_map_to_sphere(lon_topo, lat_topo)


topo_sphere=dict(type='surface',
  x=xs,
  y=ys,
  z=zs,
  surfacecolor=topo)

import plotly.graph_objs as go

titlecolor = 'white'
bgcolor = 'black'

layout = go.Layout(
  autosize=False, width=1200, height=800,
  title = '3D spherical topography map',
  titlefont = dict(family='Courier New', color=titlecolor),
  showlegend = False,
  scene = dict(
    aspectmode='manual',
    aspectratio=go.layout.scene.Aspectratio(
      x=1, y=1, z=1)),
  paper_bgcolor = bgcolor,
  plot_bgcolor = bgcolor)

plot_data=[topo_sphere]
fig = go.Figure(data=plot_data, layout=layout)
fig.show()





