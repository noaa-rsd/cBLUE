from pathlib import Path
from shapely.geometry import Polygon
from shapely.ops import transform
import pyproj
from pyproj import CRS
from functools import partial
import subprocess
import json
from osgeo import osr
import geopandas as gpd
from tqdm import tqdm


def run_console_cmd(cmd):
    process = subprocess.Popen(
        cmd.split(' '), shell=False, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.DEVNULL)
    output, error = process.communicate()
    returncode = process.poll()
    return returncode, output


bboxes = []
attributes = {
    'name': [],
    'parent': [],
    'hor_srs': [],
    'las_version': [],
    }

wgs84_epsg = 'epsg:4326'
gpkg = Path(r'z:\las_files.gpkg')

las_dir = Path(r'T:\2017')
lasses = list(las_dir.rglob('*.las'))

for las_path in tqdm(lasses):

    las = str(las_path).replace('\\', '/')
    cmd_str = 'pdal info {} --metadata'.format(las)

    try:
        metadata = run_console_cmd(cmd_str)[1].decode('utf-8')
        meta_dict = json.loads(metadata)['metadata']

        major_version = meta_dict['major_version']
        minor_version = meta_dict['minor_version']
        las_version = f'{major_version}.{minor_version}'

        hor_wkt = meta_dict['srs']['horizontal']
        hor_srs = osr.SpatialReference(wkt=hor_wkt) 
        projcs = hor_srs.GetAttrValue('projcs')

        minx = meta_dict['minx']
        miny = meta_dict['miny']
        maxx = meta_dict['maxx']
        maxy = meta_dict['maxy']

        tile_coords = [
                (minx, miny),
                (minx, maxy),
                (maxx, maxy),
                (maxx, miny)
            ]

        project = partial(
            pyproj.transform,
            pyproj.Proj(CRS.from_string(projcs)),
            pyproj.Proj(wgs84_epsg))

        bboxes.append(transform(project, Polygon(tile_coords)))

        attributes['name'].append(las_path.name)
        attributes['parent'].append(str(las_path.parent))
        attributes['hor_srs'].append(projcs)
        attributes['las_version'].append(las_version)

    except Exception as e:
        print(f'{e} - {str(las_path)}')

print('creating GeoDataFrame...')
gdf = gpd.GeoDataFrame(geometry=bboxes, crs=wgs84_epsg)
gdf.geometry = gdf.geometry.map(lambda poly: transform(lambda x, y: (y, x), poly))

gdf['name'] = attributes['name']
gdf['parent'] = attributes['parent']
gdf['hor_srs'] = attributes['hor_srs']
gdf['las_version'] = attributes['las_version']

print('writing geopackage...')
layer = f'Lidar_Proc01__2017'
gdf.to_file(gpkg, layer=layer, driver='GPKG')
