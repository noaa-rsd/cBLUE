from pathlib import Path
from shapely.geometry import Polygon
from shapely.ops import transform, unary_union
import pyproj
from pyproj import CRS
from functools import partial
import subprocess
import json
import pdal
from osgeo import osr
import geopandas as gpd
from tqdm import tqdm


dir1 = Path(r'\\ngs-s-rsd\Lidar_Proc01\2016\FL1604-TB-N-880_MarcoIsland_p\06_RIEGL_PROC\04_EXPORT\Green\04_FL1604-TB-N-880_g_gpsa_rf_ip_wsf_r_adj_qc')
lasses = dir1.rglob('*.las')
bboxes = []
attributes = {
    'path': [],
    'hor_srs': [],
    }

wgs84_epsg = 'epsg:4326'


def run_console_cmd(cmd):
    process = subprocess.Popen(
        cmd.split(' '), shell=False, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.DEVNULL)
    output, error = process.communicate()
    returncode = process.poll()
    return returncode, output


for las_path in tqdm(list(lasses)[0:20]):
    print(las_path)

    las = str(las_path).replace('\\', '/')
    cmd_str = 'pdal info {} --metadata'.format(las)

    try:
        metadata = run_console_cmd(cmd_str)[1].decode('utf-8')
        meta_dict = json.loads(metadata)['metadata']

        srs = meta_dict['srs']
        hor_wkt = srs['horizontal']
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

        poly = transform(project, Polygon(tile_coords))
        bboxes.append(poly)

        path_parts = las_path.parts

        attributes['path'].append(str(las_path))
        attributes['hor_srs'].append(projcs)

    except Exception as e:
        print(e)

geom = [unary_union(bboxes)]
gdf = gpd.GeoDataFrame(geometry=bboxes, crs=wgs84_epsg)
print(gdf)
gdf['path'] = attributes['path']
gdf['hor_srs'] = attributes['hor_srs']
print(gdf)

gpkg = Path(r'z:\las_files.gpkg')
layer = '_2016'
gdf.to_file(gpkg, layer=layer, driver='GPKG')
