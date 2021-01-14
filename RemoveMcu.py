import json
import numpy as np
from pathlib import Path
import laspy


def get_mcu(json_path):
    with open(json_path) as J: 
        meta_data = json.load(J)
    return np.float(meta_data['VDatum region MCU']) / 100


if __name__ == '__main__':

    las_dir = Path(r'D:\RSD_PROJECTS\mcu_tpu_las')
    json_dir = las_dir

    for las in las_dir.glob('*.las'):
        json_path = json_dir / las.name.replace('.las', '.json')
        out_las_path = str(las).replace('.las', '_no_mcu.las')
        print('-' * 50)
        print(out_las_path)
        in_las = laspy.file.File(las, mode="r")
        out_las = laspy.file.File(out_las_path, mode="w", header=in_las.header)    

        for field in in_las.point_format:
            print(f'writing {field.name}...')
            las_data = in_las.reader.get_dimension(field.name)
            if field.name == 'total_tvu':
                mcu = get_mcu(json_path)
                print(f'removing MCU of {mcu}...')
                print(f'with MCU:  {las_data}')
                las_data = np.sqrt(las_data**2 - mcu**2)     
                print(f'without MCU:  {las_data}')
            out_las.writer.set_dimension(field.name, las_data)
