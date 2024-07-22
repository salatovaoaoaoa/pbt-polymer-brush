import sys
import os
import subprocess
sys.path.append('/home/tpopova/prj/polymer_brush')

from math import sqrt
from math import exp
from math import pi

def generate_pore_in_files(template_pore: str = 'pore_temp.in',
                   
                   #название параметра
                   range_param: str = 'Cs',

                   N_brush: int = 200,
                   S: int = 100,
                   Cs: float = 0.001,
                   alpha: float = -0.5,
                   D: float = 400,
                   min_range_value: float = 4,
                   max_range_value: float = 4,
                   ): 
    #какие параметры меняем?
    change_param = Cs
    
    #theta
    
    l_t = S/(2*pi*D)
    theta = N_brush/l_t

    with open(template_pore, 'r') as file:
        data = file.readlines()

    # Изменяем указанные параметры
    for str in range(len(data)):
        if 'lat : 1G : n_layers' in data[str]:
            data[str] = f'lat : 1G : n_layers : {D}\n'
        elif 'mon : X0 : valence' in data[str]:
            data[str] = f'mon : X0 : valence : {-1 * alpha}\n'
        elif 'mon : A : valence' in data[str]:
            data[str] = f'mon : A : valence : {-1 * alpha}\n'
        elif 'mon : E : valence' in data[str]:
            data[str] = f'mon : E : valence : {-1 * alpha}\n'
        elif 'mol : Cl : phibulk' in data[str]:
            data[str] = f'mol : Cl : phibulk : {Cs}\n'
        elif 'mol : pol  : composition' in data[str]:
            data[str] = f'mol : pol  : composition : (X0)1(A){N_brush - 2}(E)1\n'
        elif 'mol : pol : theta' in data[str]:
            data[str] = f'mol : pol : theta : {theta}\n'
        

    # Создаем папку в зависимости от значения change_param
    folder_name = f'_pore_range_{range_param}_from_{min_range_value}_to_{max_range_value}'
    folder_name_ = folder_name.replace('.', '_')
    if not os.path.exists(folder_name_):
        os.makedirs(folder_name_)

    # Записываем изменения обратно в файл
    file_name_prefix = f'_pore_range_{range_param}_{round(change_param, 5)}'
    full_file_name = f'{file_name_prefix}_D_{round(D,2)}_N_{N_brush}_S_{round(S,2)}_theta_{round(theta,2)}.in'
    
    #избавляюсь от точек:
    last_dot_index = full_file_name.rfind('.')
    result = full_file_name[:last_dot_index].replace('.', '_') + full_file_name[last_dot_index:]
    
    new_file_path = os.path.join(folder_name_, result)
    
    
    with open(new_file_path, 'w') as file:
        file.writelines(data)

    #Считает намикс
    subprocess.call(['namics', os.path.abspath(new_file_path)])

    # Переносим посчитанные файлы
    
    # Создаем папку, где будут храниться output файлы
    folder_name_out = f'_pore_outfiles_range_{range_param}_from_{min_range_value}_to_{max_range_value}'
        
    # заменяем все точки на _
    folder_name_out__ = folder_name_out.replace('.', '_')
    os.makedirs(f'{folder_name_out__}', exist_ok=True)
    # Получаем список файлов в папке "output"
    files_in_output = os.listdir('output')
    # Путь до созданной папки
    to_folder_out = os.path.abspath(f'{folder_name_out__}')

    # Перемещаем каждый файл в созданную папку
    for file_in_output in files_in_output:
        file_path = os.path.join('output', file_in_output)
        if os.path.isfile(file_path):
            new_file_path = os.path.join(to_folder_out, file_in_output)
            os.rename(file_path, new_file_path)
            
    file_name_pro = os.path.join(to_folder_out, file_in_output)
    
    return file_name_pro