import sys
import os
import subprocess
sys.path.append('/home/tpopova/prj/polymer_brush')

import numpy as np

def generate_files(template_anneal:str = 'annealing_brush_temp.in',
                   template_quecnhed:str = 'quecnhed_brush_temp.in',
                   #Какой шаблон меняем?
                   current_name:str = 'annealing_brush_temp.in',

                   #название параметра
                   range_param:str = 'dpK',
                   range_param_q:str = 'dpK',

                   N_brush: int = 200,
                   S: int = 100,
                   pK_brush: float = 5,
                   Cs: float = 0.001,
                   delta_pK: float = -0.8,
                   min_range_value: float = 4,
                   max_range_value: float = 4,
                   alpha_z_average: float = 0.5,
                   ):  # sourcery skip: low-code-quality
    #какие параметры меняем?
    change_param = delta_pK
    change_param_q = alpha_z_average

    #pH раствора
    pH_b = pK_brush - delta_pK

    #pH для намикса
    pH_namics = (1 + 10**(pH_b))**(-1)

    #число слоев
    N_layers = N_brush+20

    theta = N_brush/S

    with open(current_name, 'r') as file:
        data = file.readlines()

    # Изменяем указанные параметры
    for str in range(len(data)):
        if 'state : H3O : alphabulk' in data[str]:
            data[str] = f'state : H3O : alphabulk : {pH_namics}\n'
        elif 'lat : flat : n_layers' in data[str]:
            data[str] = f'lat : flat : n_layers : {N_layers}\n'
        elif 'reaction : weak : pK' in data[str]:
            data[str] = f'reaction : weak : pK : {pK_brush}\n'
        elif 'mol : Na : phibulk' in data[str]:
            data[str] = f'mol : Na : phibulk : {Cs}\n'
        elif 'mol : Cl : phibulk' in data[str]:
            data[str] = f'mol : Cl : phibulk : {Cs}\n'
        elif 'mol : brush  : composition' in data[str]:
            data[str] = f'mol : brush  : composition : (X0)1(A){N_brush-2}(G)1\n'
        elif 'mol : brush : theta' in data[str]:
            data[str] = f'mol : brush : theta : {theta}\n'

        #Quenched brush
        elif 'mon : A : valence' in data[str]:
            data[str] = f'mon : A : valence : {-1 * alpha_z_average}\n'
        elif 'mon : G : valence' in data[str]:
            data[str] = f'mon : G : valence : {-1 * alpha_z_average}\n'
        elif 'mon : X0 : valence' in data[str]:
            data[str] = f'mon : X0 : valence : {-1 * alpha_z_average}\n'

    # Создаем папку в зависимости от значения current_name
    folder_name = f'_annealing_range_{range_param}_from_{min_range_value}_to_{max_range_value}' if current_name == template_anneal else f'_quecnhed_range_{range_param_q}_from_{min_range_value}_to_{max_range_value}'
    folder_name_ = folder_name.replace('.', '_')
    if not os.path.exists(folder_name_):
        os.makedirs(folder_name_)

    # Записываем изменения обратно в файл
    file_name = current_name if file.name == template_anneal else 'quecnhed_brush_temp.in'
    file_name_prefix = f'_annealing_range_{range_param}_{round(change_param, 4)}' if file.name == template_anneal else f'_quenched_range_{range_param_q}_{round(change_param_q, 4)}'
    new_file_path = os.path.join(folder_name_, f'{file_name_prefix}_pH_b_{pH_b}_Cs_{Cs}_N_{N_brush}_theta_{theta}.in').replace('.', '_', 3)
    with open(new_file_path, 'w') as file:
        file.writelines(data)

    #Считает намикс
    subprocess.call(['namics', os.path.abspath(new_file_path)])

    # Переносим посчитанные файлы
    
    # Создаем папку, где будут храниться output файлы
    folder_name_out = f'_anneal_outfiles_range_{range_param}_from_{min_range_value}_to_{max_range_value}' if current_name == template_anneal else f'_quecnhed_outfiles_range_{range_param_q}_from_{min_range_value}_to_{max_range_value}'
        
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
    
    return pH_b, file_name_pro