import os

def get_files(path_dir):
    list_dir = os.listdir(path_dir)
    list_paths_file = [ path_dir + '/' + path_file for path_file in list_dir ]
    return list_paths_file

def get_all_files():
    return get_files('inc') + get_files('src')

def extract_global_variables(path_file):
    lines = []
    with open(path_file) as f:
        while True:
            line = f.readline()
            if line:
                lines.append(line)
            else:
                break

    num_parentheses = 0
    num_brace = 0
    num_inequality = 0

    global_variables = []

    for line in lines:
        if '(' in line:
            num_parentheses = num_parentheses + 1
        if '{' in line:
            num_brace = num_brace + 1
        if '<' in line:
            num_inequality = num_inequality + 1


        if num_parentheses == 0 and num_brace == 0 and num_inequality == 0:
            if '#' not in line and 'extern' not in line and line[0] != '/' \
                and line != '' and line != '\t' \
                and '*/' not in line and '/*' not in line \
                and 'typedef' not in line and 'namespace' not in line:
                global_variables.append(line)

        if ')' in line:
            num_parentheses = num_parentheses - 1
        if '}' in line:
            num_brace = num_brace - 1
        if '>' in line:
            num_inequality = num_inequality - 1
    return global_variables

def get_global_variables(all_files):
    total_global_variables = []
    for path_file in all_files:
        total_global_variables = total_global_variables + extract_global_variables(path_file)
    return total_global_variables