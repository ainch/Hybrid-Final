import parse

def change_global_to_extern(path_file):
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

    lines_new = []

    for line in lines:
        line_new = line
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
                if ';' in line:
                    line_new = 'extern' + ' ' + line
        if ')' in line:
            num_parentheses = num_parentheses - 1
        if '}' in line:
            num_brace = num_brace - 1
        if '>' in line:
            num_inequality = num_inequality - 1
        lines_new.append(line_new)

    with open(path_file, 'w') as f:
        for line_new in lines_new:
            f.write(line_new)

def refactor():
    all_files = parse.get_all_files()
    for path_file in all_files:
        change_global_to_extern(path_file)