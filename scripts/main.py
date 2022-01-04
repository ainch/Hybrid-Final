import parse

if __name__=="__main__":
    all_files = parse.get_all_files()
    print(all_files)
    total_global_variables = parse.get_global_variables(all_files)
 
    for global_variable in total_global_variables:
        print('[',global_variable,']')