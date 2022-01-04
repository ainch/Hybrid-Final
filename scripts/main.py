import parse
import refactor

if __name__=="__main__":
    all_files = parse.get_all_files()
    total_global_variables = parse.get_global_variables(all_files)
 
    with open('out.txt', 'w') as f:
        for global_variable in total_global_variables:
            f.write(global_variable)

    refactor.refactor()