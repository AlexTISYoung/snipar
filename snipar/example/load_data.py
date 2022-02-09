import os
import shutil

def load_example_data(dest="example_data"):
    data_dir_name = "example_data"
    dir_name = os.path.dirname(os.path.realpath(__file__))
    data_dir = os.path.join(dir_name, data_dir_name)
    shutil.copytree(data_dir, dest)

if __name__ == '__main__':
    load_example_data()
