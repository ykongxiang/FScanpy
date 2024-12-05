import os
import pkg_resources

def get_test_data_path(filename: str) -> str:
    return pkg_resources.resource_filename('FScanpy', f'data/test_data/{filename}')

def list_test_data() -> list:
    data_dir = pkg_resources.resource_filename('FScanpy', 'data/test_data')
    return os.listdir(data_dir)