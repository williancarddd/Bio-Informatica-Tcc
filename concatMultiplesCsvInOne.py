import csv
import os
import pandas as pd

"""
concant multiples csv files in one
"""

def concat_csv_files():
    #all file names
    os.chdir('csv_result')
    relative_path = f"{os.getcwd()}/csv_result/"
    file_names = os.listdir('.')
    #combine all files in the list
    combined_csv = pd.concat([pd.read_csv(f) for f in file_names ])
    #export to csv
    combined_csv.to_csv( "combined_csv.csv", index=False, encoding='utf-8-sig')

if __name__ == '__main__':
    concat_csv_files()