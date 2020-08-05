#!/usr/bin/python3
# -*- coding: utf-8 -*-
#

import os, re, sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy import stats

def run():
    result_folder = 'selection_analysis_results/'
    source_data_file_list = [result_folder+'VEGFA_table.csv', result_folder+'PGF_table.csv', result_folder+'VEGFB_table.csv', result_folder+'VEGFC_table.csv', result_folder+'VEGFD_table.csv']
    #source_data_file_list = ['analysis/VEGFA_table.csv']
    protein_dictionary = {}
    for file in source_data_file_list:
        protein_dictionary[file.split('/')[1].split('_')[0]] = pd.read_csv(file)

    #data = {'Site': range(1, 370)}
    #df = pd.DataFrame (data, columns = ['Site'])
 
    data = {}
    df = pd.DataFrame (data, columns = [])

    for protein, data in protein_dictionary.items():
        df[protein] = (data['beta'] / data['alpha']).rolling(11, center = True, min_periods = 1).mean()
    
    #fig, ax = plt.subplots()
    #ax.set(yscale="log")
    #sns_plot = sns.lineplot(data = df, palette = "tab10", linewidth = 1, ax = ax)
    sns_plot = sns.lineplot(data = df, palette = "tab10", linewidth = 1)

    sns_plot.set(xlabel = 'aa position', ylabel = 'Ka/Ks (11 aa window)')

    #sns_plot.savefig('figure.png', transparet = True)
    sns_plot.figure.savefig('figure.pdf', dpi = 300)
    #sns_plot.savefig('figure.eps', orientation = 'landscape', dpi = 300)
    sns_plot.figure.savefig('figure.svg')
    #
    #
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    print(df)
    print(df.shape)

if __name__ == '__main__':
    run()
