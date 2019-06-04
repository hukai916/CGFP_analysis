"""
Try to repeat/customize the heatplot figure with Python.
Usage:
python plot.py Raw_data/sample_cluster.txt Raw_data/metagenomes.txt

sample_cluster.txt is downloaded from EFI-CGFP sever, metagenomes.txt is retrieved from there too.
Results will be stored in Plots folder.
"""

import sys
import pandas as pd
import seaborn as sb
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from collections import OrderedDict

group_label = ["anterior nares", "buccal mucosa", "supragingival plaque", "tongue dorsum", "stool", "posterior fornix"]
group_color = ['orange', 'darkblue', 'blue', 'deepskyblue', 'brown', 'green']

file_data = sys.argv[1]
file_meta = sys.argv[2]
meta_dict = OrderedDict()

def _get_groups(file_meta):
    """
    Regroup df.
    """
    for line in open(file_meta):
        colname, group = line.split(":")
        group = group.strip()
        if not group in meta_dict:
            meta_dict[group] = [colname]
        else:
            meta_dict[group].append(colname)
    return(meta_dict)

def _get_new_group(df):
    test = _get_groups(file_meta)
    df_new = pd.DataFrame()

    # Reorder df by group ids:
    for key in group_label:
        cols = test[key]
        df_tem = df[cols].copy()
        df_tem.columns = [key] * df_tem.shape[1]
        df_new = pd.concat([df_new, df_tem], axis=1)
    return(df_new)

def plot_hp(df_new, name):
    """
    Plot heatmap given df.
    """

    cbar_kws = {"ticks":[1,10e-1,10e-2,10e-3,10e-4,10e-5]}
    g = sb.heatmap(df_new + 0.0001,
                        cbar_kws=cbar_kws,
                        vmin=0.0001,
                        vmax=1,
                        norm=LogNorm(vmin=df_new.as_matrix().min(), vmax=df_new.as_matrix().max()),
                        cmap='jet',
                        xticklabels=False)
    ax = plt.gca()
    fig = plt.gcf()
    total_col_num = df_new.shape[1]
    positions = str(ax.get_position()).split()
    x_start    = float(positions[0].split('=')[1][:-1])
    x_size     = float(positions[2].split('=')[1][:-1]) - x_start

    #x_size  = 0.62
    #x_start = 0.125

    plt.title("Abundances for colored ssn; Median Method")

    # Below add grouping according to metagenome postions:
    for i in range(len(group_label)):
        num = sum([1 for x in list(df_new.columns) if x == group_label[i]])
        gap  = (x_size) * num/total_col_num
        fig.add_axes([x_start,0.088,gap,0.02], facecolor=group_color[i])
        x_start = x_start + gap

        # Below removes ticks from plot:
        plt.tick_params(
            axis='both',       # changes apply to the both x and y-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            left=False,        # ticks along the left edge are off
            labelbottom=False,
            labelleft=False)

        # Below removes box plot frame:
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        plt.xlabel('\n'.join(group_label[i].split()))
    fig.savefig('Plots/' + name + '.svg')
    plt.close()

def boxplot_input(df_new):
    """
    Convert df into boxplot format.
    """
    start = 0
    data_list = []
    for i in range((len(group_label))):
        end = start + sum([1 for x in list(df_new.columns) if x == group_label[i]])
        data = df_new.iloc[:, start:end].values.tolist()[0]
        start = end
        data_list.append(data)
    return(data_list)

def plot_boxplot(box_input, name):
    fig, ax = plt.subplots()
    ax.boxplot(box_input, labels=['\n'.join(x.split()) for x in group_label])
    fig.savefig('Plots/Cluster_' + name + '.svg')
    plt.close()


df = pd.read_csv(file_data, sep='\t', header=0)
df = df.rename(columns={"Feature \ Sample": "Cluster Number"})
df = df.set_index("Cluster Number")

list1 = ['1', '4', '5', '6', '13', '18', '24', '25', '28', '31', '32', '36']
list2 = ['S6', 'S7', 'S9', 'S13', 'S17', 'S19', 'S24', 'S32', 'S34', 'S37',
         'S39', 'S43', 'S47', 'S48', 'S58', 'S60', 'S62']

# Plot for heatmaps:
nameList = ['Cluster', 'Singleton', 'Combined']
i = 0
for data in [list1, list2, list1 + list2]:
    df_selected = df[df.index.isin(data)]
    df_new = _get_new_group(df_selected)
    plot_hp(df_new, nameList[i])
    i = i + 1

# Plot for boxplots:
for cluster in list1 + list2:
    df_selected = df[df.index.isin([cluster])]
    df_new = _get_new_group(df_selected)
    box_input = boxplot_input(df_new)
    plot_boxplot(box_input, cluster)
