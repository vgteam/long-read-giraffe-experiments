#!/usr/bin/env python3
"""Plot a line plot version of a cumulative histogram.

Accepts a column name with a min/max/step.
Data is found via TSV config file.
Each line in the config file has the format:
<statistics TSV>\t<format index>\t<label>

Format indexes 0-7 use COLORS with linestyle -
Next 8-15 use linestyle --, 16-23 :
Don't use more than that, seriously
"""

import argparse  # for command-line argument parsing
import csv  # for reading TSV files
import numpy as np  # for numerical operations
import matplotlib.pyplot as plt  # for plotting
from typing import List  # for type hints

# https://davidmathlogic.com/colorblind/
COLORS = ['#E69F00', '#56B4E9', '#009E73', '#0072B2', 
          '#D55E00', '#CC79A7', '#F0E442', '#000000']

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot cumulative histograms as line plots")
    parser.add_argument("--column-name", "-c", type=str,
                        help="Column name to plot")
    parser.add_argument("--min-val", "-m", type=float,
                        help="Minimum value for histogram bins")
    parser.add_argument("--max-val", "-M", type=float,
                        help="Maximum value for histogram bins")
    parser.add_argument("--step", "-s", type=float,
                        help="Step size for histogram bins")
    parser.add_argument("--title", "-t", type=str, help="Title for the plot")
    parser.add_argument("--output-file", "-o", type=str, help="Output file")
    parser.add_argument("config", type=str, help="Path to the config file")
    return parser.parse_args()

def read_config(config_file: str, column_name: str) -> List[tuple]:
    """Read and use a config file

    Parameters
    ----------
    config_file : str
        Path to the config file.
    column_name : str
        Name of the column to read from the TSV files.
    
    Returns
    -------
    List[tuple]
        Each tuple contains ([values], color, linestyle, label).
    """
    config = []
    with open(config_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) != 3:
                raise ValueError("Config file must have three columns: "
                                 "<file>\\t<format i>\\t<label>")
            values = [float(row[column_name])
                      for row in csv.DictReader(parts[0], delimiter='\t')]
            format_index = parts[1]

            # Interpret format index (see file docstring)
            try:
                format_index = int(parts[1])
                if 0 <= format_index <= 7:
                    linestyle = '-'
                elif 8 <= format_index <= 15:
                    linestyle = '--'
                elif 16 <= format_index <= 23:
                    linestyle = ':'
                else:
                    raise ValueError()
            except ValueError:
                raise ValueError('Format index (col 2) must be an integer 0-23')
            color = COLORS[format_index % 8]
            config.append((values, color, linestyle, parts[2]))
    return config

def plot_cumulative_line(to_plot: List[tuple], column_name: str,
                         bins: np.ndarray, title: str,
                         output_file: str) -> None:
    """Plot a cumulative histogram as a line plot.

    Parameters
    ----------
    to_plot : List[tuple]
        List of tuples containing ([values], color, linestyle, label).
    column_name : str
        Name of column for labels & image name.
    bins : np.ndarray
        Histogram bins.
    title : str
        Title for the plot.
    output_file : str
        Output file name.
    """
    plt.figure(figsize=(8, 4.5))

    bin_centers = (bins[:-1] + bins[1:]) / 2
    for values, color, linestyle, label in to_plot:
        counts, _ = np.histogram(values, bins=bins)
        if len(values) == 0:
            print(f"Warning: No values to plot for label '{label}'")
        else:
            cum_counts = np.cumsum(counts[::-1])[::-1] / len(values) * 100
            plt.plot(bin_centers, cum_counts, color=color, 
                    linestyle=linestyle, label=label)
    
    plt.grid()
    plt.title(title)

    plt.xlabel(column_name, fontsize=13)
    plt.ylabel('% reads ≥', fontsize=13)
    plt.legend(loc='bottom left')

    plt.savefig(output_file)
    plt.close()

if __name__ == '__main__':
    args = parse_args()

    max_bin = args.max_val + args.step * 1 # Include last bin as well
    bins = np.arange(args.min_val, max_bin, args.step)
    config = read_config(args.config, args.column_name)

    plot_cumulative_line(config, args.column_name, bins, 
                         args.title, args.output_file)
