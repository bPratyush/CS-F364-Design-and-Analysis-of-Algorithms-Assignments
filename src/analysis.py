import re
import matplotlib.pyplot as plt
import time
import numpy as np

def parse_output(file):
    with open(file, 'r') as f:
        lines = f.readlines()
    
    cliques = []
    for line in lines:
        if "Maximal Clique:" in line:
            clique = list(map(int, re.findall(r'\d+', line)))
            cliques.append(clique)
    
    largest_clique_size = max(len(clique) for clique in cliques)
    total_cliques = len(cliques)
    size_distribution = [len(clique) for clique in cliques]
    
    return largest_clique_size, total_cliques, size_distribution

def plot_histogram(data, title, xlabel, ylabel, filename):
    if all(isinstance(i, int) for i in data):
        bins = range(1, max(data) + 2)
    else:
        bins = np.linspace(min(data), max(data), 30)
    
    plt.hist(data, bins=bins, edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(filename)
    plt.clf()

datasets = ['output.txt']  # Add all output files here
execution_times = []
largest_cliques = []
total_cliques = []
size_distributions = []

for dataset in datasets:
    start_time = time.time()
    largest_clique_size, total_cliques_count, size_distribution = parse_output(dataset)
    end_time = time.time()
    
    execution_times.append(end_time - start_time)
    largest_cliques.append(largest_clique_size)
    total_cliques.append(total_cliques_count)
    size_distributions.extend(size_distribution)

# Plot execution time histogram
plot_histogram(execution_times, 'Execution Time', 'Time (s)', 'Frequency', 'execution_time_histogram.png')

# Plot largest clique size histogram
plot_histogram(largest_cliques, 'Largest Clique Size', 'Clique Size', 'Frequency', 'largest_clique_size_histogram.png')

# Plot total number of maximal cliques histogram
plot_histogram(total_cliques, 'Total Number of Maximal Cliques', 'Number of Cliques', 'Frequency', 'total_cliques_histogram.png')

# Plot distribution of different size cliques histogram
plot_histogram(size_distributions, 'Distribution of Different Size Cliques', 'Clique Size', 'Frequency', 'clique_size_distribution_histogram.png')