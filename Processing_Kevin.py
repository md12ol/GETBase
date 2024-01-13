import math
import os
from operator import itemgetter
import numpy as np
import copy

import matplotlib.pyplot as plt
import numpy as np
# Need to install on linux before use, and include as packages in PyCharm
from graphviz import Graph  # Network visualizer (nodes and edges)

raw_data_path = "./Output/"
figure_path = "./Figures/"
main_filename = "best01.dat"
num_runs = 25
lower_better = False
precision = 6
col_width = 10 + precision
num_nodes = 256


# Go through the best##.dat file and gathers all important results from the file.
def get_data(file_path: str):
    # Generate empty lists to hold all the pertinent information from the file.
    fits = []  # Fitness values
    nets_num_edges = []  # Number of edges in each network
    nets_edge_lists = []  # Edge lists for the generated networks
    one_net_lists = []  # Collector for the edge lists of one network

    # Flags for current position in the file while processing it.
    graph_flag = False  # Graph metadata flag
    edge_lists_flag = False  # Graph's edge lists flag

    with open(file_path) as f:  # Open file best##.dat
        lines = f.readlines()  # Read the file line-by-line and place in a list of lines
        for line in lines:
            if line.__contains__("Best Fitness: "):  # We found a fitness value!
                line = line.split(": ")[1]  # i.e. line = "20.8\n"
                # Removes all the whitespace and newline/carriage return from the right side of the string.
                line = line.rstrip()  # i.e. line = "20.8"
                fits.append(float(line))
                if len(one_net_lists) > 0:  # We found the end of the previous Graph's edge lists!
                    edge_lists_flag = False
                    nets_edge_lists.append(one_net_lists[:-1])  # Remove blank line list from edge lists
                    pass
                pass
            elif line.__contains__("Graph"):  # We found the Graph!  Time to process the Graph's metadata!
                graph_flag = True
                pass
            elif line.__contains__("W Hist:"):  # The next line contains the first of the edge lists
                edge_lists_flag = True
                # Could get the weight histogram here, but we don't want it right now.
                graph_flag = False
                one_net_lists = []
                pass
            elif graph_flag:  # Process the Graph's metadata!
                if line.__contains__("Edges: "):  # Get the number of edges in the network!
                    nets_num_edges.append(int(line.split(": ")[1].rstrip()))
                    pass
                pass
            elif edge_lists_flag:  # Process the Graph's edge lists!
                if line.rstrip() != '':  # Check if the line is empty
                    one_net_lists.append([int(i) for i in line.rstrip().split(" ")])  # Generate one node's adjacency list
                else:
                    one_net_lists.append([])
                pass
            pass
        nets_edge_lists.append(one_net_lists)
        pass

    # data[run_num][0] = run number
    # data[run_num][1] = run's fitness
    # data[run_num][2] = number of edges in the best network for this run
    # data[run_num][3][:][:] = edge lists for the best network for this run
    data = [[idx + 1, fits[idx], nets_num_edges[idx], nets_edge_lists[idx]] for idx in range(num_runs)]
    data.sort(key=itemgetter(1))  # Ascending
    if not lower_better:
        data.reverse()  # Descending
        pass

    return data


# Converts the adjacency lists of a network to a list of all the edges (pairs of nodes) in the network.
def adj_lists_to_list_of_edges(node_lists):
    edge_lists = []
    for fr in range(num_nodes):  # For each node
        for to in node_lists[fr]:  # For each node adjacent to this node
            edge_lists.append([fr, to])  # A pair representing one edge is added
            pass
        pass
    return edge_lists


# Generates a network with the provided list of edges.
def make_graph(list_of_edges: [], out_file: str):
    g = Graph(engine="sfdp")
    e_cout = 0

    # Graph attributes
    g.graph_attr.update(dpi='300', size="6,6", splines='true', overlap='false')
    g.node_attr.update(color='black', shape='point', fixedsize='true', width='0.01')
    g.edge_attr.update(color='black', penwidth='0.25')

    # Add the nodes
    for n in range(num_nodes):
        if n == 0:
            g.node(str(n), color='red')
            pass
        else:
            g.node(str(n))
        pass

    # Add the edges
    for idx, d in enumerate(list_of_edges):
        if d[0] < d[1]:  # Not directed, thus we don't want to repeat edges (i.e. edge 2 -> 3 AND 3 -> 2 exists)
            g.edge(str(d[0]), str(d[1]), color='black')
            e_cout += 1
            pass
        pass
    g.render(filename=out_file, directory=figure_path, cleanup=True, format='png')
    del g
    print("Made network: " + out_file + " with " + str(e_cout) + " Edges")
    pass


# Generates a boxplot with the lists of results provided in data.
def make_boxplot(data: [], out_path: str):
    data = copy.deepcopy(data)
    plt.style.use("seaborn-v0_8")
    plt.rc('xtick', labelsize=10)
    plt.rc('ytick', labelsize=10)

    # Create a figure object
    f = plt.figure()
    f.set_figheight(4.5)
    f.set_figwidth(10)
    plot = f.add_subplot(111)

    # Calculate x values and generate the plot.
    xs = [i + 1 for i in range(len(data))]
    vp = plot.violinplot(data, xs, showmedians=True, widths=0.85)

    # Make it fancy!
    for pc in vp["bodies"]:
        pc.set_facecolor("#5770DB")
        pc.set_linewidth(1)
        pc.set_edgecolor("black")
        pass
    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians'):
        vc = vp[partname]
        vc.set_linewidth(1)
        vc.set_alpha(1)
        vc.set_color("#5770DB")
        pass

    plot.grid(visible="True", axis="y", which='minor', color="white", linewidth=0.5)
    plt.minorticks_on()

    f.suptitle("Experiment Results", fontsize=14)
    plot.set_xlabel("Experiment", fontsize=12)
    plt.xticks([x + 1 for x in range(len(data))])

    f.subplots_adjust(bottom=0.1, top=0.93, left=0.03, right=0.99)
    f.savefig(out_path + "Boxplot.png", dpi=500)

    plt.close()
    pass


# Generates a table "experiment_table.dat" that indicates which experiment number is associated with which folder.
def make_exp_table(folder_names):
    with open("./Figures/experiment_table.dat", "w") as f:
        f.write(str("Exp Labels").ljust(col_width))
        f.write("Folder Name")
        f.write("\n")
        for idx, fold in enumerate(folder_names):
            f.write(str("EXP" + str(idx + 1)).ljust(col_width))
            f.write(fold)
            f.write("\n")
            pass
    pass


def main():
    print("START")
    folder_names = os.listdir(raw_data_path)
    make_exp_table(folder_names)

    # Fill the data list with data from all experiments in the raw_data_path folder.
    data = []
    for fold in folder_names:
        dat = get_data(raw_data_path + fold + "/" + main_filename)
        data.append(dat)
        pass

    # data now looks like this:
    # data[exp_num][run_idx][0] = run number
    # data[exp_num][run_idx][1] = run's fitness
    # data[exp_num][run_idx][2] = number of edges in run's best network
    # data[exp_num][run_idx][3][:][:] = edge lists for the best network for this run

    # We need just the fitness values for each experiment.  This will put them in data2:
    data2 = []
    exp_dat = []
    for dat in data:
        for run in dat:
            exp_dat.append(run[1])
            pass
        data2.append(exp_dat)
        exp_dat = []
        pass
    # data2[exp_num][:] = fitness values for all the runs in experiment exp_num

    make_boxplot(data2, figure_path)

    # Now, generate the best network for each experiment.
    for idx, dat in enumerate(data):
        make_graph(adj_lists_to_list_of_edges(dat[0][3]), "EXP" + str(idx + 1).zfill(2) + "_best_network")
        pass

    print("END")
    pass


main()
