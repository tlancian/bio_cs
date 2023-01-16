import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import argparse
import peeling

def degree_distribution(G, subgraph):
    
    degrees = {elem:0 for elem in subgraph}
    
    for i in range(len(subgraph)):
        for j in range(i+1, len(subgraph)):
            try:
                degrees[subgraph[i]] += G[subgraph[i]][subgraph[j]]
                degrees[subgraph[j]] += G[subgraph[i]][subgraph[j]]
            except TypeError:
                pass
    return list(degrees.values())

dataset = "transcr"
classes = ["CPTAC", "TCGA"]


parser = argparse.ArgumentParser()

parser.add_argument('g1', help='Graph 1', type=str)
parser.add_argument('g2', help='Graph 1', type=str)
parser.add_argument('d', help='difference graph', type=str)
parser.add_argument('C', help='C', type=float)

args = parser.parse_args()


diff_graph, atlas_diff = peeling.graph_handler(args.d)
results = peeling.C_peeling(diff_graph, args.C)


sol = [atlas_diff[node] for node in results[0]]

#Dens_A
G_A = peeling.graph_handler(args.g1)
degrees_A = degree_distribution(G_A, sol)

#Dens_B
G_B = peeling.graph_handler(args.g2)
degrees_B = degree_distribution(G_B, sol)

degrees = pd.DataFrame([["degree"]*(len(degrees_A)+len(degrees_B)),
                    ["class_1"]*len(degrees_A)+["class_2"]*len(degrees_B),
                    degrees_A+degrees_B]).transpose()

degrees.columns = ["type","Subtype","Degree"]
degrees.type = degrees.type.astype("category")
degrees.Subtype = degrees.Subtype.astype("category")
degrees.Degree = degrees.Degree.astype("float64")


deg = sns.violinplot(x="type", y="Degree", hue="Subtype", split=True, data=degrees, inner = "quartile", palette = ["lightseagreen", "red"])
deg.set_title("class_1 - class_2 Contrast Subgraph")
deg.set(xticklabels=[], xticks = [], xlabel=None)
plt.show()
