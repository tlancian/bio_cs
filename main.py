import peeling
import pickle
import argparse
import os


#Input

parser = argparse.ArgumentParser()

parser.add_argument('d', help='Path to Graph', type=str)
parser.add_argument('C', help='C', type=float)

args = parser.parse_args()


#Execution

diff_graph, atlas = peeling.graph_handler(args.d)
results = peeling.C_peeling(diff_graph, args.C)

print("Graph: {}".format(args.d))
print("Solution: {}".format([atlas[node] for node in results[0]]))
print("Solution's Density: {}".format(results[1]))
