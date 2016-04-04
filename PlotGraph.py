# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from networkx import *
import pickle, os, sys

#arg: <input filename> <output name>
def plot_graph(G, output_name):
	monitor2(G, output_name)

def monitor2(graph, output_name):
	print "Edges : %d, nodes : %d."%(graph.number_of_edges(), graph.number_of_nodes())
	print "Connected : %s"%is_connected(graph)
	fig = plt.gcf()
	fig.set_size_inches(30,30)
	draw_spring(graph)
	plt.savefig("best_graph_"+output_name+".pdf", bbox_inches='tight')
	plt.clf()
	
def draw_graph(n):
	draw_spring(n)
	plt.show()

with open(sys.argv[1], 'rb') as input:
	graph = pickle.load(input)

print graph

plot_graph(graph, sys.argv[2])