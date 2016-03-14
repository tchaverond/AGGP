# -*- coding: utf-8 -*-

import random
import copy
import matplotlib.pyplot as plt
from networkx import *
import heapq
import collections


class Simulation:
	"""
	Default Simulation constructor
	"""
	def __init__(self):
		# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
		# -__-__-__-__-__-__-__-__-__-__-                 Attributes                 -__-__-__-__-__-__-__-__-__-__- #
		# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
		# Nombre de graphes
		self.nb_graphs = 10
		# Nombre de sommets de chaque graphe
		self.nb_nodes = 10
		# Liste des ponderations des scores
		self.score_weights = [1, 1, 1]
		# Probas pour chaque type de mutation
		self.prob_cross_mutation = 0.5
		self.prob_mutation = 0.5
		self.prob_insertion = 0.5
		self.prob_deletion = 0.5
		# Coefficient de calcul de la fécondité
		self.coef_fertility = 1
		# Facteur maximum du nombre d'aretes initiales
		self.coef_ini_edges = 6
		# Verification du nombre maximal d'aretes
		self.coef_ini_edges = min( int(.5*(self.nb_nodes - 1)), self.coef_ini_edges)
		# Liste de genes (= liste des graphes)
		self.genome = [self._generate_graph(self.nb_nodes) for i in xrange(self.nb_graphs)]

	# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
	# -__-__-__-__-__-__-__-__-__-__-                  Methods                   -__-__-__-__-__-__-__-__-__-__- #
	# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #


	# -__-__-__-__-__-__-__-__-__-__-           Initilization Methods            -__-__-__-__-__-__-__-__-__-__- #
	"""
	Create a random graph with N nodes
	"""
	def _generate_graph(self, N):
		graph = dense_gnm_random_graph(N, random.randint(N, N * self.coef_ini_edges), seed=None)
		while not is_connected(graph) :
			graph = dense_gnm_random_graph(N, N+2, seed=None)
		return graph


	# -__-__-__-__-__-__-__-__-__-__-               Score Methods                -__-__-__-__-__-__-__-__-__-__- #


	"""
	Computes and returns the overall score of the networkx object G, using the weights associated to each component (all passed as parameters)
	"""
	def global_score(self, G) :

		return self.score_weights[0]*eval_degree_distrib(G) + self.score_weights[1]*eval_clustering_coef(G) + self.score_weights[2]*eval_aspl(G)



	"""
	Evaluates the degree distribution of the networkx object G passed as parameter,
	compares it to the power law k^-2.25 by computing the Smirnov statistic, and returns this statistic (score)
	"""
	def eval_degree_distrib(self, G) :

		all_degrees = degree_centrality(G).values()

		# Fonction de répartition (en s’assurant que les éléments sont dans le bon ordre -> voir a posteriori)
		list1 = collections.Counter(all_degrees).values()

		# Fonction de répartition théorique
		list2 = []
		for k in xrange(1,max(all_degrees)) : 
			list2.append(k**(-2.25))

		list2 = list2/sum(list2)

		# Statistique du test de Smirnov
		score = heapq.nlargest(1, (((a - b), a, b) for a, b in zip(list1, list2)))

		return score



	"""
	Evaluates the clustering coefficient (depending of the degree k) of the networkx object G passed as parameter,
	compares it to 1/k by computing the Smirnov statistic, and returns this statistic (score)
	"""
	def eval_clustering_coef(self, G) :

		# Calcul coefficients clustering
		clust_coeffs = clustering(G).values()

		# Fonction de répartition (en s’assurant que les éléments sont dans le bon ordre -> voir a posteriori)
		list1 = collections.Counter(clust_coeffs).values()

		# Fonction de répartition théorique
		list2 = []
		for k in xrange(1,max(clust_coeffs)) :
			list2.append(float(1.0/k))

		list2 = list2/sum(list2)

		# Statistique du test de Smirnov
		score = heapq.nlargest(1, (((a - b), a, b) for a, b in zip(list1, list2)))

		return score



	"""
	Evaluates the average length of the shortest path (aspl) of the networkx object G passed as parameter,
	compares it (relative gap) to log(log(N)) (N number of nodes), and returns this statistic (score)
	"""
	def eval_aspl(self, G) :

		l_avg = average_shortest_path_length(G)
		N = number_of_nodes(G)
		score = abs(l_avg - log(log(N))) / log(log(N))
		return score




	# -__-__-__-__-__-__-__-__-__-__-             Genetic Algorithm              -__-__-__-__-__-__-__-__-__-__- #


	"""
	Removes a random edge
	"""
	def remove_random_edge(self, indiv):
		index = random.randint(0,indiv.number_of_edges())
		indiv.remove_edge(indiv.edges[index])
		return 0

	"""
	Adds an edge between to random nodes
	"""
	def add_random_edge(self, indiv):
		index1 = random.randint(0,indiv.number_of_nodes())
		index2 = random.randint(0,indiv.number_of_nodes())
		indiv.add_edge(indiv.nodes[index1],indiv.nodes[index2])
		return 0

	"""
	Mutates punctually the graph, it changes the edges of a vertice chosen randomly
	"""
	def mutate(self, graph_list):
		for indiv in xrange(0,len(graph_list)):
			if random.random(0,1) < self.prob_mutation :
				index_edge = random.randint(0,indiv.number_of_edges)
				index_node = random.randint(0,indiv.number_of_edges)
				if random.random(0,1) < 0.5 :
					indiv.add_edge(indiv.edges[index][1], index.nodes[index_node])
				else :
					indiv.add_edge(index.nodes[index_node], indiv.edges[index][2])
				indiv.remove_edge(indiv.edges[index])
			if random.random(0,1) < self.prob_insertion:
				add_random_edge(indiv)
			if random.random(0,1) < self.prob_deletion:
				remove_random_edge(indiv)
		return 0

	"""
	Crossing - over
	"""
	def cross_mutate(self, graph_list):
		for indiv in xrange(0,len(graph_list)):
			if random.random(0,1) < self.prob_cross_mutation:
				other = indiv
				while other == indiv:
					other = random.randint(0,len(graph_list))
				rand1 = random.randint(0,indiv.number_of_edges())
				rand2 = random.randint(0,indiv.number_of_edges())
				while rand1 == rand2:
					rand2 = random.randint(0,indiv.number_of_edges())
				start = min(rand1,rand2)
				end = max(rand1,rand2)
				to_switch_indiv = indiv.edges[start:end]
				to_switch_other = other.edges[start:end]
				indiv.remove_edges_from(to_switch_indiv)
				other.add_edges_from(to_switch_indiv)
				other.remove_edges_from(to_switch_other)
				indiv.add_edges_from(to_switch_other)
		return 0



# -__-__-__-__-__-__-__-__-__-__-           Visualization Methods            -__-__-__-__-__-__-__-__-__-__- #
	"""
	Draw the n-th graph of the genome
	"""
	def draw_graph(self, n):
		draw_circular(self.genome[n])
		plt.show()

# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
# -__-__-__-__-__-__-__-__-__-__-                    Main                    -__-__-__-__-__-__-__-__-__-__- #
# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #


S = Simulation()
S.draw_graph(0)
print S.genome[0].number_of_edges()
print "Hello world !"



# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
