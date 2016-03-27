# -*- coding: utf-8 -*-

import random
import copy
import matplotlib.pyplot as plt
from networkx import *
import heapq
import collections
import math
import profile
import numpy as np


class Simulation:
	"""
	Default Simulation constructor
	"""
	def __init__(self):
		# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
		# -__-__-__-__-__-__-__-__-__-__-                 Attributes                 -__-__-__-__-__-__-__-__-__-__- #
		# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
		# Nombre de graphes
		self.nb_graphs = 150
		# Nombre de sommets de chaque graphe
		self.nb_nodes = 100
		# Liste des ponderations des scores
		self.score_weights = [1, 3, 0.5]
		# Probas pour chaque type de mutation
		self.prob_cross_mutation = 0.06
		self.prob_mutation = 0.22
		self.prob_insertion = 0.22
		self.prob_deletion = 0.22
		self.prob_swap = 0.22
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
			graph = dense_gnm_random_graph(N, N * self.coef_ini_edges, seed=None)
		return graph


	# -__-__-__-__-__-__-__-__-__-__-               Score Methods                -__-__-__-__-__-__-__-__-__-__- #


	"""
	Computes and returns the overall score of the networkx object G, using the weights associated to each component (all passed as parameters)
	"""
	def global_score(self, G) :

		return self.score_weights[0]*self.eval_degree_distrib(G) + self.score_weights[1]*self.eval_clustering_coef(G) + self.score_weights[2]*self.eval_aspl(G)



	"""
	Evaluates the degree distribution of the networkx object G passed as parameter,
	compares it to the power law k^-2.25 by computing the Smirnov statistic, and returns this statistic (score)
	"""
	def eval_degree_distrib(self, G) :

		all_degrees = degree_centrality(G).values()
		max_degree = int((number_of_nodes(G)-1)*max(all_degrees))

		# Fonction de répartition (en s’étant assuré que les éléments sont dans le bon ordre (OrderedDict))
		dict1 = collections.OrderedDict(sorted(collections.Counter(all_degrees).items()))
		list1 = dict1.values()

		norm = sum(list1)
		for i in xrange(len(list1)) :
			list1[i] = float(list1[i])/norm

		# Fonction de répartition théorique
		list2 = []
		for k in xrange(1,1+max_degree) : 
			list2.append(k**(-2.25))

		norm = sum(list2)
		for i in xrange(len(list2)) :
			list2[i] = list2[i]/norm


		# Statistique du test de Smirnov
		score = heapq.nlargest(1, (abs(a - b) for a, b in zip(list1, list2)))

		return score[0]



	"""
	Evaluates the clustering coefficient (depending of the degree k) of the networkx object G passed as parameter,
	compares it to 1/k by computing the Smirnov statistic, and returns this statistic (score)
	"""
	def eval_clustering_coef(self, G) :

		# Calcul coefficients clustering
		clust_coeffs = clustering(G).values()
		max_clust_coeff = int((number_of_nodes(G)-1)*max(clust_coeffs))

		# Fonction de répartition (en s’étant assuré que les éléments sont dans le bon ordre (OrderedDict))
		dict1 = collections.OrderedDict(sorted(collections.Counter(clust_coeffs).items()))
		list1 = dict1.values()

		norm = sum(list1)
		for i in xrange(len(list1)) :
			list1[i] = float(list1[i])/norm

		# Fonction de répartition théorique
		list2 = []
		for k in xrange(1,1+max_clust_coeff) :
			list2.append(float(1.0/k))

		norm = sum(list2)
		for i in xrange(len(list2)) :
			list2[i] = list2[i]/norm

		# Statistique du test de Smirnov
		score = heapq.nlargest(1, (abs(a - b) for a, b in zip(list1, list2)))
		
		# if max_clust_coeff is equal to 0 (can happen with only a few edges), list2 is empty, and so is score
		if score != [] :
			return score[0]
		else :
			return (self.nb_nodes + 1)



	"""
	Evaluates the average length of the shortest path (aspl) of the networkx object G passed as parameter,
	compares it (relative gap) to log(log(N)) (N number of nodes), and returns this statistic (score)
	"""
	def eval_aspl(self, G) :

		l_avg = average_shortest_path_length(G)
		N = number_of_nodes(G)
		score = abs(l_avg - math.log(math.log(N))) / math.log(math.log(N))
		return score

	"""
	Compute score for each graph in the list 
	"""
	def compute_all_score(self, graph_list):
		return [self.global_score(graph_list[i]) for i in xrange(len(graph_list))]

	"""
	Sort graph in list by their score
	"""
	def sort_by_score(self, graph_list):
		Tab = dict()
		res = []
		for i in xrange(len(graph_list)):
			Tab[self.global_score(graph_list[i])] = []
		for i in xrange(len(graph_list)):
			Tab[self.global_score(graph_list[i])].append(graph_list[i])
		sorted_tab_key = sorted(Tab.keys())
		for i in sorted_tab_key:
			for a in Tab[i]:
				res.append(a)
		return res

	# -__-__-__-__-__-__-__-__-__-__-             Genetic Algorithm              -__-__-__-__-__-__-__-__-__-__- #


	"""
	Removes a random edge
	"""
	def remove_random_edge(self, indiv):
		index = random.randint(0,indiv.number_of_edges()-1)
		indiv.remove_edge(indiv.edges()[index][0],indiv.edges()[index][1])
		return 0

	"""
	Adds an edge between to random nodes
	"""
	def add_random_edge(self, indiv):
		index1 = random.randint(0,indiv.number_of_nodes()-1)
		index2 = random.randint(0,indiv.number_of_nodes()-1)
		indiv.add_edge(indiv.nodes()[index1],indiv.nodes()[index2])
		return 0

	"""
	Mutates punctually the graph, it changes the edges of a vertice chosen randomly
	"""
	def mutate(self, graph):
		# Swap
		if random.random() < self.prob_swap :
			nx.double_edge_swap(graph, nswap=1, max_tries = 20)
		# Point mutation
		if random.random() < self.prob_mutation :
			index_edge = random.randint(0,graph.number_of_edges()-1)
			index_node = random.randint(0,graph.number_of_nodes()-1)
			if random.random() < 0.5 :
				graph.add_edge(graph.edges()[index_edge][0], graph.nodes()[index_node])
			else :
				graph.add_edge(graph.nodes()[index_node], graph.edges()[index_edge][1])
			graph.remove_edge(*graph.edges()[index_edge])
		# Insertion
		if random.random() < self.prob_insertion:
			self.add_random_edge(graph)
		# Deletion
		if random.random() < self.prob_deletion:
			self.remove_random_edge(graph)
		return graph

	"""
	Crossing - over
	"""
	def cross_mutate(self, graph_list):
		for indiv in graph_list :
			if random.random() < self.prob_cross_mutation:
				other = indiv
				while other == indiv:
					other = graph_list[random.randint(0,len(graph_list)-1)]
				rand1 = random.randint(0,min(indiv.number_of_edges(), other.number_of_edges())-1)
				rand2 = random.randint(0,min(indiv.number_of_edges(), other.number_of_edges())-1)
				while rand1 == rand2:
					rand2 = random.randint(0,indiv.number_of_edges()-1)
				start = min(rand1,rand2)
				end = max(rand1,rand2)
				to_switch_indiv = indiv.edges()[start:end]
				to_switch_other = other.edges()[start:end]
				# Take out common edges
				common = list(set(to_switch_other) ^ set(to_switch_indiv))
				to_switch_indiv = [i for i in to_switch_indiv if i in common]
				to_switch_other = [i for i in to_switch_other if i in common]
				indiv.remove_edges_from(to_switch_indiv)
				other.add_edges_from(to_switch_indiv)
				other.remove_edges_from(to_switch_other)
				indiv.add_edges_from(to_switch_other)
		return 0
	
	"""
	Create a new generation
	"""
	def new_generation(self):

		# Get current graph list and scores
		graph_list = self.genome
		score_list = self.compute_all_score(self.genome)

		# Create an array that will contain the new genome
		new_graph_list = [None]*len(self.genome)

		# Compute Fecondity
		F = [None]*len(self.genome)
		best_indiv = None
		best_score = 0
		for i in xrange(len(graph_list)):
			F[i] = math.exp(self.coef_fertility*(1/score_list[i]))
			if (best_score < (1/score_list[i]) ):
				best_score = 1/score_list[i]
				best_indiv = graph_list[i]
		t = sum(F)
		F = [F[i]/t for i in xrange(len(F))]

		# Reproduce
		nb_desc = np.random.multinomial(self.nb_graphs - 1, F)
		new_graph_index = 0
		for i in xrange(len(graph_list)):
			for j in xrange(nb_desc[i]):
				if j+1 == nb_desc[i] :
					new_graph_list[new_graph_index] = self.mutate(graph_list[i])
				else :
					new_graph_list[new_graph_index] = self.mutate(copy.deepcopy(graph_list[i]))
				new_graph_index += 1
		new_graph_list[new_graph_index] = copy.deepcopy(best_indiv)

		# Crossing over
		self.cross_mutate(new_graph_list)

		# Verify if graphs are connexe
		for i in xrange(len(new_graph_list)):
			while not is_connected(new_graph_list[i]):
				p = random.random()
				p_added = 0
				for j in xrange(len(F)):
					p_added += F[j]
					if p_added >= p:
						new_graph_list[i] = self.mutate(copy.deepcopy(self.genome[j]))
						break
		self.genome = new_graph_list
		return 0

# -__-__-__-__-__-__-__-__-__-__-           Visualization Methods            -__-__-__-__-__-__-__-__-__-__- #
	"""
	Draw the n-th graph of the genome
	"""
	def draw_graph(self, n):
		draw_spring(self.genome[n])
		plt.show()


	def monitor(self, n):
		print self.genome[n].number_of_edges()
		print is_connected(self.genome[n])
		self.draw_graph(n)

	def monitor2(self, graph):
		print graph.number_of_edges(), graph.number_of_nodes()
		print is_connected(graph)
		draw_spring(graph)
		plt.show()


	def monitor_score(self) :
		score = [[], [], [], []]
		for g in self.genome:
			score[0].append(self.score_weights[0]*self.eval_degree_distrib(g))
			score[1].append(self.score_weights[1]*self.eval_clustering_coef(g))
			score[2].append(self.score_weights[2]*self.eval_aspl(g))
			score[3].append(self.global_score(g))
		print "eval_degree_distrib:"
		print min(score[0]), max(score[0]), sum(score[0])/len(score[0])
		print "eval_clustering_coef:"
		print min(score[1]), max(score[1]), sum(score[1])/len(score[1])
		print "eval_aspl:"
		print min(score[2]), max(score[2]), sum(score[2])/len(score[2])
		print "total:"
		print min(score[3]), max(score[3]), sum(score[3])/len(score[3])

# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
# -__-__-__-__-__-__-__-__-__-__-                    Main                    -__-__-__-__-__-__-__-__-__-__- #
# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #


S = Simulation()

def test(S) :
	for i in xrange(101):
		if (i %10) == 0 :
			print i, len(S.genome), sum([g.number_of_edges() for g in S.genome])/len(S.genome)
			S.monitor_score()
		S.new_generation()
#S.draw_graph(0)
#S.draw_graph(49)
#S.monitor_score()
test(S)
#S.draw_graph(0)
#S.draw_graph(49)
#S.monitor_score()
#profile.run("test(S)")



# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
