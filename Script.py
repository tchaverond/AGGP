# -*- coding: utf-8 -*-

import random
import copy
import matplotlib
from networkx import *
import heapq
import collections


# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
# -__-__-__-__-__-__-__-__-__-__-                  Methods                   -__-__-__-__-__-__-__-__-__-__- #
# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #



# -__-__-__-__-__-__-__-__-__-__-               Score Methods                -__-__-__-__-__-__-__-__-__-__- #


"""
Computes and returns the overall score of the networkx object G, using the weights associated to each component (all passed as parameters)
"""
def global_score(G,score_weights) :

	return score_weights[0]*eval_degree_distrib(G) + score_weights[1]*eval_clustering_coef(G) + score_weights[2]*eval_aspl(G)



"""
Evaluates the degree distribution of the networkx object G passed as parameter,
compares it to the power law k^-2.25 by computing the Smirnov statistic, and returns this statistic (score)
"""
def eval_degree_distrib(G) :

	all_degrees = degree_centrality(G).values()

	# Fonction de répartition (en s’assurant que les éléments sont dans le bon ordre -> voir a posteriori)
	list1 = collections.Counter(all_degrees).values()

	# Fonction de répartition théorique
	list2 = []
	for k in xrange(1,max(all_degrees)) : 
		list2.append(k**(-2.25))

	list2 = list2/sum(list2)

	# Statistique du test de Smirnov
	score = heapq.nlargest(1, ((a - b), a, b) for a, b in zip(list1, list2))

	return score



"""
Evaluates the clustering coefficient (depending of the degree k) of the networkx object G passed as parameter,
compares it to 1/k by computing the Smirnov statistic, and returns this statistic (score)
"""
def eval_clustering_coef(G) :

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
	score = heapq.nlargest(1, ((a - b), a, b) for a, b in zip(list1, list2))

	return score



"""
Evaluates the average length of the shortest path (aspl) of the networkx object G passed as parameter,
compares it (relative gap) to log(log(N)) (N number of nodes), and returns this statistic (score)
"""
def eval_aspl(G) :

	l_avg = average_shortest_path_length(G)
	N = number_of_nodes(G)
	score = abs(l_avg - log(log(N))) / log(log(N))
	return score




# -__-__-__-__-__-__-__-__-__-__-             Genetic Algorithm              -__-__-__-__-__-__-__-__-__-__- #






# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
# -__-__-__-__-__-__-__-__-__-__-                    Main                    -__-__-__-__-__-__-__-__-__-__- #
# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #


print "Hello world !"



# -__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__-__- #
