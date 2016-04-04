# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import os, sys

#arg: <input filename> <prefix output filename> <windowed mode: 0 or 1> <score 1 weight> <score 2 weight> <score 3 weight>
scores_filename = sys.argv[1]
prefix = sys.argv[2] 
windows_plot = int(sys.argv[3]) 
score_weights = [float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])]
output_name = "end"

f = open(scores_filename, "r")
line = f.readline()
scores_1_mean = []
scores_1_var = []
scores_2_mean = []
scores_2_var = []
scores_3_mean = []
scores_3_var = []
global_score_mean = []
global_score_var = []
global_score_min = []
while line != "":
	scores_list = line.split()
	scores_1_temp = []
	scores_2_temp = []
	scores_3_temp = []
	global_score_temp = []
	for i in scores_list:
		tmp_score = i.split("/")
		scores_1_temp.append(score_weights[0]*float(tmp_score[0]))
		scores_2_temp.append(score_weights[1]*float(tmp_score[1]))
		scores_3_temp.append(score_weights[2]*float(tmp_score[2]))
		global_score_temp.append(float(tmp_score[3]))
	scores_1_mean.append(sum(scores_1_temp)/len(scores_1_temp))
	scores_1_var.append(np.std(scores_1_temp))
	scores_2_mean.append(sum(scores_2_temp)/len(scores_2_temp))
	scores_2_var.append(np.std(scores_2_temp))
	scores_3_mean.append(sum(scores_3_temp)/len(scores_3_temp))
	scores_3_var.append(np.std(scores_3_temp))
	global_score_mean.append(sum(global_score_temp)/len(global_score_temp))
	global_score_var.append(np.std(global_score_temp))
	global_score_min.append(min(global_score_temp))
	line = f.readline()
p1 = plt.plot(scores_1_mean, label = 'Mean of the score on degree distribution', color = 'blue')
p2 = plt.plot(scores_2_mean, label = 'Mean of the score on clustering coefficient', color = 'green')
p3 = plt.plot(scores_3_mean, label = 'Mean of the score on average shortest path length', color = 'red')
plt.fill_between(range(len(scores_1_mean)), np.array(scores_1_mean)+np.array(scores_1_var), np.array(scores_1_mean)-np.array(scores_1_var), facecolor = 'blue', alpha = 0.5)
plt.fill_between(range(len(scores_2_mean)), np.array(scores_2_mean)+np.array(scores_2_var), np.array(scores_2_mean)-np.array(scores_2_var), facecolor = 'green', alpha = 0.5)
plt.fill_between(range(len(scores_3_mean)), np.array(scores_3_mean)+np.array(scores_3_var), np.array(scores_3_mean)-np.array(scores_3_var), facecolor = 'red', alpha = 0.5)
plt.xlabel('Generation')
plt.ylabel('Score')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(18, 18)
plt.savefig(prefix+"scores_evolution_"+output_name+".pdf", bbox_inches = 'tight')
axes = plt.gca()
axes.set_xlim([0,len(scores_1_mean)])
plt.savefig(prefix+"scores_evolution_"+output_name+".pdf", bbox_inches='tight')
if (windows_plot):
	plt.show()
else :
	plt.clf()
p4 = plt.plot(global_score_mean, label = 'Mean of the weighted average of the scores', color = 'blue')
plt.fill_between(range(len(global_score_mean)), np.array(global_score_mean)+np.array(global_score_var), np.array(global_score_mean)-np.array(global_score_var), facecolor = 'blue', alpha = 0.5)
p5 = plt.plot(global_score_min, label = 'Min of the weighted average of the scores', color = 'green')
plt.xlabel('Generation')
plt.ylabel('Score')
plt.legend()
fig = plt.gcf()
fig.set_size_inches(18, 18)
axes = plt.gca()
axes.set_xlim([0,len(scores_1_mean)])
plt.savefig(prefix+"global_score_evolution_"+output_name+".pdf", bbox_inches='tight')
if (windows_plot):
	plt.show()
else :
	plt.clf()
f.close()