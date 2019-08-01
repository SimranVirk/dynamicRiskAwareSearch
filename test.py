import networkx as nx 
import heapq
import numpy as np
from scipy.special import erfinv, erfc
import scipy.integrate
import copy
from collections import defaultdict
import generate_graphs as gen
import matplotlib.pyplot as plt
import math

threshold = 0.6

class PathObject:
	def __init__(self, path, mean ,var):
		self.path = path
		self.mean = mean
		self.var = var

	def __lt__(self, other):
		"""
		Less than function. Returns whether self is less than another path object
		"""

		adomb = dom(self, other)
		bdoma = dom(other, self)

		if adomb:
			return True
		elif bdoma:
			return False
		else:
			return self.mean < other.mean

	#***
	#def an equals method


def dom(p1, p2):
	"""
	Returns whether pathobject p1 dominates path object p2
	"""
	rhs = p2.mean + math.sqrt(2 * (p1.var + p2.var)) * erfinv(1 - 2 * threshold)
	return p1.mean < rhs

def betterPathToNode(closed, p):
	last = p.path[-1]
	for i in closed[last]:
		if dom(i, p):
			return True
	return False

def betterPathToGoal(closed, goal, p):
	for i in closed[goal]:
		if dom(i, p):
			return True
	return False

def constructSubgraph(graph, closed, goal, start):
	R = nx.Graph()
	#contains a map of nodes to paths from them to the goal as well as the mean/var of those paths
	#if closed is like g in a*, path_sets is like h in a*
	path_sets = defaultdict(set)

	R.add_node(start)
	for p in closed[goal]:
		curm = p.mean
		curv = p.var
		path = p.path

		path_sets[path[0]].add(p)

		for i in range(1, len(path)):
			if i < len(path) - 1:
				curm -= graph.edges[path[i-1], path[i]]['mean']
				curv -= graph.edges[path[i-1], path[i]]['var']

				path_sets[path[i]].add(PathObject(path[i:], curm, curv))

			if path[i] not in R.nodes:
				R.add_node(path[i])

			if (path[i - 1], path[i]) not in R.edges:
				R.add_edge(path[i - 1], path[i], 
					mean = graph.edges[path[i - 1], path[i]]['mean'], 
					var = graph.edges[path[i - 1], path[i]]['var'])


	final_sets = {}
	for k,v in path_sets.items():
		final_sets[k] = list(v)
	return (R, final_sets)		

def get_successors(graph, node):
	print(graph.successors(node))
	return graph.successors(node)

def get_mean(graph, node1, node2):
	return graph.edges[node1, node2]['mean']

def get_var(graph, node1, node2):
	return graph.edges[node1, node2]['var']

#=======================================================================
#Functions for compare path sets i.e messy integration

def prob_comparator(x, paths, node1, node2, cost1, cost2):
	paths1 = paths[node1]
	paths2 = paths[node2]

	#calculate f(x, node1)
	f_sub = 1.0
	for i in range(len(paths1)):
		mean_i = paths1[i].mean
		var_i = paths1[i].var
		f_sub *= 0.5 * erfc((x - mean_i - cost1) * 1.0 / math.sqrt(2 * var_i))

	f = 1 - f_sub

	#calculate g(x, node2)
	g = 0
	for j in range(len(paths2)):
		mean_j = paths2[j].mean
		var_j = paths2[j].var
		factor1 = math.exp(-1 * ((x - cost2 - mean_j)/math.sqrt(2 * var_j))**2)/(math.sqrt(2 * math.pi * var_j))
		factor2 = 1.0
		for i in range(len(paths2)):
			if i != j:
				mean_i = paths2[i].mean
				var_i = paths2[i].var
				factor2 *= 0.5 * erfc((x - mean_i - cost2)/math.sqrt(2 * var_i))

		g += factor1 * factor2

	return f * g


def comparePathSets(g, end, current, nbrs, paths):

	if len(nbrs) == 0:
		print("ERROR 0 neighbors")
		return None

	if len(nbrs) == 1:
		return nbrs[0]

	#should compare the concrete cost of the end with the expected cost of another path to it. for now just go to end - this will cause problems
	if end in nbrs:
		return end

	else:

		#find the actual costs of the edges
		costs = {}
		for nbr in nbrs:
			sample = np.random.normal(g.edges[current, nbr]['mean'], g.edges[current, nbr]['var'])
			if sample < 0:
				sample = 0

			costs[nbr] = sample

		min_node = nbrs[0]
		idx = 1
		while idx < len(nbrs):
			if costs[min_node] == 0:
				break
			elif costs[nbrs[idx]] == 0:
				min_node = nbrs[idx]
				break

			else:
				prob = scipy.integrate.quad(lambda x: prob_comparator(x, paths, min_node, nbrs[idx], costs[min_node], costs[nbrs[idx]]), -np.inf, np.inf)

				if prob[0] < 0.5:
					min_node = nbrs[idx]

				idx += 1

		return min_node



#=======================================================================



def ragspath(graph, start, end, threshold):

	#Phase 1 - pruning
	open_paths = []
	heapq.heapify(open_paths)
	closed = defaultdict(list)

	path0 = PathObject([start], 0, 0)
	heapq.heappush(open_paths, path0)

	while len(open_paths) > 0:
		cur = heapq.heappop(open_paths)
		
		cur_node = cur.path[-1]
		closed[cur_node].append(cur)

		if cur_node != end:

			for succ in graph.neighbors(cur_node):
				if succ not in cur.path:
					#path is acyclic

					newm = cur.mean + get_mean(graph, cur_node, succ)
					newv = cur.var + get_var(graph, cur_node, succ)
					new_path = copy.deepcopy(cur.path)
					new_path.append(succ)

					p = PathObject(new_path, newm, newv)

					if not betterPathToNode(closed, p):
						heapq.heappush(open_paths, p)

		if len(open_paths) > 0 and betterPathToGoal(closed, end, heapq.nsmallest(1, open_paths)[0]):
			break

	if len(closed[end]) == 0:
		print("No path to goal")
		return []

	#Make the graph
	G_nd, path_sets = constructSubgraph(graph, closed, end, start)

	final = []
	current = start

	while current != end:
		final.append(current)

		nbrs = list(set([p.path[1] for p in path_sets[current] if p.path[1] not in final]))

		next_node = comparePathSets(G_nd, end, current, nbrs, path_sets)
		current = next_node

	final.append(end)

	return final




for i in range(20):
	test = gen.GraphGenerator(10, 10, 15)
	graph, ns = test.gen_graph(20, 30, 30)
	nx.draw_networkx(graph, ns, edgelist = graph.edges)
	# plt.show()

	print(nx.astar_path(graph, 0, 19))
	print(ragspath(graph, 0, 19, 0.8))


