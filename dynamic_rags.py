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


def dom(p1, p2):
	"""
	Returns whether pathobject p1 dominates path object p2
	"""
	rhs = p2.mean + math.sqrt(2 * (p1.var + p2.var)) * erfinv(1 - 2 * threshold)
	return p1.mean < rhs

#intelligent updates
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


class DRAGS:

	def __init__(self, graph, start, end, threshold):
		self.g = graph
		self.start = start
		self.end = end
		self.threshold = threshold

		self.opened = []
		heapq.heapify(self.opened)

		self.closed = defaultdict(list)
	
		path0 = PathObject([start], 0, 0)
		heapq.heappush(opened, path0)

		self.updated = 1
		self.path_sets = {}


	def betterPathToNode(self, p):
	last = p.path[-1]
	for i in self.closed[last]:
		if dom(i, p):
			return True
	return False

	def betterPathToGoal(self, p):
		for i in self.closed[self.end]:
			if dom(i, p):
				return True
		return False

	def deleteFromClosed(self, u, v):
		"""
		 Delete everythin in closed that contains edge u,v

		This needs to be optimized

		"""
		q = [v]
		
		while len(q) > 0:
			cur = q[0]
			for obj in closed[cur]:
				if u in obj.path and v in obj.path and obj.path.index(u) + 1 == obj.path.index(v):
					closed[obj.path[-1]].remove(obj)
					for succ in self.g.successors(obj.path[-1]):
						if succ not in q:
							q.append(succ)	
						

	def deleteFromOpen(self, u, v):
		"""
		Delete everything in opened that contains edge u,v
		"""
		for obj in self.opened:
			if u in obj.path and v in obj.path and obj.path.index(u) + 1 == obj.path.index(v):
				self.opened.remove(obj)

	def update_edge(self, u, v, newm, newv):
		pass


	def getBetterPathstoNode(self, node, p):
		out = []
		for obj in self.closed[node]:
			if dom(p, obj):
				out.append(obj)

		return out

	def prune_graph(self):
		while len(self.opened) > 0:
			cur = heapq.heappop(self.opened)
		
			cur_node = cur.path[-1]

			self.closed[cur_node].append(cur)
			for path in self.getBetterPathstoNode(cur_node, cur):
				#This should never run w/o memoization

				#remove paths from closed and open
				self.deleteFromClosed(path.path[-2], path.path[-1])
				self.deleteFromOpened(path.path[-2], path.path[-1])

			if cur_node != end:

				for succ in self.g.neighbors(cur_node):
					if succ not in cur.path:
						#path is acyclic

						newm = cur.mean + get_mean(self.g, cur_node, succ)
						newv = cur.var + get_var(self.g, cur_node, succ)
						new_path = copy.deepcopy(cur.path)
						new_path.append(succ)

						p = PathObject(new_path, newm, newv)

						if not self.betterPathToNode(p):
							#***check if this thing is not already in closed - do i need to do this?
							heapq.heappush(self.opened, p)

			if len(self.opened) > 0 and self.betterPathToGoal(heapq.nsmallest(1, open_paths)[0]):
				break

		if len(closed[end]) == 0:
			print("No path to goal")
			return []

	#===============================================================================================
	# Take step loop
	#===============================================================================================

	def constructSubgraph(self):
		#contains a map of nodes to paths from them to the goal as well as the mean/var of those paths
		#if closed is like g in a*, path_sets is like h in a*
		path_sets = defaultdict(set)
		nd_edges = set([])

		for p in self.closed[self.end]:
			curm = p.mean
			curv = p.var
			path = p.path

			path_sets[path[0]].add(p)

			for i in range(1, len(path)):
				if i < len(path) - 1:
					curm -= self.g.edges[path[i-1], path[i]]['mean']
					curv -= self.g.edges[path[i-1], path[i]]['var']

					path_sets[path[i]].add(PathObject(path[i:], curm, curv))

				nd_edges.add((path[i-1], path[i]))


		self.path_sets = {}
		for k,v in path_sets.items():
			self.path_sets[k] = list(v)

		return 

	#=======================================================================
	#Functions for compare path sets i.e messy integration

	def prob_comparator(self, x, node1, node2, cost1, cost2):
		paths1 = self.path_sets[node1]
		paths2 = self.path_sets[node2]

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


	def comparePathSets(self, current, nbrs):

		if len(nbrs) == 0:
			print("ERROR 0 neighbors")
			return None

		if len(nbrs) == 1:
			return nbrs[0]

		#should compare the concrete cost of the end with the expected cost of another path to it. for now just go to end - this will cause problems
		if self.end in nbrs:
			return self.end

		else:

			#find the actual costs of the edges
			costs = {}
			for nbr in nbrs:
				sample = np.random.normal(self.g.edges[current, nbr]['mean'], self.g.edges[current, nbr]['var'])
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
					prob = scipy.integrate.quad(lambda x: self.prob_comparator(x, min_node, nbrs[idx], costs[min_node], costs[nbrs[idx]]), -np.inf, np.inf)

					if prob[0] < 0.5:
						min_node = nbrs[idx]

					idx += 1

			return min_node



	#=======================================================================

	def take_step(self):
		#check if need to make path sets
		if self.updated:
			self.getGoodPathSubsets()
			self.updated = False

		final = []
		current = self.start

		while current != self.end:
			final.append(current)
			# takestep()

			nbrs = list(set([p.path[1] for p in self.path_sets[current] if p.path[1] not in final]))

			next_node = self.comparePathSets(current, nbrs)
			current = next_node

		final.append(self.end)

		# print_graph_2colors(graph, ns, nd_edges, final, 'r', 'b')

		return final










