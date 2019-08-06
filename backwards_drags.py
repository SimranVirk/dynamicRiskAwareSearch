import networkx as nx 
import heapq
import numpy as np
from scipy.special import erfinv, erfc
import scipy.integrate
import copy
from collections import defaultdict
import generate_graphs as gen
# import matplotlib.pyplot as plt
import math
import time
import test

threshold = 0.75

# def print_graph_2colors(graph, ns, nd_edges, final, color1, color2):
#     print("printing")
#     nx.draw_networkx(graph, ns, edgelist = graph.edges)
#     path = []
#     for i in range(1, len(final)):
#         path.append((final[i - 1], final[i]))
#     nx.draw_networkx_edges(graph, ns, edgelist = list(nd_edges), edge_color = color1)
#     nx.draw_networkx_edges(graph, ns, edgelist = path, edge_color = color2)
#     plt.axis('off')
#     plt.show(block=True)

def dom(p1, p2):
	"""
	Returns whether pathobject p1 dominates path object p2
	"""
	rhs = p2.mean + math.sqrt(2 * (p1.var + p2.var)) * erfinv(1 - 2 * threshold)
	return p1.mean < rhs

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

	def __init__(self, graph, ns, start, end, threshold):
		self.g = graph
		self.start = start
		self.end = end
		self.threshold = threshold
		self.vert_locs = ns #for drawing the graph

		#open and closed
		self.opened = []
		heapq.heapify(self.opened)
		path0 = PathObject([end], 0, 0)
		heapq.heappush(self.opened, path0)

		self.closed = defaultdict(list)

		#take_step structures
		self.nd_edges = []
		self.path_sets = {}	

		#whether to re-generate above two structures
		self.updated = 1

		#the path that is taken
		self.current_pos = self.start
		self.path = [self.start]
		self.cost = 0


	def get_mean(self, node1, node2):
		return self.g.edges[node1, node2]['mean']

	def get_var(self, node1, node2):
		return self.g.edges[node1, node2]['var']

	def setDomPath(self, path, s):
		for x in s:
			if dom(x, path):
				return True
		return False


	def betterPathToNode(self, p):
		last = p.path[-1]
		return self.setDomPath(p, self.closed[last])

	def betterPathToGoal(self, p):
		#since it is backwards search, the goal is the current position
		return self.setDomPath(p, self.closed[self.current_pos])

	#==========================================================================
	#updation functions

	def deleteFromClosed(self, u, v):
		"""
		 Delete everythin in closed that contains edge u,v
		"""

		#*** this is making the assumption that we only get updated information at the
		#edges we touch

		removals = []
		for i in range(len(self.closed[u])):
			obj = self.closed[u][i]
			if obj.path[-2] == v:
				removals.append(obj)

		for obj in removals:
			self.closed[u].remove(obj)
		
		return
						

	def deleteFromOpen(self, u, v):
		"""
		Delete everything in opened that contains edge u,v
		"""
		#*** this removes things from list while iterating through it
		removals = []
		for obj in self.opened:
			if u in obj.path and v in obj.path:
				if obj.path.index(v) + 1 == obj.path.index(u):
					removals.append(obj)


		for obj in removals:
			self.opened.remove(obj)



	def getWorsePaths(self, p, s):
		out = []
		for obj in s:
			if dom(p, obj):
				out.append(obj)

		return out

	def update_edge(self, u, v, newm, newv):
		self.deleteFromClosed(u, v)
		self.deleteFromOpen(u, v)
		temp = self.closed[u]

		for p_v in self.closed[v]:
			path = copy.deepcopy(p_v.path)
			path.append(u)
			mean = p_v.mean + self.get_mean(u,v)
			var = p_v.var + self.get_var(u, v)
			p_u = PathObject(path, mean, var)

			
			if self.setDomPath(p_u, temp):
				#new path is bad
				continue
			else:
				#new path is good
				for p in self.getWorsePaths(p_u, temp):
					self.deleteFromClosed(p.path[-2], p.path[-1])
					self.deleteFromOpen(p.path[-2], p.path[-1])

				temp.append(p_u)
				heapq.heappush(self.opened, p_u)

		self.closed[u] = temp

		return 

	#=======================================================================================

	def getBetterPathstoNode(self, node, p):
		out = []
		for obj in self.closed[node]:
			if dom(p, obj):
				out.append(obj)

		return out

	def prune_graph(self):
		i = 0
		while len(self.opened) > 0 and not self.betterPathToGoal(heapq.nsmallest(1, self.opened)[0]):
			i += 1
			cur = heapq.heappop(self.opened)
		
			cur_node = cur.path[-1]

			self.closed[cur_node].append(cur)

			#Removing this check because of assumption that edge changes are
			#discovered when the edge is reached
			# for path in self.getBetterPathstoNode(cur_node, cur):
			# 	# print("bad")
			# 	#This should never run w/o memoization

			# 	#remove paths from closed and open
			# 	self.deleteFromClosed(path.path[-2], path.path[-1])
			# 	self.deleteFromOpen(path.path[-2], path.path[-1])

			if cur_node != self.current_pos:

				for pred in self.g.predecessors(cur_node):
					if pred not in cur.path:
						#path is acyclic

						newm = cur.mean + self.get_mean(pred, cur_node)
						newv = cur.var + self.get_var(pred, cur_node)
						new_path = copy.deepcopy(cur.path)
						new_path.append(pred)

						p = PathObject(new_path, newm, newv)

						if not self.betterPathToNode(p):
							#***check if this thing is not already in closed - do i need to do this?
							heapq.heappush(self.opened, p)


		print("prune graph ", i, " iterations")

		if len(self.closed[self.current_pos]) == 0:
			print("No path to goal")
			return []

	#===============================================================================================
	# Take step loop
	#===============================================================================================

	def printPathsets(self):
		for k,v in self.path_sets.items():
			print(k)
			for path in v:
				print(path.path, ", ", path.mean, ", ", path.var)

		return

	def getGoodPathSubsets(self):
		#contains a map of nodes to paths from them to the goal as well as the mean/var of those paths
		#if closed is like g in a*, path_sets is like h in a*
		path_sets = defaultdict(set)
		nd_edges = set([])

		temp = copy.deepcopy(self.closed[self.start])
		for obj in temp:
			path = copy.deepcopy(obj.path)
			path.reverse()
			obj.path = path

		for p in temp:
			curm = p.mean
			curv = p.var
			path = p.path

			for i in range(1, len(path)):
				if i == len(path) - 1:
					path_sets[path[i - 1]].add(PathObject(path[i:], 0, 0.1))
				else:
					curm -= self.g.edges[path[i-1], path[i]]['mean']
					curv -= self.g.edges[path[i-1], path[i]]['var']

					if curv <= 0:
						curv = 0.1
						print("bad")

					path_sets[path[i - 1]].add(PathObject(path[i:], curm, curv))

				nd_edges.add((path[i-1], path[i]))


		self.path_sets = {}
		for k,v in path_sets.items():
			self.path_sets[k] = list(v)

		self.nd_edges = nd_edges

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
			self.cost += max(0, np.random.normal(self.g.edges[current, nbrs[0]]["mean"], 
				self.g.edges[current, nbrs[0]]["var"]))
			return nbrs[0]

		#should compare the concrete cost of the end with the expected cost of another path to it. for now just go to end - this will cause problems
		if self.end in nbrs:
			self.cost += max(0, np.random.normal(self.g.edges[current, self.end]["mean"], 
				self.g.edges[current, self.end]["var"]))
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

			self.cost += costs[min_node]
			return min_node



	#=======================================================================

	def take_step(self):
		# time3 = time.time()
		#check if need to make path sets
		if self.updated:
			self.getGoodPathSubsets()
			self.updated = False

		if self.current_pos == self.end:
			print("reached goal")
			return None

		#all the next nodes in the paths starting at current pos, excluding the ones already in the
		#path taken
		nbrs = list(set([p.path[0] for p in self.path_sets[self.current_pos] if p.path[0] not in self.path]))

		next_node = self.comparePathSets(self.current_pos, nbrs)

		self.current_pos = next_node
		self.path.append(next_node)

		return next_node

	#=============================================================================
	def run(self, gg):
		self.prune_graph()
		time1 = time.time()

		node = self.take_step()
		print("next path node ", node)

		update_ctr = 0

		while node != self.end:
			updates = []

			for edge in gg.changes.keys():
				if edge[0] == node:
					updates.append(edge)


			if len(updates) > 0:
				update_ctr += 1

			for edge in updates:
				print("updating edge ", edge)
				meanvar = gg.changes[edge]
				self.g.edges[edge[0], edge[1]]["mean"] = meanvar[0]
				self.g.edges[edge[0], edge[1]]["var"] = meanvar[1]

				self.update_edge(edge[0], edge[1], meanvar[0], meanvar[1])
				self.updated = True

			if len(updates) > 0:
				self.prune_graph()

			node = self.take_step()
			#see if any nbr edges are different
			print("next path node ",node)


		total = time.time() - time1
		print("\n\n================")
		print("time: ", total, "updates: ", update_ctr, "cost: ", self.cost)

		return total, update_ctr, self.cost


	def run_replan(self, gg):
		lst = test.run(self.g, self.start, self.end, gg)
		print("done")
		return lst

		



if __name__ == "__main__":

	t1 = 0
	u1 = 0
	t2 = 0
	u2 = 0
	c1 = 0
	c2 = 0
	# for i in range(10):

	# nx.draw_networkx(graph, ns, edgelist = graph.edges)
	# plt.show()

	better = 0
	bettercost  = 0
	for i in range(10):
		gg = gen.GraphGenerator(10, 10, 20)
		graph, ns = gg.gen_graph(30, 40, 40)	
		d = DRAGS(graph, ns, 0, 29, 0.75)
		gg.getEdgeChanges(30)

		l1 = d.run(gg)
		l2 = d.run_replan(gg)

		t1 += l1[0]
		u1 += l1[1]

		t2 += l2[0]
		u2 += l2[1]

		c1 += l1[2]
		c2 += l2[2]

		if l1[0] < l2[0]:
			better += 1

		if l1[2] < l2[2]:
			bettercost += 1

	print("============")
	print("stats - t1: ", t1,", u1: ",u1, "c1 ", c1)
	print("t2: ", t2, "u2: ", u2, "c2: ", c2)
	print("better ", better)
	print("better cost ", bettercost)

	# tuple1 = d.take_step()




	#plan backwards -- from goal to start

	#edge updates are found by following the path and discovering that the input
	# was incorrect - this models seeing a dynamic/uncertain environment

	#in info gain - would be a bunch of edges nearby the edge just traversed (i.e a bunch of 
	#edges that a correlated with the edge just traversed)

	#want to compare a certain amount of updates with replanning on updates
	#want faster compute/comparable paths




