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
import time
import test
import csv

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

def testdom(m1, v1, m2, v2):
	rhs = m2 + math.sqrt(2 * (v1 + v2)) * erfinv(1 - 2 * threshold)
	return m1 < rhs

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

		#metrics
		self.path_expansions = 0
		self.paths_considered = 0
		self.heap_ops = 0
		self.closed_accesses = 0

		#open and closed
		self.opened = []
		heapq.heapify(self.opened)
		path0 = PathObject([end], 0, 0)
		heapq.heappush(self.opened, path0)
		self.heap_ops += 1

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
		self.mean_cost = 0
		self.var_cost = 0

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


		#rmove from the goal
		removals = []
		for p in self.closed[self.current_pos]:
			if u in p.path and v in p.path:
				if p.path.index(u) == p.path.index(v) + 1:
					removals.append(p)
		for p in removals:
			self.closed[self.current_pos].remove(p)

		#remove from u
		removals = []
		for p in self.closed[u]:
			if len(p.path) < 2:
				continue
			if p.path[-2] == v:
				removals.append(p)

		for p in removals:
			self.closed[u].remove(p)			

		# #remove from nbrs
		nodes = list(self.g.predecessors(u))
		seen = []
		i = 0
		while len(nodes) > 0:
			i += 1

			cur = nodes.pop(0)
			seen.append(cur)
			removals = []


			for p in self.closed[cur]:
				if u in p.path and v in p.path:
					if p.path.index(u) == p.path.index(v) + 1:
						removals.append(p)

			for p in removals:
				self.closed[cur].remove(p)						

			if len(removals) > 0:
				for node in self.g.predecessors(cur):
					if node not in seen:
						nodes.append(node)

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

	
	def update_edge(self, u, v):
		# print("updating", u, v)

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
					self.deleteFromClosed(p.path[-1], p.path[-2])
					self.deleteFromOpen(p.path[-1], p.path[-2])

				temp.append(p_u)
				heapq.heappush(self.opened, p_u)
				self.heap_ops += 1

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
		opened = self.opened
		current_pos = self.current_pos
		closed = self.closed

		i = 0
		while len(opened) > 0: 
			cur = heapq.heappop(opened)		
			self.heap_ops += 1

			cur_node = cur.path[-1]

			if not self.betterPathToNode(cur):
				closed[cur_node].append(cur)
			else:
				continue

			i += 1
			if cur_node != current_pos:

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
							heapq.heappush(opened, p)
							self.heap_ops += 1

			if len(opened) > 0 and self.betterPathToGoal(heapq.nsmallest(1, opened)[0]):
				break

		print("prune graph ", i, " iterations")
		self.path_expansions += i

		if len(closed[current_pos]) == 0:
			print("No path to goal")
			return []

		self.closed = closed
		self.opened = opened
		# print(len(self.closed[self.current_pos]), "paths")


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

		#***sometimes this messes up and has multiple of the same path
		temp = copy.deepcopy(list(set(self.closed[self.current_pos])))
		self.paths_considered += len(temp)

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

					curm -= self.get_mean(path[i - 1], path[i])
					curv -= self.get_var(path[i - 1], path[i])


					if curv <= 0:
						curv = 0.1
						print("bad")
						print("path ", p.path, "edge for next iter (should be last edge)", path[i], path[i + 1])
						print(self.g.edges[path[i], path[i + 1]]['var'])

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
			print("ERROR 0 neighbors for node ", current)
			return None

		if len(nbrs) == 1:
			self.mean_cost += self.g.edges[current, nbrs[0]]["mean"] 
			self.var_cost += self.g.edges[current, nbrs[0]]["var"]
			self.cost += max(0, np.random.normal(self.g.edges[current, nbrs[0]]["mean"], 
				self.g.edges[current, nbrs[0]]["var"]))
			return nbrs[0]

		#should compare the concrete cost of the end with the expected cost of another path to it. for now just go to end - this will cause problems
		if self.end in nbrs:
			self.mean_cost += self.g.edges[current, self.end]["mean"] 
			self.var_cost += self.g.edges[current, self.end]["var"]
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

			self.mean_cost += self.g.edges[current, min_node]["mean"]
			self.var_cost += self.g.edges[current, min_node]["var"]
			self.cost += costs[min_node]
			return min_node



	#=======================================================================

	def take_step(self):
		#check if need to make path sets
		if self.updated:
			self.getGoodPathSubsets()
			self.updated = False

		if self.current_pos == self.end:
			print("reached goal")
			return None

		#all the next nodes in the paths starting at current pos, excluding the ones already in the
		#path taken

		# self.printPathsets()
		nbrs = list(set([p.path[0] for p in self.path_sets[self.current_pos] if p.path[0] not in self.path]))

		next_node = self.comparePathSets(self.current_pos, nbrs)

		self.current_pos = next_node
		self.path.append(next_node)

		return next_node

	#=============================================================================
	def run(self, gg, ):
		time0 = time.time()
		x = self.prune_graph()
		memo =len(self.closed[self.current_pos]) 
		print(len(self.closed[self.current_pos]), "paths")
		if x == []:
			return [-1 for i in range(6)]
		time1 = time.time()

		for path in self.closed[self.current_pos]:
			print(path.path)



		node = self.take_step()
		# print("next path node ", node)

		update_ctr = 0

		while node != self.end:
			updates = []

			#2 step lookahead
			for edge in gg.changes.keys():
				if edge[0] in self.g.successors(node):
					updates.append(edge)

					meanvar = gg.changes[edge]
					self.g.edges[edge[0], edge[1]]["mean"] = meanvar[0]
					self.g.edges[edge[0], edge[1]]["var"] = meanvar[1]


			if len(updates) > 0:
				update_ctr += 1
				self.updated = True


			print(len(updates), "updates")
			for edge in updates:
				gg.changedEdge(edge)
				self.update_edge(edge[0], edge[1])


			if len(updates) > 0:
				self.prune_graph()

			node = self.take_step()
			#see if any nbr edges are different
			# print("next path node ",node)


		time_end = time.time()
		total_time = time.time() - time0
		replan_time = time.time() - time1

		print("\n\n================")
		print("replan time: ", replan_time, "updates: ", update_ctr, "cost: ", self.cost)


		print("heap", self.heap_ops, "paths ", self.path_expansions)

		return total_time, replan_time, self.cost, self.mean_cost, self.var_cost, update_ctr, self.path_expansions, self.heap_ops, self.paths_considered, memo


	def run_replan(self, gg, graph):
		rags = test.RAGS(graph, self.vert_locs, self.start, self.end)
		lst = rags.run(gg)
		print("done")
		return lst

		



if __name__ == "__main__":

	num_nodes = 30
	total_num_edges = 0
	num_changes = 5
	iters = 20
	start_goal_disconnected = 0


	# 1 -- dynamic rags, 2 -- replan from scratch
	#total time, replan time, number of updates, cost, mean and variance of cost
	t1 = 0
	rt1 = 0
	u1 = 0
	c1 = 0
	mc1 = 0
	vc1 = 0

	t2 = 0
	rt2 = 0
	u2 = 0
	c2 = 0
	mc2 = 0
	vc2 = 0

	better = 0
	bettercost  = 0

	with open('results/updated/final_%d_%d.csv' % (num_nodes, num_changes), 'a', newline = '') as csvfile:

		fieldnames = ['num_nodes', 'num_edges', 'threshold', 'num_edges_changed', 
	'd_totaltime', 'd_replantime', 'd_cost', 'd_meancost', 'd_varcost', 'd_numupdates',
	'd_paths_expanded', 'd_heap_ops', 'd_paths_considered',
	'r_totaltime', 'r_replantime', 'r_cost', 'r_meancost', 'r_varcost', 'r_numupdates',
	'r_paths_expanded', 'r_heap_ops', 'r_paths_considered']
		writer = csv.DictWriter(csvfile, fieldnames = fieldnames)

		writer.writeheader()



		for i in range(iters):
			print("iteration ", i)
			gg = gen.GraphGenerator(10, 10)
			graph = gg.load_graph(num_nodes, i)
			ns = []
			# nx.draw_networkx_edges(graph, ns, edgelist = graph.edges, arrowsize = 2)
			# plt.axis('off')
			# plt.show()


			num_edges = len(graph.edges)
			total_num_edges += num_edges
			print(num_edges, " edges")

			graph_copy = copy.deepcopy(graph)
			d = DRAGS(graph_copy, ns, 0, num_nodes - 1, 0.75)
			gg.getEdgeChanges(int(num_edges / num_changes))
			print(num_edges/ num_changes, " edge changes")

			l2 = d.run_replan(gg, graph)
			if l2[0] == -1:
				start_goal_disconnected += 1
				continue

			l1 = d.run(gg)

			t1 += l1[0]
			rt1 += l1[1]
			u1 += l1[5]

			t2 += l2[0]
			rt2 += l2[1]
			u2 += l2[5]

			c1 += l1[2]
			c2 += l2[2]


			writer.writerow({'num_nodes': num_nodes, 'num_edges' : total_num_edges, 
					'threshold': threshold, 'num_edges_changed': num_changes, 
					'd_totaltime': l1[0], 'd_replantime': l1[1], 'd_cost': l1[2], 'd_meancost': l1[3], 'd_varcost': l1[4], 'd_numupdates': l1[5],
			'd_paths_expanded': l1[6], 'd_heap_ops': l1[7], 'd_paths_considered': l1[8],
			'r_totaltime': l2[0], 'r_replantime': l2[1], 'r_cost': l2[2], 'r_meancost': l2[3], 'r_varcost': l2[4], 'r_numupdates': l2[5],
			'r_paths_expanded': l2[6], 'r_heap_ops': l2[7], 'r_paths_considered': l2[8]})



			if l1[1] < l2[1]:
				better += 1

			if l1[3] < l2[3]:
				bettercost += 1

	print("============")
	print("stats - rt1: ", rt1, "t1", t1 , "u1: ",u1, "c1 ", c1)
	print("rt2: ", rt2, "t2", t2, "u2: ", u2, "c2: ", c2)
	print("better ", better)
	print("better cost ", bettercost)

	# tuple1 = d.take_step()

	


