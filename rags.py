from collections import defaultdict
import copy
import heapq
import math
import time
import networkx as nx 
import numpy as np
from scipy.special import erfinv, erfc
import scipy.integrate
import matplotlib.pyplot as plt


threshold = 1.0

def print_graph_2colors(graph, ns, nd_edges, color1):
    print("printing")
    nx.draw_networkx_edges(graph, ns, edgelist=graph.edges, arrowsize=2, node_size=100)
    path = []
    # for i in range(1, len(final)):
    #     path.append((final[i - 1], final[i]))
    nx.draw_networkx_edges(graph, ns, edgelist=list(nd_edges), arrowsize=2, edge_color=color1, width=1.5, node_size=100)
    # nx.draw_networkx_edges(graph, ns, edgelist = path, edge_color = color2)
    plt.axis('off')
    plt.show(block=True)

def print_graph_3colors(graph, ns, nd_edges, color1, final, color2):
    print("printing")
    nx.draw_networkx_edges(graph, ns, edgelist=graph.edges, arrowsize=2, node_size=100)
    path = []
    for i in range(1, len(final)):
        path.append((final[i - 1], final[i]))
    nx.draw_networkx_edges(graph, ns, edgelist=list(nd_edges), arrowsize=2, edge_color=color1, width=1.5, node_size=100)
    nx.draw_networkx_edges(graph, ns, edgelist=path, edge_color=color2, arrowsize=2, node_size=100, width=2)
    plt.axis('off')
    plt.show(block=True)

class PathObject:
    def __init__(self, path, mean, var):
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



def dom(p1, p2):
    """
    Returns whether pathobject p1 dominates path object p2
    """
    rhs = p2.mean + math.sqrt(2 * (p1.var + p2.var)) * erfinv(1 - 2 * threshold)
    return p1.mean < rhs

class RAGS:

    def __init__(self, graph, vertex_locs, start, end, thresh):
        self.g = graph
        self.current_pos = start
        self.end = end
        global threshold
        threshold = thresh

        self.ver_locs = vertex_locs #the locations of the vertices - this is onlyfor drawing the graph

        #metrics
        self.paths_expanded = 0
        self.paths_considered = 0
        self.heap_ops = 0
        self.closed_accesses = 0

        #open and closed
        self.opened = []
        self.closed = defaultdict(list)

        #take_step structures
        self.nd_edges = []
        self.path_sets = {} 

        self.updated = 1

        #the path that is taken
        self.path = [self.current_pos]
        self.cost = 0
        self.mean_cost = 0
        self.var_cost = 0

        self.first_plan = True
        self.nd_edges = []


    def init_strs(self):
        self.opened = []
        heapq.heapify(self.opened)
        path0 = PathObject([self.end], 0, 0)
        heapq.heappush(self.opened, path0)
        self.heap_ops += 1

        self.closed = defaultdict(list)


    def betterPathToNode(self, p):
        self.closed_accesses += 1
        last = p.path[-1]
        for i in self.closed[last]:
            if dom(i, p):
                return True
        return False

    def betterPathToGoal(self, p):
        self.closed_accesses += 1
        for i in self.closed[self.current_pos]:
            if dom(i, p):
                return True
        return False

    def constructSubgraph(self):
        # R = nx.Graph()
        #contains a map of nodes to paths from them to the goal as well as the mean/var of those paths
        #if closed is like g in a*, path_sets is like h in a*
        path_sets = defaultdict(set)
        raw_paths = defaultdict(list)
        nd_edges = set([])

        # R.add_node(start)

        #if backwards -
        temp = copy.deepcopy(self.closed[self.current_pos])
        for obj in temp:
            path = copy.deepcopy(obj.path)
            path.reverse()
            obj.path = path

        #if forwards - 
        # temp = closed[goal]

        self.paths_considered += len(temp)
        for p in temp:
            curm = p.mean
            curv = p.var
            path = p.path


            for i in range(1, len(path)):
                if i == len(path) - 1:
                    if path[i:] not in raw_paths[path[i - 1]]:
                        path_sets[path[i - 1]].add(PathObject(path[i:], 0, 0.1))
                        raw_paths[path[i - 1]].append(path[i:])
                else:   
                    curm -= self.get_mean(path[i - 1], path[i])
                    curv -= self.get_var(path[i - 1], path[i])

                    if curv <= 0:
                        curv = 0.1
                        print("bad")

                    if path[i:] not in raw_paths[path[i - 1]]:
                        path_sets[path[i - 1]].add(PathObject(path[i:], curm, curv))
                        raw_paths[path[i - 1]].append(path[i:])
                # if path[i] not in R.nodes:
                    # R.add_node(path[i])

                # if (path[i - 1], path[i]) not in R.edges:
                nd_edges.add((path[i-1], path[i]))
                    # R.add_edge(path[i - 1], path[i], 
                    #   mean = self.get_mean(path[i - 1], path[i]), 
                    #   var = self.get_var(path[i - 1], path[i]))


        final_sets = {}
        for k,v in path_sets.items():
            final_sets[k] = list(v)

        self.path_sets = final_sets
        # print_graph_2colors(self.g, self.ver_locs, nd_edges, 'r')

    

    def get_successors(self, node):
        print(self.g.successors(node))
        return self.g.successors(node)

    def get_mean(self, node1, node2):
        return self.g.edges[node1, node2]['mean']

    def get_var(self, node1, node2):
        return self.g.edges[node1, node2]['var']

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
            self.printPathsets()
            print(self.path)
            print("ERROR 0 neighbors")
            return None

        if len(nbrs) == 1:

            self.mean_cost += self.get_mean(current, nbrs[0])
            self.var_cost += self.get_var(current, nbrs[0])
            self.cost += max(0, np.random.normal(self.get_mean(current, nbrs[0]), 
                self.get_var(current, nbrs[0])))
            return nbrs[0]

        #should compare the concrete cost of the end with the expected cost of another path to it. for now just go to end - this will cause problems
        if self.end in nbrs:
            self.mean_cost += self.get_mean(current, self.end)
            self.var_cost += self.get_var(current, self.end)
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

    def backward_phase1(self):
        self.init_strs()

        i = 0   

        while len(self.opened) > 0:
            cur = heapq.heappop(self.opened)
            self.heap_ops += 1
            
            cur_node = cur.path[-1]
            if not self.betterPathToNode(cur):
                self.closed[cur_node].append(cur)
            else:
                continue

            i += 1
            if cur_node != self.current_pos:

                for pred in self.g.predecessors(cur_node):
                    if pred not in cur.path:
                        #path is acyclic

                        newm = cur.mean + self.get_mean(pred, cur_node)
                        newv = cur.var + self.get_var(pred,  cur_node)
                        new_path = copy.deepcopy(cur.path)
                        new_path.append(pred)

                        p = PathObject(new_path, newm, newv)

                        if not self.betterPathToNode(p):
                            heapq.heappush(self.opened, p)
                            self.heap_ops += 1

            if len(self.opened) > 0 and self.betterPathToGoal(heapq.nsmallest(1, self.opened)[0]):
                break

        if len(self.closed[self.current_pos]) == 0:
            print("No path to goal")
            return []

        print("prune graph ",i)
        self.paths_expanded += i


    def takestep(self):

        if self.updated:
            self.constructSubgraph()
            self.updated = False

        if self.current_pos == self.end:
            print("reached goal")
            return None

        nbrs = list(set([p.path[0] for p in self.path_sets[self.current_pos] if p.path[0] not in self.path]))
        
        next_node = self.comparePathSets(self.current_pos, nbrs)

        self.current_pos = next_node
        self.path.append(next_node)

        return next_node



    def printPathsets(self):
        for k,v in self.path_sets.items():
            print(k)
            for path in v:
                print(path.path, ", ", path.mean, ", ", path.var)

        return

    def run(self, gg):

        time0 = time.time()
        x = self.backward_phase1()
        if x == []:
            return [-1 for i in range(8)]

        print(len(self.closed[self.current_pos]), "paths")
        time1 = time.time()

        #Make the graph

        node = self.takestep()
        # print("next path node ", node)

        update_ctr = 0

        while node != self.end:

            updates = []
            for edge in gg.changes1.keys():
                if edge[0] in self.g.successors(self.current_pos):
                    # print("updating bc of edge ", edge)
                    meanvar = gg.changes1[edge]
                    self.g.edges[edge[0], edge[1]]["mean"] = meanvar[0]
                    self.g.edges[edge[0], edge[1]]["var"] = meanvar[1]

                    updates.append(edge)

            for edge in updates:
                gg.changed_edge_r(edge)

            if len(updates) > 0:
                self.updated = True
                update_ctr += 1
                closed = self.backward_phase1()
            
            node = self.takestep()
            print("next path node ", node)


        # print_graph_3colors(self.g, self.ver_locs, self.nd_edges, 'r', self.path, 'b')

        end_time = time.time()
        total_time = end_time - time0
        replan_time = end_time - time1
        print("\n\n================")
        print("replan time: ", replan_time, "updates: ", update_ctr, "cost: ", self.cost)


        print("heap", self.heap_ops, "paths ", self.paths_expanded)
        return total_time, replan_time, self.cost, self.mean_cost, self.var_cost, update_ctr, self.paths_expanded, self.heap_ops, self.paths_considered, self.closed_accesses

