from collections import defaultdict
import copy
import heapq
import math
import time
import matplotlib.pyplot as plt
import numpy as np
import rags
from scipy.special import erfinv, erfc
import scipy.integrate


threshold = 1.0

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

    #***
    #def an equals method


class DRAGS:

    def __init__(self, graph, ns, start, end, thresh):
        self.g = graph
        self.start = start
        self.end = end
        global threshold
        threshold = thresh
        self.vert_locs = ns #the locations of the vertices - this is onlyfor drawing the graph

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

    def get_closed(self, node):
        self.closed_accesses += 1
        return self.closed[node]

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
        return self.setDomPath(p, self.get_closed(last))

    def betterPathToGoal(self, p):
        #since it is backwards search, the goal is the current position
        return self.setDomPath(p, self.get_closed(self.current_pos))

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
        current_closed = self.get_closed(self.current_pos)
        for p in current_closed:
            if u in p.path and v in p.path:
                if p.path.index(u) == p.path.index(v) + 1:
                    removals.append(p)
        for p in removals:
            self.closed[self.current_pos].remove(p)

        #remove from u
        removals = []
        u_closed = self.get_closed(u)
        for p in u_closed:
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

            cur_closed = self.get_closed(cur)
            for p in cur_closed:
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
            self.heap_ops += 1
            self.opened.remove(obj)


    def getWorsePaths(self, p, s):
        out = []
        for obj in s:
            if dom(p, obj):
                out.append(obj)

        return out

    
    def update_edge(self, u, v):
        initial = [p.path for p in self.closed[u]]
        self.deleteFromClosed(u, v)
        self.deleteFromOpen(u, v)


        temp = self.get_closed(u)
        removed = False

        v_closed = self.get_closed(v)
        for p_v in v_closed:
            path = copy.deepcopy(p_v.path)
            path.append(u)
            mean = p_v.mean + self.get_mean(u,v)
            var = p_v.var + self.get_var(u, v)
            p_u = PathObject(path, mean, var)

            
            if self.setDomPath(p_u, temp):
                #new path is bad
                if not removed and path in initial:
                    removed = True
                continue
            else:
                #new path is good
                for p in self.getWorsePaths(p_u, temp):
                    self.deleteFromClosed(p.path[-1], p.path[-2])
                    self.deleteFromOpen(p.path[-1], p.path[-2])

                temp.append(p_u)
                heapq.heappush(self.opened, p_u)
                self.heap_ops += 1


        #see if any paths need to be added
        paths_u = [p.path for p in temp]
        potential_additions = []
        heapq.heapify(potential_additions)

        if removed:
            for node in self.g.successors(u):
                closed_node = self.get_closed(node)
                for path in closed_node:
                    check = copy.deepcopy(path.path)
                    check.append(u)
                    if check not in paths_u:
                        mean = path.mean + self.get_mean(u, node)
                        var = path.var + self.get_var(u, node)
                        p_u = PathObject(check, mean, var)

                        heapq.heappush(potential_additions, p_u)

        while len(potential_additions) > 0:
            p_next = heapq.heappop(potential_additions)
            if self.setDomPath(p_next, temp):
                break
            else:
                temp.append(p_next)
                heapq.heappush(self.opened, p_next)


        self.closed[u] = temp




        return 

    #=======================================================================================

    def getBetterPathstoNode(self, node, p):
        out = []
        node_closed = self.get_closed(node)
        for obj in node_closed:
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
        print("print path sets")
        for k,v in self.path_sets.items():
            print(k)
            for path in v:
                print(path.path, ", ", path.mean, ", ", path.var)

        return

    def getGoodPathSubsets(self):
        path_sets = defaultdict(set)
        raw_paths = defaultdict(list)
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
                    if path[i:] not in raw_paths[path[i - 1]]:
                        path_sets[path[i - 1]].add(PathObject(path[i:], 0, 0.1))
                        raw_paths[path[i - 1]].append(path[i:])
                else:

                    curm -= self.get_mean(path[i - 1], path[i])
                    curv -= self.get_var(path[i - 1], path[i])


                    if curv <= 0:
                        curv = 0.1
                        print("bad")
                        print("path ", p.path, "edge for next iter (should be last edge)", path[i], path[i + 1])
                        print(self.g.edges[path[i], path[i + 1]]['var'])

                    if path[i:] not in raw_paths[path[i - 1]]:
                        path_sets[path[i - 1]].add(PathObject(path[i:], curm, curv))
                        raw_paths[path[i - 1]].append(path[i:])

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
            self.printPathsets()
            print(self.path)
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
    def run(self, gg):
        time0 = time.time()
        x = self.prune_graph()
        print(len(self.closed[self.current_pos]), "paths")
        if x == []:
            return [-1 for i in range(6)]
        time1 = time.time()

        node = self.take_step()
        print("next path node ", node)

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
                gg.changed_edge(edge)
                self.update_edge(edge[0], edge[1])


            if len(updates) > 0:
                self.prune_graph()

            node = self.take_step()
            #see if any nbr edges are different
            print("next path node ",node)


        time_end = time.time()
        total_time = time.time() - time0
        replan_time = time.time() - time1

        print("\n\n================")
        print("replan time: ", replan_time, "updates: ", update_ctr, "cost: ", self.cost)


        print("heap", self.heap_ops, "paths ", self.path_expansions)

        return total_time, replan_time, self.cost, self.mean_cost, self.var_cost, update_ctr, self.path_expansions, self.heap_ops, self.paths_considered, self.closed_accesses


    def run_replan(self, gg, graph):

        r = rags.RAGS(graph, self.vert_locs, self.start, self.end, threshold)
        lst = r.run(gg)
        print("done")
        return lst

        



