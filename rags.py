# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 15:10:42 2019

@author: simra

Given a graph with a start and an exits and prob dist, use
RAGS to find the best path to the exit

Requires - a graph with each edge labelled with a mean and a variance
variance cannot be 0 -- causes div by zero error

Main function to call to find path -- find_path

things to add
-dynamic normal distributions for each edge -- more interesting -- gaussian process

multi agent - sub swarms
-exploration - ?
-decentralized planning and consolidation -- more of a demo


"""

from igraph import *
import networkx as nx
from scipy.special import erfinv, erfc
import scipy.integrate
import math
from math import sqrt
import matplotlib.pyplot as plt
import copy
from collections import defaultdict
import numpy as np
from heapq import heappush, heappop, heapify, nsmallest


global D_thresh

# =============================================================================
# =============================================================================
# # Test Graphs
# =============================================================================
# =============================================================================

G = nx.Graph()
G.add_nodes_from([0,1,2,3,4,5])
G.add_edges_from([(0,1), (0,2), (0,3), (1,2), (3,5),(2,4), (4,5)])
nx.set_edge_attributes(G, {(0,1): {'mean': 2, 'var':1},
                           (0,2): {'mean': 1, 'var':4},
                           (0,3): {'mean': 3, 'var':2},
                           (1,2): {'mean': 1, 'var':5},
                           (3,5): {'mean': 1, 'var':7},
                           (2,4): {'mean': 2, 'var':1},
                           (4,5): {'mean': 1, 'var':2}})



#Constructing an igraph
G1 = Graph()
G1.add_vertices(range(6))
edges = [(0,1), (0,2), (0,3), (1,2), (3,5),(2,4), (4,5)]
G1.add_edges(edges)
meansg1 = [2,1,3,1,1,2,1]
varsg1 = [1,1,2,1,1,2,1]
for i,e in enumerate(G1.es):
    e["mean"] = meansg1[i]
    e["var"] = varsg1[i] 


T = nx.Graph()
T.add_nodes_from(range(7))
T.add_edges_from([(0,1), (0,3), (0,2), (2,1), (1,4), (1,5), (4,7),(5,7), (3,6), (6,7)])
# nx.set_edge_attributes(T, {(0,1): {'mean': 5, 'var':5},
#                            (0,2): {'mean': 2, 'var':1},
#                            (0,3): {'mean': 4, 'var':4},
#                            (3,6): {'mean': 5, 'var':5},
#                            (6,7): {'mean': 7, 'var':0.1},
#                            (1,4): {'mean': 5, 'var':5},
#                            (4,7): {'mean': 6, 'var':1},
#                            (2,1): {'mean': 3, 'var':1},
#                            (1,5): {'mean': 8, 'var':8},
#                            (5,7): {'mean': 5, 'var':2}})
nx.set_edge_attributes(T, {(0,1): {'mean': 5, 'var':14},
                           (0,2): {'mean': 2, 'var':1},
                           (0,3): {'mean': 4, 'var':6},
                           (3,6): {'mean': 5, 'var':5},
                           (6,7): {'mean': 7, 'var':0.1},
                           (1,4): {'mean': 5, 'var':8},
                           (4,7): {'mean': 6, 'var':1},
                           (2,1): {'mean': 3, 'var':1},
                           (1,5): {'mean': 8, 'var':10},
                           (5,7): {'mean': 5, 'var':2}})


S = nx.Graph()
S.add_nodes_from(range(6))
S.add_edges_from([(0,1), (1,2), (2,5), (0,3), (3,5), (0,4), (4,5)])
nx.set_edge_attributes(S, {(0,1): {'mean': 4, 'var': 1},
                            (1,2): {'mean': 3, 'var': 1},
                            (2,5): {'mean': 5, 'var': 3},
                            (0,3): {'mean': 8, 'var': 3},
                            (3, 5): {'mean':6, 'var': 2},
                            (0,4): {'mean': 8, 'var': 1},
                            (4,5): {'mean': 9, 'var':1}})




# =============================================================================
# #Testing note -- 
# 
# Just keeping a record of sepcific numbers so I don't have to do these calculations again
# This graph tests phase 1 mainly. all paths have low variance so it is not very useful for testing
# phase 2 which is more probabilistic
# 
# The above graph (T) tests the following things
# - The open heap retruns the correct smallest element -- print each path as it is popped and manually verify
# - If we have seen a shorter(dominating) path to a current nbr, then the new path to it is not put on open
#     - here, we see the path 0 - 1 - 2 before 0 - 2. The new 0 - 2 path should trigger the non trivial case of 
#     nonDombySeen and return that it is in fact dominated by a seen path and should not be added to open queue
#     - if we change the weight of the 0 - 2 edge to 1.5 then we get the opposite result from nonDombySeen and it
#     is placed in the open queue
# -Ending prematurely
#     - with threshodl 0.65, the sumtraction factor in dominance calc for the first two paths seen to the end 
#     (with means 12, 14 and variances 3,3) is -1.62. in this case the lighter path dominates and there is only one
#     path in the non dominated portion
#     
#     - increasing threshold to 0.85 makes the subtraction factor -4.38 so the second path (14, 3) is included but a
#     third path with (18,3) is dominated by 12 and so we end the while loop prematurely.
# =============================================================================


H = nx.Graph()
H.add_nodes_from([0,1,2,3,4,5])
H.add_edges_from([(0,1), (0,2), (0,3), (1,5), (2,5),(3,4), (4,5)])
nx.set_edge_attributes(H, {(0,1): {'mean': 3, 'var':0.5},
                                   (1,5): {'mean': 3, 'var':1},
                                   (0,2): {'mean': 3, 'var':1},
                                   (2,5): {'mean': 5, 'var':1},
                                   (0,3): {'mean': 3, 'var':3},
                                   (3,4): {'mean': 2, 'var':10},
                                   (4,5): {'mean': 2, 'var':1}})


# =============================================================================
# #some graph printing functions
# =============================================================================
def print_graph(graph):
    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos)
    nx.draw_networkx_edges(graph, pos)
    nx.draw_networkx_labels(graph, pos)
    edge_labels = nx.get_edge_attributes(graph,'mean')
    nx.draw_networkx_edge_labels(graph, pos, labels = edge_labels)
    plt.axis('off')
    plt.show(block=True)  

# print_graph(S)

def print_graph_color(graph, nd_edges, color):
    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos)
    nx.draw_networkx_edges(graph, pos)
    nx.draw_networkx_labels(graph, pos)
    nx.draw_networkx_edges(graph, pos, edgelist = list(nd_edges), edge_color = color)
    plt.axis('off')
    plt.show(block=True) 


def print_graph_2colors(graph, nd_edges, final, color1, color2):
    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos)
    nx.draw_networkx_edges(graph, pos)
    nx.draw_networkx_labels(graph, pos)   
    path = []
    for i in range(1, len(final)):
        path.append((final[i - 1], final[i]))
    nx.draw_networkx_edges(graph, pos, edgelist = list(nd_edges), edge_color = color1)
    nx.draw_networkx_edges(graph, pos, edgelist = path, edge_color = color2)
    plt.axis('off')
    plt.show(block=True)

def print_graph_horizon(graph, nd_edges, final, subset, color1, color2, color3):
    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos)
    nx.draw_networkx_edges(graph, pos)
    nx.draw_networkx_labels(graph, pos)   
    path = []
    for i in range(1, len(final)):
        path.append((final[i - 1], final[i]))
    nx.draw_networkx_edges(graph, pos, edgelist = list(nd_edges), edge_color = color1)
    
    rest = []
    for i in range(1, len(subset)):
        rest.append((subset[i - 1], subset[i]))    

    nx.draw_networkx_edges(graph, pos, edgelist = path, edge_color = color2)
    nx.draw_networkx_edges(graph, pos, edgelist = rest, edge_color = color3)

    plt.axis('off')
    plt.show(block=True)    

    
"""
Custom object to put in a heap with a custom comparator function (__lt__)
This allows me to make the comparator based on Dominance as described in paper
"""
class MyHeapObj:
    def __init__(self, path, path_mean, path_var):
        self.path = path
        self.path_mean = path_mean
        self.path_var = path_var
        
    def __lt__(self, other):
#        print("hello from __lt__!")
        boolean = nonDom(self.path_mean, self.path_var, other.path_mean, other.path_var)
        
        #tie breaking of neither dominates the other ***
        boolean2 = nonDom(other.path_mean, other.path_var, self.path_mean, self.path_var)
        if boolean == boolean2:
            return self.path_mean < other.path_mean
        
        return not boolean
 
    

"""
Inputs: mean_a, var_a - mean and variance of path a
        mean_b, var_b - mean and variance of path b
Returns: True - a does not dominate b; else False
"""
def nonDom(mean_a, var_a, mean_b, var_b):
    rhs = mean_b + math.sqrt(2 * (var_a + var_b)) * erfinv(1 - 2 * D_thresh)
    if mean_a < rhs:
        return False
    else:
        return True #b is not dominated by a
    
               
# ***make a change here? -- since we are modelling targets, risk may be better associated with a node than an edge

# =============================================================================
# These get the mean and the variance of an edge
# =============================================================================
def get_mean(graph, node1, node2):
    return graph.edges[node1, node2]['mean']

def get_var(graph, node1, node2):
    return graph.edges[node1, node2]['var'] 
        

"""
makes sure that the new path under consideration is not dominated by
a path seen before that ends with the same node

Inputs: path - the path under consideration
        end_node- the last node in the path
        closed - a hashmap of nodes to paths to them (along with the mean and variance of the path)
Returns: True if there is not a path to end_node in closed that dominates path; False otherwise
"""
def nonDombySeen(path, end_node, closed):
    if end_node not in closed:
        return True
    else:
        #check that no seen path to this dominates this path
        # print("non trivial")
        nonDomflag = True
        for cpath in closed[end_node]:
            if not nonDom(cpath.path_mean, cpath.path_var, path.path_mean, path.path_var):
                nonDomflag = False
                break
        # print("no dom is")
        # print(nonDomflag)
        return nonDomflag
 
"""
makes sure that the path is not dominated by any of the paths in the path set

Inputs: path - the path to check against path set
        path_set - set of seen paths
        
Returns: True if there is a path in path_set that dominates path; False otherwise
"""
def DombyPaths(path, path_set):
    for nd_path in path_set:
        if not nonDom(nd_path.path_mean, nd_path.path_var, path.path_mean, path.path_var):
            # print("going to break while - ")
            # print(path.path)
            # print(nd_path.path_mean, nd_path.path_var, path.path_mean, path.path_var)
            return True
    return False

"""
Inputs - a graph, and the set of non dominated edges in that graph (i.e, paths that have
a high probability of giving you a low cost)

Returns -  a graph of the non dominated portion and the path sets - i.e a hashmap
of nodes to all the non dominated paths to the goal
(and the set of non dominated edges -- this is for printing the original graph with the non dominated portion highlighted)
"""
def getNdPortion(graph, p_nd):
    R = nx.Graph()
    path_sets = defaultdict(list)
    nd_edges = set([])
    for path_obj in p_nd:
        path = path_obj.path
        mean = path_obj.path_mean
        var = path_obj.path_var
        
        if path[0] not in R.nodes:
            R.add_node(path[0])
            
        new_mean = mean 
        new_var = var
        path_sets[path[0]].append(MyHeapObj(path, new_mean, new_var))    
            
        for i in range(1, len(path)):
            
            #add the rest of the path to the end to the path sets
            if i < len(path) - 1:
                new_mean = new_mean - graph.edges[path[i - 1], path[i]]['mean']
                new_var = new_var - graph.edges[path[i - 1], path[i]]['var']
                
                ## inefficient ??
                if path[i:] not in [x.path for x in path_sets[path[i]]]:
                    path_sets[path[i]].append(MyHeapObj(path[i:], new_mean, new_var)) 
            
            
            
            if path[i] not in R.nodes:
                R.add_node(path[i])
                
            #check if edge not already in graph?
            if (path[i-1], path[i]) not in R.edges:
                nd_edges.add((path[i-1], path[i]))
                R.add_edge(path[i - 1], path[i], mean = graph.edges[path[i - 1], path[i]]['mean'], var = graph.edges[path[i - 1], path[i]]['var'])
                        
#    print_graph_color(graph, nd_edges, 'r')      
    return (R, path_sets, nd_edges)


# =============================================================================
# Functions for calculating probability of paths giving us lower cost
# =============================================================================
def d(x, c_0, mean, var):
   # return (x - c_0 - mean) * 1.0/(math.sqrt(2 * var))
    return (x - c_0 - mean) * 1.0/(math.sqrt(2) * var)

def mult_term(x, path_set, edgecost1, edgecost2):
    mult = 1
    for i in range(len(path_set)):
        curr_path = path_set[i]
        mult *= 0.5 * erfc(d(x, edgecost1 - edgecost2, curr_path.path_mean, curr_path.path_var))
    return mult

def f(x, path_set, edgecost1, edgecost2):
    """
    x - the variable to integrate over
    node - the node which we are comparing
    path_sets - the set of paths from this node to the end
    edgecost - the concrete edge cost from the current node in the path search to this node
    """
    test = 1 - mult_term(x, path_set, edgecost1, edgecost2)
    return test
 
def g(x, path_set, edgecost):
    ret = 0
    for i in range(len(path_set)):
        b = copy.deepcopy(path_set)
        b.pop(i)
        factor1 = mult_term(x, b, edgecost, 0)
        factor2 = math.exp(-1.0 * (d(x, edgecost, 0, path_set[i].path_var))**2)/(math.sqrt(2 * math.pi)* path_set[i].path_var )

        # factor2 = math.exp(-1.0 * (d(x, edgecost, path_set[i].path_mean, path_set[i].path_var))**2)/(math.sqrt(2 * math.pi)* path_set[i].path_var )
        # factor2 = math.exp(-1.0 * (d(x, edgecost, path_set[i].path_mean, path_set[i].path_var))**2)/(math.sqrt(2 * math.pi* path_set[i].path_var ))
        ret += factor1 * factor2
        
    return ret

def test_compare(edgecost1, edgecost2, paths1, paths2):
    prob = scipy.integrate.quad(lambda x: f(x,  paths1, edgecost1) * g(x, paths2, edgecost2), -np.inf, np.inf)
    # print(prob)
    return None

def prob_comparator(x, paths, node1, node2, cost1, cost2):
    f_sub = 1.0
    for i in range(len(paths[node1])):
        mean_i = paths[node1][i].path_mean
        var_i = paths[node1][i].path_var
        f_sub *= 0.5 * erfc((x - mean_i - cost1) * 1.0/math.sqrt(2 * var_i))
              
    f = 1 - f_sub

    g = 0
    for j in range(len(paths[node2])):
        mean_j = paths[node2][j].path_mean
        var_j = paths[node2][j].path_var
        factor1 = math.exp(-1 * ((x - cost2 - mean_j)/math.sqrt(2 * var_j))**2)/(math.sqrt(2 * math.pi * var_j))
        factor2 = 1.0
        for i in range(len(paths[node2])):
            if i != j:
                mean_i = paths[node2][i].path_mean
                var_i = paths[node2][i].path_var
                factor2 *= 0.5 * erfc((x - mean_i - cost2)/math.sqrt(2 * var_i))

        g += factor1 * factor2

    return f * g


def comparePathSet(graph, current, nbrs, paths, end):
    """
    graph: a graph of non dominated paths
    nbrs: the nodes from which to consider paths to the end
    paths: a hashmap from nodes to a list of paths to the end along with the mean and variance
    
    Returns: the node to take to increase the likelihood of a low cost path
    """  
    #how to get the paths dict?
    # print(current)
    # print(nbrs)
    
    #only one neighbor or the end node is a neighbor -- don't do comparisons
    if len(nbrs) == 1:
        sample = np.random.normal(graph.edges[current, nbrs[0]]['mean'], graph.edges[current, nbrs[0]]['var'])
        if sample < 0:
            # print("ERROR prone")
            sample = 0

        return nbrs[0], sample
    elif end in nbrs:

        #*** this actually might give a bad path
        sample = np.random.normal(graph.edges[current, end]['mean'], graph.edges[current, end]['var'])        
        if sample < 0:
            # print("ERROR prone")
            sample = 0

        return end, sample
    else:
        costs = {}
        for nbr in nbrs:
            sample = np.random.normal(graph.edges[current, nbr]['mean'], graph.edges[current, nbr]['var'])
            if sample < 0:
                # print("ERROR prone")
                sample = 0                       
            costs[nbr] = sample
    
        #want to find the smallest by pairwise comparisons
        min_node = nbrs[0]
        idx = 1
        while idx < len(nbrs):
            #compare node1 and node2
            
            #probability that the cost of the current min node will be less than the cost of the ith nbr
            if costs[min_node] == 0:
                idx += 1
                continue
            elif costs[nbrs[idx]] == 0:
                min_node = nbrs[idx]
                idx +=1
                continue
            # prob = scipy.integrate.quad(lambda x: f(x,  paths[min_node], costs[min_node], costs[nbrs[idx]]) * g(x, paths[nbrs[idx]], costs[nbrs[idx]]), -np.inf, np.inf)
            prob = scipy.integrate.quad(lambda x: prob_comparator(x, paths, min_node, nbrs[idx], costs[min_node], costs[nbrs[idx]]), -np.inf, np.inf)


            #This should never happen
            # if prob[0]<0 or prob[0] >1:
            #     print("ERROR -- probability is " + str(prob))
            #     #prints for testing
            #     print("")
            #     print("comparing")
            #     print(str(min_node) + " " + str(nbrs[idx]))
            #     print(str(costs[min_node]) + " " + str(costs[nbrs[idx]]))
            #     for path in paths[min_node]:
            #         print(str(path.path) +  " " + str(path.path_mean ) + " " + str(path.path_var))
            #     for path in paths[nbrs[idx]]:
            #         print(str(path.path) +  " " + str(path.path_mean ) + " " + str(path.path_var))
            #     return None            

            if prob[0] < 0.5:
                min_node = nbrs[idx]

            idx += 1
        return min_node, costs[min_node]
       


def find_path(graph, start, end, threshold, horizon = None):
    """
        Input: networkx graph with edge prob dists and a start and end node
        
        Output: sequence of nodes to give least cost, least risk path from start to end
    """   
    global D_thresh
    D_thresh = threshold
    
    #Phase 1 - construct non dominated part of graph
    opened = []
    heapify(opened)
    
    closed = defaultdict(list)
    p_nd = set([])
    
    path_s = MyHeapObj([start], 0, 0)
    heappush(opened, path_s)
        
    while len(opened) > 0:
        path_cur = heappop(opened)
        # print(path_cur.path)
        last_node = path_cur.path[-1]
        closed[last_node].append(path_cur)
        
        if last_node == end:
            #Found a non dominated path to the end node
            # print("adding path to p_nd")
            p_nd.add(path_cur) 
        else:
            #Add all the acyclic paths to the open queue if they are not dominated by a seen
            #path to the node and 
            for nbr in graph.neighbors(last_node):
                if nbr not in path_cur.path:
                    #the path is acyclic
                    
                    #construct new path object
                    new_mean = path_cur.path_mean + get_mean(graph, nbr, last_node)
                    new_var = path_cur.path_var + get_var(graph, nbr, last_node)
                    new_path = copy.deepcopy(path_cur.path)
                    new_path.append(nbr)
                    path_n = MyHeapObj(new_path, new_mean, new_var)
                    
                    if nonDombySeen(path_n, nbr, closed):
                        heappush(opened, path_n)

        #if the lightest path to be considered is dominated by a path to the end then
        #end this sweep of the graph
        if len(opened)>0 and DombyPaths(nsmallest(1, opened)[0], p_nd):
            break

    #get the non dominated graph
    G_nd, path_sets, nd_edges = getNdPortion(graph, p_nd)
    
    # print('Phase 2 - construct opt path')
    #prints path_Sets
    # for node, paths in path_sets.items():
    #     for heapobj in paths:
            # print(heapobj.path, heapobj.path_mean, heapobj.path_var)                        
    
    final = []
    actual_cost = 0
    mean_cost = 0
    current = start
    rest = []
    if horizon != None:
        final.append(start)
        #get the rags path up to some horizon and then do a* on means
        for i in range(horizon):
            # print(i)
            # Get all the neighbors from the non dom paths that dont cause a cycle
            nbrs = [path.path[1] for path in path_sets[current] if path.path[1] not in final]
            
            next_node, seen_cost = comparePathSet(G_nd, current, list(set(nbrs)), path_sets, end)

            actual_cost += seen_cost
            mean_cost += G_nd.edges[current, next_node]['mean']
            current = next_node
            final.append(current)

            if current == end:
                break
        print(final)

        #get the a* path for the rest
        # if final[-1] != end and to_end:
        #     rest_path = nx.astar_path(G_nd, final[-1], end, weight = 'mean')
        #     # print(rest_path)
        #     # print_graph_horizon(graph, nd_edges, final, rest_path, 'r', 'b', 'g')
        #     final.extend(rest_path[1:])
            # print(final)
        rest = nx.astar_path(G_nd, final[-1], end, weight = 'mean')

        # else:
        #     print_graph_2colors(graph, nd_edges, final, 'r', 'b')


    else:
        while current != end:
            final.append(current)
            
            # Get all the neighbors from the non dom paths that dont cause a cycle
            nbrs = [path.path[1] for path in path_sets[current] if path.path[1] not in final]
            
            next_node, seen_cost = comparePathSet(G_nd, current, list(set(nbrs)), path_sets, end)

            actual_cost += seen_cost
            mean_cost += G_nd.edges[current, next_node]['mean']
            current = next_node


        final.append(end)

        print_graph_2colors(graph, nd_edges, final, 'r', 'b')
    

    
    #prints the final graph
    # print_graph_color(graph, nd_edges, 'r')

    return final, rest, actual_cost, mean_cost

# ===========================================================================
# ===========================================================================
# #THESE WRAPPERS ARE NO LONGER REQUIRED
# #because map.py now has a networkx graph object stored

# def find_path_wrapper(graph, start, end, threshold, horizon = None):
#     """
#         Converts an igraph representation to a networkx representation
#         and gives this as input to the find path function which executes the RAGS algorithm
#     """

#     #convert to networkx
#     nxG = nx.Graph(graph.get_edgelist())
#     for e in graph.es:
#         nxG.edges[e.source, e.target]["mean"] = e["mean"]
#         nxG.edges[e.source, e.target]["var"] = e["var"]

#     path, rest, actual_cost, mean_cost = find_path(nxG, start, end, threshold, horizon)
#     print(str(actual_cost) + "  " + str(mean_cost))
#     return path, rest


# def find_partial_path_wrapper(graph, start, end, threshold, horizon = None):
#     """
#         Converts an igraph representation to a networkx representation
#         and gives this as input to the find path function which executes the RAGS algorithm
#     """

#     #convert to networkx
#     nxG = nx.Graph(graph.get_edgelist())
#     for e in graph.es:
#         nxG.edges[e.source, e.target]["mean"] = e["mean"]
#         nxG.edges[e.source, e.target]["var"] = e["var"]

#     path, actual_cost, mean_cost = find_path(nxG, start, end, threshold, horizon, False)
#     print(str(actual_cost) + "  " + str(mean_cost))
#     return path

# ===========================================================================
# ===========================================================================
    

#===============================================
#TESTING

# for i in range(100):
#     print(i)
#     print(find_path(T, 0, 7, 0.85, 3))
    # print(find_path(S, 0, 5, 0.99, 2))
    # print(find_path_wrapper(G1, 0, 5, 0.65))
# print("rags.py")
# D_thresh = 0.5
# print(nonDom(4,2, 5,2))

# def ftest(x):
#     return (1 - 0.5 * erfc((x - 9.3106)/math.sqrt(2))) * (math.exp(-(((x - 21.999)/2)**2))/math.sqrt(2 * math.pi * 2))

# print(scipy.integrate.quad(lambda x: ftest(x),-np.inf, np.inf))


