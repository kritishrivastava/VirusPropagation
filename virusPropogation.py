# Implementation for Project 5 (Option 1): Virus Propagation on Static Networks

# Project Team Members:
# 1. Manjusha Trilochan Awasthi (mawasth)
# 2. Kriti Shrivastava (kshriva)
# 3. Rachit Thirani (rthiran)

import networkx as nx
import numpy as np
import os
import matplotlib.pyplot as plt
import random
from operator import itemgetter
import scipy

def readGraph(graphName):
    #Input: name of the file from which graph has to be read
    #Output: creates and returns a networkx graph
    #Read file and create a networkx graph
    graph = nx.Graph()
    filePath = os.getcwd()+ "/" +graphName
    f = open(filePath, 'r')
    #Skipping the first row containing number of nodes and edges
    next(f)
    for line in f:
        line = line.split()
        graph.add_edge(int(line[0]),int(line[1]))
    return graph

def calculateEffectiveStrength(b , d, lambda1):
    #Input: transmission probability b, healing probability d and highest eigenvalue for adjacency matrix lambda1
    #Output: effective strength s for the SIS Virus Propagation Model, min transmission probability and max healing probability for epidemic
    #Calculating Cvpm (constant that depends on the virus propagation model) using the formula Cvpm = b/d for SIS model
    Cvpm = b/d
    #Calculating strength using the formula s = lambda1 * CVPM
    strength = lambda1 * Cvpm
    print("Effective strength: ",np.real(strength))
    if strength > 1:
        print("Since Effective strength is greater than 1, the infection will result in an epidemic.")
    else:
        print("Since Effective strength is less than 1, the virus will die quickly.")
    #Find minimum transmission probability (b) that results in a networkwide epidemic i.e strength > 1
    min_beta = d / lambda1
    print("Minimum transmission probability for epidemic: ",np.real(min_beta))
    #Find maximum healing probability (d) that results in a networkwide epidemic i.e. strength > 1
    max_delta = b * lambda1
    print("Maximum healing probability for epidemic: ",np.real(max_delta))


def plotSvsB(b, d, lambda1, part):
    #Input: transmission probability b, healing probability d and highest eigenvalue for adjacency matrix lambda1
    #Output: Plots a graph for Relationship between Virus Strength and Transmission probability (beta)
    betas = np.arange(0, b+0.001, 0.001)
    strengths = [lambda1*(beta/d) for beta in betas]
    plt.figure(figsize=(12, 6))
    plt.plot(betas, strengths, "-o")
    plt.title("Relationship between Virus Strength and Transmission probability (beta)")
    plt.xlabel('Beta')
    plt.ylabel('Effective Strength')
    plotFile = "output/" + "betaVSstrength_"+ part +".png"
    #os.makedirs(os.path.dirname(plotFile), exist_ok=True)
    plt.savefig(plotFile, bbox_inches='tight')


def plotSvsD(b, d, lambda1, part):
    # Input: transmission probability b, healing probability d and highest eigenvalue for adjacency matrix lambda1
    # Output: Plots a graph for Relationship between Virus Strength and Healing probability (delta)
    deltas = np.arange(0.01, d+0.1, 0.1)
    strengths = [lambda1 * (b / delta) for delta in deltas]
    plt.figure(figsize=(12, 6))
    plt.plot(deltas, strengths, "-o")
    plt.title("Relationship between Virus Strength and Healing probability (delta)")
    plt.xlabel('Delta')
    plt.ylabel('Effective Strength')
    plotFile = "output/" + "deltaVSstrength_"+ part +".png"
    #os.makedirs(os.path.dirname(plotFile), exist_ok=True)
    plt.savefig(plotFile, bbox_inches='tight')

def sis_vpm(graph,transmission,healing):
    # Input: netwworkx graph, transmission probability and healing probability
    # Output: Fraction of infected nodes at time steps ranging from 0 to 100
    nodes = graph.nodes()
    #Fraction of nodes which are infected
    fraction= list()
    #Find initially infected nodes
    infected_nodes = set(np.random.choice(nodes, int(len(nodes)/10), replace = False))
    #Considering that 1/10th of the nodes are infected initially at time step 0
    fraction.append(len(infected_nodes) / len(nodes))
    #Find number of infected nodes at the remaining 99 time steps.
    for i in range(99):
        #Finding the nodes from the currently infected nodes which have healed
        healed_nodes=set(np.random.choice(list(infected_nodes), int(healing*len(infected_nodes))))
        #Removing the healed nodes from the list of infected nodes
        infected_nodes=infected_nodes.difference(healed_nodes)
        #Finding neighbors of the infected nodes which are at a risk of infection
        for node in infected_nodes:
            neighbors=graph.neighbors(node)
            #From the list of neighbors of the infected node, removing the nodes which are already infected
            healthy_nodes = set(neighbors).difference(infected_nodes)
            #If the infected node has some healthy neighbors, find ones which get infected
            if(len(healthy_nodes)) > 0:
                #Adding the newly infected nodes to the list of infected nodes
                infected_nodes=infected_nodes.union(set(np.random.choice(neighbors, int(transmission*len(neighbors)))))
        #Finding the fraction of infected nodes at every time step
        fraction.append(float(len(infected_nodes))/float(len(nodes)))
    return fraction

def plotInfectionvsTime(graph, b, d, part):
    # Input: networkx graph, transmission probability b, healing probability d
    # Output: Plots a graph between the average fraction of infected nodes vs time
    infectionFraction = list()
    #Perform 10 simulations of virus propagation
    for i in range(10):
        #Calculate the fraction of infected nodes for SIS virus propagation model
        f = sis_vpm(graph, b, d)
        infectionFraction.append(f)
    #Calucate average of fraction of infected nodes at every time step
    avgInfectionFraction = [float(sum(col)) / len(col) for col in zip(*infectionFraction)]
    #Plot graph for avg infection fraction vs time
    plt.figure(figsize=(12, 6))
    plt.plot(avgInfectionFraction)
    plt.xlabel("Time")
    plt.ylabel("Average fraction of Population Infected")
    plotFile = "output/" + "infectionVStime_" + part + ".png"
    plt.savefig(plotFile, bbox_inches='tight')


def plotSvsk(b, d, graph, policy):
    # Input: transmission probability, healing probability, networkx graph and the policy number
    # Output: Plots a graph between the virus strength and the number of vaccines
    strengths = []
    kStart = 0
    kEnd = nx.number_of_nodes(graph)
    step = 500
    rangeOfK = [k for k in range(kStart,kEnd,step)] 
    for k in rangeOfK:
        if policy == 'A':
            immunizedGraph = immunizeUsingPolicyA(graph.copy(), k) 
        elif policy == 'B':
            immunizedGraph = immunizeUsingPolicyB(graph.copy(), k) 
        elif policy == 'C':
            immunizedGraph = immunizeUsingPolicyC(graph.copy(), k)  
        elif policy == 'D':
            immunizedGraph = immunizeUsingPolicyD(graph.copy(), k)
        lambda1 = scipy.sparse.linalg.eigs( nx.to_numpy_matrix(immunizedGraph), k=1, which='LM')[0]
        Cvpm = b/d
        strength = np.real(lambda1 * Cvpm)
        if strength < 1:
	        print("Effective strength is less than 1, the virus will die quickly for k :", k)
        strengths.append(strength)
    plt.figure(figsize=(12, 6))
    plt.plot(rangeOfK, strengths, "-o")
    plt.title("Relationship between Virus Strength and number of vaccines needed to prevent epidemic")
    plt.xlabel('Number of Vaccines: k')
    plt.ylabel('Virus Strength: s')
    plt.axhline(y=1, linewidth=2, color="r")
    plt.text(kStart+50, 1.2 , 'Epidemic threshold i.e s=1')
    plotFile = "output/" + "kVSstrength_policy"+ policy +".png"
    #os.makedirs(os.path.dirname(plotFile), exist_ok=True)
    plt.savefig(plotFile, bbox_inches='tight')


def immunizeUsingPolicyA(graph, k):
    # Immunize and remove random k nodes from the network
    kRandomNodes = random.sample(range(0, nx.number_of_nodes(graph)), k)
    graph.remove_nodes_from(kRandomNodes)
    return graph

    
def immunizeUsingPolicyB(graph, k):
    # Immunize and remove k highest degree nodes from the network
    kHighestDegreeNodes = [node for node,degree in sorted(graph.degree_iter(),key=itemgetter(1),reverse=True)][:k]
    graph.remove_nodes_from(kHighestDegreeNodes)
    return graph
  

def immunizeUsingPolicyC(graph, k):
    # Immunize and remove the highest degree node iteratively(k times) from the network
    for i in range(k):
        highestDegreeNode = max(graph.degree_iter(), key=lambda node_degree: node_degree[1])[0]
        graph.remove_node(highestDegreeNode)
    return graph


def immunizeUsingPolicyD(graph, k):
    # Immunize and remove the nodes corresponding to the k largest values in the eigen vector of the largest eigen value
    largestEigenValue, largestEigenVector = scipy.sparse.linalg.eigs( nx.to_numpy_matrix(graph), k=1, which='LM', return_eigenvectors=True)
    absValuesWithIndex = []
    for index, value in enumerate(largestEigenVector):
        absValuesWithIndex.append((index, abs(value)))
        kCorrespondingNodes = [absValueWithIndex[0] for absValueWithIndex in sorted(absValuesWithIndex, key = lambda absValueWithIndex : absValueWithIndex[1], reverse = True)][:k]
    graph.remove_nodes_from(kCorrespondingNodes)
    return graph


if __name__ == "__main__":
    #Read file and create a graph
    filename = "static.network"
    graph = readGraph(filename)

    #Find largest eigenvalue of the adjacency matrix of the graph
    #adjacencySpectrum = nx.adjacency_spectrum(graph)
    lambda1 = scipy.sparse.linalg.eigs( nx.to_numpy_matrix(graph), k=1, which='LM')[0]

    #Calculate the effective strength of the virus on the static contact network
    ### Part 1 for b1 and d1
    # Setting the value for Transmission probability(b) and healing probability(d)
    b = 0.20
    d = 0.70
    print("-------------  Results for beta1(0.20) and delta1(0.70)  -------------")
    calculateEffectiveStrength(b, d, lambda1)
    #Plot graph for strength vs b (d fixed)
    plotSvsB(b, d, lambda1, "part1")
    #Plot graph for strength vs d (b fixed)
    plotSvsD(b, d, lambda1, "part1")

    ### Part 2 for b2 and d2
    # Setting the value for Transmission probability(b) and healing probability(d)
    b = 0.01
    d = 0.60
    print("\n\n-------------  Results for beta1(0.01) and delta1(0.60)  -------------")
    calculateEffectiveStrength(b, d, lambda1)
    # Plot graph for strength vs b (d fixed)
    plotSvsB(b, d, lambda1, "part2")
    # Plot graph for strength vs d (b fixed)
    plotSvsD(b, d, lambda1, "part2")

    ### 2) Virus propagation simulation across network
    ### Part 1 for b1 and d1
    b = 0.20
    d = 0.70
    part = "part1"
    plotInfectionvsTime(graph, b, d, part)
    ### Part 2 for b2 and d2
    b = 0.01
    d = 0.60
    part = "part2"
    plotInfectionvsTime(graph, b, d, part)


    ### 3) Immunization Policy for preventing the spread of epidemic
    # 3.c) Psuedo-code for each immunization policy
    k = 200
    immunizedGraphA = immunizeUsingPolicyA(graph.copy(), k)
    immunizedGraphB = immunizeUsingPolicyB(graph.copy(), k)
    immunizedGraphC = immunizeUsingPolicyC(graph.copy(), k)
    immunizedGraphD = immunizeUsingPolicyD(graph.copy(), k)
 
    # 3.d) Calculating effective strength of virus on immunized network for each immunization policy
    b = 0.20
    d = 0.70
    print("\n\n-------------  Calculating effective strength(s) of virus on network immunized using policy A :  -------------")  
    lambda1 = scipy.sparse.linalg.eigs( nx.to_numpy_matrix(immunizedGraphA), k=1, which='LM')[0]
    calculateEffectiveStrength(b, d, lambda1)

    print("\n\n-------------  Calculating effective strength(s) of virus on network immunized using policy B :  -------------")
    lambda1 = scipy.sparse.linalg.eigs( nx.to_numpy_matrix(immunizedGraphB), k=1, which='LM')[0]
    calculateEffectiveStrength(b, d, lambda1)

    print("\n\n-------------  Calculating effective strength(s) of virus on network immunized using policy C :  -------------")
    lambda1 = scipy.sparse.linalg.eigs( nx.to_numpy_matrix(immunizedGraphC), k=1, which='LM')[0]
    calculateEffectiveStrength(b, d, lambda1)

    print("\n\n-------------  Calculating effective strength(s) of virus on network immunized using policy D :  -------------")
    lambda1 = scipy.sparse.linalg.eigs( nx.to_numpy_matrix(immunizedGraphD), k=1, which='LM')[0]
    calculateEffectiveStrength(b, d, lambda1)

    # 3.e) Estimating minimum number of vaccines needed to prevent epidemic for each immunization policy
    print("\n\n-------------  Finding minimum number of vaccines needed to prevent epidemic on network immunized using policy A :  -------------")
    plotSvsk(b, d, graph.copy(), "A")
    print("\n\n-------------  Finding minimum number of vaccines needed to prevent epidemic on network immunized using policy B :  -------------")
    plotSvsk(b, d, graph.copy(), "B")
    print("\n\n-------------  Finding minimum number of vaccines needed to prevent epidemic on network immunized using policy C :  -------------")
    plotSvsk(b, d, graph.copy(), "C")
    print("\n\n-------------  Finding minimum number of vaccines needed to prevent epidemic on network immunized using policy D :  -------------")
   # plotSvsk(b, d, graph.copy(), "D")

    # 3.f) Plot the average fraction of infected nodes at each time step for immunized networks
    b = 0.20
    d = 0.70
    part = "policyA"
    plotInfectionvsTime(immunizedGraphA, b, d, part)
    part = "policyB"
    plotInfectionvsTime(immunizedGraphB, b, d, part)
    part = "policyC"
    plotInfectionvsTime(immunizedGraphC, b, d, part)
    part = "policyD"
    plotInfectionvsTime(immunizedGraphD, b, d, part)


 
