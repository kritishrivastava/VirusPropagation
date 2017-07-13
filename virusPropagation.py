# Implementation for Project 5 (Option 1): Virus Propagation on Static Networks

# Project Team Members:
# 1. Manjusha Trilochan Awasthi (mawasth)
# 2. Kriti Shrivastava (kshriva)
# 3. Rachit Thirani (rthiran)

import networkx as nx
import numpy as np
import os
import matplotlib.pyplot as plt

def readGraph(graphName):
    #Input: name of the file from which graph has to be read
    #Output: creates and returns a networkx graph
    #Read file and create a networkx graph
    graph = nx.Graph()
    filePath = os.getcwd()+ "\\" +graphName
    f = open(filePath, 'r')
    #Skipping the first row containing number of nodes and edges
    next(f)
    for line in f:
        line = line.split()
        graph.add_edge(int(line[0]),int(line[1]))
    return graph

def calculateEffectiveStrength(b , d, lambda1):
    #Input: transmission probability β, healing probability δ and highest eigenvalue for adjacency matrix lambda1
    #Output: effective strength s for the SIS Virus Propagation Model, min transmission probability and max healing probability for epidemic
    #Calculating Cvpm (constant that depends on the virus propagation model) using the formula Cvpm = b/d for SIS model
    Cvpm = b/d
    #Calculating strength using the formula s = λ1 · CVPM
    strength = lambda1 * Cvpm
    print("Effective strength: ",np.real(strength))
    if strength > 1:
        print("Since Effective strength is greater than 1, the infection will result in an epidemic.")
    else:
        print("Since Effective strength is less than 1, the virus will die quickly.")
    #Find minimum transmission probability (β) that results in a networkwide epidemic i.e strength > 1
    min_beta = d / lambda1
    print("Minimum transmission probability for epidemic: ",np.real(min_beta))
    #Find maximum healing probability (δ) that results in a networkwide epidemic i.e. strength > 1
    max_delta = b * lambda1
    print("Maximum healing probability for epidemic: ",np.real(max_delta))


def plotSvsB(b, d, lambda1, part):
    #Input: transmission probability β, healing probability δ and highest eigenvalue for adjacency matrix lambda1
    #Output: Plots a graph for Relationship between Virus Strength and Transmission probability (beta)
    betas = np.arange(0, b+0.001, 0.001)
    strengths = [lambda1*(beta/d) for beta in betas]
    plt.figure(figsize=(12, 6))
    plt.plot(betas, strengths, "-o")
    plt.title("Relationship between Virus Strength and Transmission probability (beta)")
    plt.xlabel('Beta')
    plt.ylabel('Effective Strength')
    plotFile = "output/" + "betaVSstrength_"+ part +".png"
    os.makedirs(os.path.dirname(plotFile), exist_ok=True)
    plt.savefig(plotFile, bbox_inches='tight')


def plotSvsD(b, d, lambda1, part):
    # Input: transmission probability β, healing probability δ and highest eigenvalue for adjacency matrix lambda1
    # Output: Plots a graph for Relationship between Virus Strength and Healing probability (delta)
    deltas = np.arange(0.01, d+0.1, 0.1)
    strengths = [lambda1 * (b / delta) for delta in deltas]
    plt.figure(figsize=(12, 6))
    plt.plot(deltas, strengths, "-o")
    plt.title("Relationship between Virus Strength and Healing probability (delta)")
    plt.xlabel('Delta')
    plt.ylabel('Effective Strength')
    plotFile = "output/" + "deltaVSstrength_"+ part +".png"
    os.makedirs(os.path.dirname(plotFile), exist_ok=True)
    plt.savefig(plotFile, bbox_inches='tight')


if __name__ == "__main__":
    #Read file and create a graph
    filename = "static.network"
    graph = readGraph(filename)

    #Find largest eigenvalue of the adjacency matrix of the graph
    adjacencySpectrum = nx.adjacency_spectrum(graph)
    lambda1 = max(adjacencySpectrum)

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