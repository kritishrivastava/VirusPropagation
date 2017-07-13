# Project 5: Virus Propagation on static networks
Project members:
1. Manjusha Trilochan Awasthi (mawasth)
2. Kriti Shrivastava (kshriva)
3. Rachit Thirani (rthiran)

The goal of this project is to analyze the propagation of a virus in a static network and prevent a networkwide epidemic. Tasks performed: 
1. Analysis and understanding a virus propagation model.
2. Calculation of the effective strength of a virus.
3. Simulation of the propagation of a virus in a network.
4. Implementation of immunization policies to prevent a virus from spreading across a network.

The implementation was based on the findings from the paper "B. Aditya Prakash, Deepayan Chakrabarti, Michalis Faloutsos, Nicholas Valler, and Christos Falou tsos. Got the Flu (or Mumps)? Check the Eigenvalue! arXiv:1004.0060 [physics.socph], 2010".


Python version used: Python 3.6.0

OS used: Windows 10 64bit (RAM:16GB)

Python libraries needed:

1. Networkx: install using the command "pip install networkx" or follow the instructions mentioned here https://networkx.github.io/documentation/development/install.html
2. Numpy: install using the command "pip install numpy" or follow the intructions mentioned here https://docs.scipy.org/doc/numpy/user/install.html
3. Matplotlib: install using the command "pip install matplotlib" or follow the intructions mentioned here 
https://matplotlib.org/users/installing.html



Instructions to run the program:
1. Go the the directory containing the python script.
2. Give the command "python virusPropagation.py" 


Dataset used: static.network


Output of the program:
1. Numeric value of the effective strength (s) of the virus on the static contact network provided (static.network)
2. Answer to the question: "Will the infection spread across the network (i.e., result on an epidemic), or will it die quickly?"
3. Graph plot 1: Plot between the transmission probability(beta) and effective strength of the virus. The plot "staticGraph_betaVSstrength.png" is generated in the directory containing the python script under /output/. 
4. Graph plot 2: Plot between the Healing probability(delta) and effective strength of the virus. The plot "staticGraph_deltaVSstrength.png" is generated in the directory containing the python script under /output/. 
5. The minimum transmission probability (β) that results in a networkwide epidemic.
6. The maximum healing probability (δ) that results in a networkwide epidemic.
