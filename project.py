import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt

# Load the PPI network data
df = pd.read_csv('PathLinker_2018_human-ppi-weighted-cap0_75.txt', sep='\t', header=None, dtype={2: object})
df.columns = ['tail', 'head', 'edge_weight', 'edge_type']

# Create a directed graph
G = nx.from_pandas_edgelist(df, source='tail', target='head', edge_attr='edge_weight', create_using=nx.DiGraph())

# Function to find the shortest path between two proteins
def shortest_path(protein1, protein2):
    if protein1 not in G or protein2 not in G:
        return "Proteins not found in the interactome"
    try:
        path = nx.shortest_path(G, source=protein1, target=protein2, weight='edge_weight')
        return path
    except nx.NetworkXNoPath:
        return "No path found between the proteins"

# Function to list all directly connected proteins
def connected_proteins(protein):
    if protein not in G:
        return "Protein not found in the interactome"
    neighbors = list(G.neighbors(protein))
    return neighbors

# Function to draw a histogram for the degrees of a set of proteins
def degree_histogram(proteins):
    degrees = [G.degree(p, None) for p in proteins]
    plt.hist(degrees, bins=max(degrees)+1, edgecolor='black', alpha=0.7)
    plt.xlabel('Degree')
    plt.ylabel('Frequency')
    plt.title('Degree Histogram')
    plt.show()

# Function to rank proteins based on their degrees
def rank_proteins(proteins):
    degrees = [(p, G.degree(p)) for p in proteins]
    degrees.sort(key=lambda x: x[1], reverse=True)
    return degrees

# Function to map protein UniProt ID to its gene name
def map_id_to_gene(protein_id):
    # This function needs to be implemented based on how you can map the protein ID to gene name
    pass

# Convert the weighted graph to an unweighted graph and save it as an adjacency matrix
G_unweighted = nx.Graph(G)
adj_matrix = nx.adjacency_matrix(G_unweighted)
print("Adjacency Matrix:\n", adj_matrix.toarray())

# Example usage
protein1_id = 'PROT1'
protein2_id = 'PROT2'

print("Shortest Path:", shortest_path(protein1_id, protein2_id))
print("Connected Proteins:", connected_proteins(protein1_id))

# Test degree histogram
proteins_to_test = [protein1_id, protein2_id]
degree_histogram(proteins_to_test)

# Test ranking proteins
ranked_proteins = rank_proteins(proteins_to_test)
print("Ranked Proteins:", ranked_proteins)
