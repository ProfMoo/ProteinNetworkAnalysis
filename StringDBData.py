import networkx as nx
import matplotlib.pyplot as plt

'''
Enter the following file names ensuring that they are saved within the 
same folder as this code file. Three files are the following:

- A text file with the symbols of the core proteins for the disease
- The tsv file from STRING exports: "simple tabular text output"
- The network coordinates text file from STRING

If desired, this part may be altered such that the user enters the file name in
in the command window by using the input command.
'''

## Creates a list of all the core nodes to be used later for comparison
def create_list_of_core_nodes(core_file):
    core_lines = core_file.readlines()
    core_nodes = []
    for c in core_lines:
        c = c.strip('\n')
        core_node = c
        core_nodes.append(core_node)

    return core_nodes

'''
Identifies important info about nodes, ENSP ids, and connected nodes and 
sorts them into useful lists.
- nodes_in_network: list of all the nodes in the extended network
- edge_list: list of all the node pairs connected by an edge in the network
- id_list: list of all the nodes in the network with their respective ENSP ID's
'''
def get_interaction_data(interaction_lines):
    interaction_lines = interaction_lines.readlines()
    nodes_in_network = []
    edge_list = []
    id_list = []
    for i in interaction_lines:
        interaction = i.split('\t')
        node1 = interaction[0]
        node2 = interaction[1]
        ensp1 = interaction[4].strip("9606.")
        ensp2 = interaction[5].strip("9606.")
        pair = (node1, node2)
        id1 = (node1, ensp1)
        id2 = (node2, ensp2)
        if pair not in edge_list:
            edge_list.append(pair)
        if node1 not in nodes_in_network:
            nodes_in_network.append(node1)
            id_list.append(id1)
        if node2 not in nodes_in_network:
            nodes_in_network.append(node2) 
            id_list.append(id2) 

    return edge_list, id_list, nodes_in_network

## Counts the number of core nodes that appear in the connected extended network   
def get_num_core_nodes(nodes_in_network):
    numcore = 0    
    for nodes in nodes_in_network:
        if nodes in core_nodes:
            numcore +=1

    return numcore

'''
Allows user to obtain a description of all nodes in the network as well 
as their coordinates and colors. The description includes the protein's symbol,
name, and function. This information can be printed. The annotation is all the 
information separated into a list. The prot_name is the full name of the
protein.
'''
def pull_description_of_protein(coordinate_file):
    coor_lines = coordinate_file.readlines()
    detail_list = []
    for j in coor_lines:
        description = j.split('\t')
        node = description[0]
        annotation = description[4]
        annotation = annotation.split(';')
        prot_name = annotation[0]
        details = (node, prot_name)
        detail_list.append(details)

    return detail_list

def create_networkx_graph(edge_list, nodes_in_network):
    ## Recreates the STRING network using NetworkX
    G = nx.Graph()
    G.add_edges_from(edge_list)
    G.edges()

    ## Caluclates the degree for each node and the average degree in the network
    deg = G.degree(nodes_in_network)
    deglist = []
    for k,v in deg:
        deglist.append(v)
    avg_degree = sum(deglist) / len(deglist)

    ## Calculates the clustering coefficient for each node in the network
    cluster = (nx.clustering(G))
    avg_cluster = (nx.average_clustering(G))
    for dict_value in cluster:
        for k, v in cluster.items():
            cluster[k] = round(v, 4) 

    ## Calculated the betweenness coefficient for each node in the network
    bet = nx.betweenness_centrality(G)
    betlist = []
    for (k,v) in bet.items():
        betlist.append(v)
        #print(k,v)
    avg_bet = sum(betlist) / len(betlist)

    ## Prints a list of all cliques of a certain size or larger within the network
    cliques = list(nx.find_cliques(G))
    i = 1
    for j in cliques:
        if len(j) > 5:
            #print("Clique", str(i) + ':', j)
            i += 1

    ## Creates a list of data for each node with includes clustering coefficient and deg
    node_data = []
    for n in G:
        nlist = [n, nx.clustering(G,n), nx.degree(G,n)]
        node_data.append(nlist)

    print("The average degree for nodes in this network is:", avg_degree)
    print("The average clustering coefficient in this network is:", avg_cluster)
    print("The average betweenness centrality in this network is:", avg_bet, "\n")

    return G

def draw_graph(G):
    ##Created network image
    nx.draw(G, with_labels = True, font_size = 8)

    ## Displays network image
    plt.show()

if __name__ == "__main__":
    core_file = open("PCOSCore.txt", "r")
    interaction_lines = open("PCOS920.tsv", "r")
    coordinate_file = open("PCOS920nc.txt", "r")

    core_nodes = create_list_of_core_nodes(core_file)

    edge_list, id_list, nodes_in_network = get_interaction_data(interaction_lines)
    
    numcore = get_num_core_nodes(nodes_in_network)

    detail_list = pull_description_of_protein(coordinate_file)

    G = create_networkx_graph(edge_list, nodes_in_network)

    draw_graph(G)

    ## Prints a summary of all network information
    print("The number of total nodes in the extended network is:", len(nodes_in_network))
    print("The number of total edges in the extended network is:", len(edge_list))
    print("The number of core nodes in this network is:", numcore)