import csv

def graph_handler(G):
    G_reader = csv.reader(open(G, "r"), delimiter= "\t")

    #Atlas Generator

    nodes_map = dict()
    id_node = 0

    for line in G_reader:

        if line[0] not in nodes_map.keys():
            nodes_map[line[0]] = id_node
            id_node += 1
        

        if line[1] not in nodes_map.keys():
            nodes_map[line[1]] = id_node
            id_node += 1
    

    G_reader = csv.reader(open(G, "r"), delimiter= "\t")


    dict_graph = {node: [0]*len(nodes_map) for node in range(id_node)}
    
    for line in G_reader:

        dict_graph[nodes_map[line[0]]][nodes_map[line[1]]] = float(line[2])
        dict_graph[nodes_map[line[1]]][nodes_map[line[0]]] = float(line[2])


    inv_map = {v: k for k, v in nodes_map.items()}

    return dict_graph, inv_map









def C_peeling(G, C):
    
    
    nodes = set(G.keys())
    V = len(nodes)

    #Degree Dictionary
    d = {k: [sum(x for x in v if x > 0), abs(sum(x for x in v if x < 0))] for k,v in G.items()}
    
    #Initial Solution
    densest_sub = list(d.keys())
    densest_avgdeg = sum([deg[0]-deg[1] for deg in d.values()])/len(d.keys())
    
    while len(d.keys()) > 2:
        
        #Compute Node to Remove
        node_to_remove = min(d.items(), key = lambda x : C*x[1][0] - x[1][1])[0]            
        
        
        #Update Degree of every Node
        for node in d:
            
            deg_to_remove = G[node][node_to_remove]
            
            if deg_to_remove > 0:
                d[node][0] -= deg_to_remove
            else:
                d[node][1] += deg_to_remove
        
        #Remove the Node
        del d[node_to_remove]
        
        #Check if solution has improved
        curr_avgdeg = sum([deg[0]-deg[1] for deg in d.values()])/len(d.keys())
        
        if curr_avgdeg > densest_avgdeg:
            densest_avgdeg = curr_avgdeg
            densest_sub = set(d.keys())
    
    return densest_sub, densest_avgdeg
























































