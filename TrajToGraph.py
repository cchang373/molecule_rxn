from ase.io import read
import networkx as nx

def remove_lower_metal(traj_file_path, out_file_path, cut_off):
    """
    remove the lower layer metal
    """
    traj_file = read(traj_file_path)
    indices = [a.index for a in traj_file if a.z <= cut_off]
    del traj_file[indices]
    traj_file.write(out_file_path)
    return

#remove_lower_metal('../structure.traj','structure_e.traj',32.5)

def locate_adsorbate(traj_file_path, metal):
    traj_file = read(traj_file_path)
    metal_indices = [a.index for a in traj_file if a.symbol == metal]
    metal_positions = [traj_file.get_positions()[i] for i in metal_indices]
    x_min = min([position[0] for position in metal_positions])
    y_min = min([position[1] for position in metal_positions])
    x_max = max([position[0] for position in metal_positions])
    y_max = max([position[1] for position in metal_positions])

    non_metal_indices = [a.index for a in traj_file if a.symbol != metal]
    for i in non_metal_indices:
        x,y,z = traj_file.get_positions()[i]
        if x < x_min and y < y_min:
            return 'front_left'
        elif x > x_max and y < y_min:
            return 'front_right'
        elif x < x_min and y > y_max:
            return 'back_left'
        elif x > x_max and y > y_max:
            return 'back_right'

    positions = [traj_file.get_positions()[i] for i in non_metal_indices]
    if True in [position[0] < x_min and position[1] > y_min and position[1] < y_max for position in positions]:
        return 'left'
    elif True in [position[0] > x_max and position[1] > y_min and position[1] < y_max for position in positions]:
        return 'right'
    elif True in [position[0] < x_max and position[0] > x_min and position[1] > y_max for position in positions]:
        return 'back'
    elif True in [position[0] < x_max and position[0] > x_min and position[1] < y_min for position in positions]:
        return 'front'
    return 'inside'

def graph_traj(traj_file_path, out_file_path, adsorbate_location, metal, z, Cmetal_cutoff, Ometal_cutoff, Hmetal_cutoff):
    traj_file = read(traj_file_path+'@-1')
    traj_file.set_constraint()
    non_metal_index = [a.index for a in traj_file if a.symbol != metal]
    length = len(traj_file)
    
    traj_file_repeat = traj_file.repeat((2,2,1))

    if adsorbate_location == 'front_left':
        non_metal_index_repeat = [3*length+i for i in non_metal_index]
    elif adsorbate_location == 'back_left':
        non_metal_index_repeat = [2*length+i for i in non_metal_index]
    elif adsorbate_location == 'front_right':
        non_metal_index_repeat = [length+i for i in non_metal_index]
    elif adsorbate_location == 'back_right':
        non_metal_index_repeat = non_metal_index

    elif adsorbate_location == 'left':
        non_metal_index_repeat = [2*length+i for i in non_metal_index]
    elif adsorbate_location == 'right' or adsorbate_location == 'back' or adsorbate_location == 'inside':
        non_metal_index_repeat = non_metal_index
    elif adsorbate_location == 'front':
        non_metal_index_repeat = [length+i for i in non_metal_index] 
    
    non_metal_remove = [a.index for a in traj_file_repeat if (a.symbol != metal and a.index not in non_metal_index_repeat)]
    del traj_file_repeat[non_metal_remove]
    metal_lower = [a.index for a in traj_file_repeat if a.z < z]
    del traj_file_repeat[metal_lower]
    
    metal_index = [a.index for a in traj_file_repeat if a.symbol == metal]
    non_metal_l_index = [a.index for a in traj_file_repeat if a.symbol != metal]
    metal_keep = []
    for i in non_metal_l_index:
        elm = traj_file_repeat[i].symbol
        #if elm == 'H':
         #   continue
        for j in metal_index:
            distance = traj_file_repeat.get_distance(i,j)
            if (elm == 'C' and distance < Cmetal_cutoff) or (elm == 'O' and distance < Ometal_cutoff) or (elm == 'H' and distance < Hmetal_cutoff):
                metal_keep.append(j)
    metal_remove = [i for i in metal_index if i not in metal_keep]
   
    del traj_file_repeat[metal_remove]
    traj_file_repeat.write(out_file_path)
    return
#s=locate_adsorbate('test/H_test.traj','Rh')
#print(s)
#graph_traj('test/H_test.traj','test/H_test_r.traj',s,'Rh',16.7,2.9,2.75,2.3)

def adsorbate_graph(traj_file_path,out_file_path, metal, CC_cutoff, CO_cutoff, CH_cutoff, OH_cutoff, Cmetal_cutoff, Ometal_cutoff, Hmetal_cutoff):
    """
    generate graph for adsorbate and connected metal atoms
    C-C, C-O, C-H, O-H C-metal, O-metal bond cutoff distance
    """
    traj_file = read(traj_file_path)
    
    non_metal_index = [a.index for a in traj_file if a.symbol != metal]
    
    G = nx.Graph()

    for i in range(len(non_metal_index)):
        for j in range(i+1,len(non_metal_index)):
            elm_1 = traj_file[non_metal_index[i]].symbol
            position_1 = traj_file.get_positions()[non_metal_index[i]]
            elm_2 = traj_file[non_metal_index[j]].symbol
            position_2 = traj_file.get_positions()[non_metal_index[j]]
            distance = traj_file.get_distance(non_metal_index[i], non_metal_index[j])
            
            if ((elm_1 == 'C' and elm_2 == 'H') or (elm_1 == 'H' and elm_2 == 'C')) and (distance <= CH_cutoff):
                node_1 = elm_1 + str(position_1)
                node_2 = elm_2 + str(position_2)
                G.add_edge(node_1, node_2)
            elif ((elm_1 == 'C' and elm_2 == 'O') or (elm_1 == 'O' and elm_2 == 'C')) and (distance <= CO_cutoff):
                node_1 = elm_1 + str(position_1)
                node_2 = elm_2 + str(position_2)
                G.add_edge(node_1, node_2)
            elif ((elm_1 == 'O' and elm_2 == 'H') or (elm_1 == 'H' and elm_2 == 'O')) and (distance <= OH_cutoff):
                node_1 = elm_1 + str(position_1)
                node_2 = elm_2 + str(position_2)
                G.add_edge(node_1, node_2)
            elif elm_1 == 'C' and elm_2 == 'C' and distance <= CC_cutoff:
                node_1 = elm_1 + str(position_1)
                node_2 = elm_2 + str(position_2)
                G.add_edge(node_1, node_2)
   
    metal_index = [a.index for a in traj_file if a.symbol == metal]
    """
    #generate graph_1.gpickle with all possible bonds
    cutoff = {'C':Cmetal_cutoff, 'O':Ometal_cutoff, 'H':Hmetal_cutoff}
    for i in range(len(non_metal_index)):
        elm = traj_file[non_metal_index[i]].symbol
        position = traj_file.get_positions()[non_metal_index[i]]
        node = elm + str(position)
        #if elm == 'H':
         #   continue
        for j in range(len(metal_index)):
            distance = traj_file.get_distance(non_metal_index[i], metal_index[j])
            if distance <= cutoff[elm]:
                node_metal = metal + str(traj_file.get_positions()[metal_index[j]])
                G.add_edge(node, node_metal)
    """
           
    #generate graph.gpickle with permitted chemical bonds
    saturate = {'C':4, 'O':2, 'H':1}
    cutoff = {'C':Cmetal_cutoff, 'O':Ometal_cutoff, 'H':Hmetal_cutoff}
    for i in range(len(non_metal_index)):
        elm = traj_file[non_metal_index[i]].symbol
        position = traj_file.get_positions()[non_metal_index[i]]
        node = elm + str(position)
        if len(G.edges(node)) == saturate[elm]:
            continue
        num_bond = len(G.edges(node))
        unsaturate = saturate[elm] - num_bond
        
        for j in range(len(metal_index)):
            distance = traj_file.get_distance(non_metal_index[i], metal_index[j])
            if unsaturate == 0:
                break
            if distance <= cutoff[elm]:
                unsaturate -= 1
                node_metal = metal + str(traj_file.get_positions()[metal_index[j]]) 
                G.add_edge(node, node_metal)
    
    
    nx.write_gpickle(G,out_file_path)
    return

