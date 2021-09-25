
#recursive bond breaking example (methanol)
from recursive_bond_breaking import *

mol_dict = generate_mol_dict('CO','smi') #create molecule dictionary from SMILES notation
G_mol = mol_graph(mol_dict) #generate molecule graph
all_groups = subgroup([G_mol]) #generate all possible intermediates by recursive bond cleavage 

#save the intermediates graph to .gpickle file
for i, group in enumerate(all_groups): 
    nx.write_gpickle(group, './%i.gpickle' % i)
all_group_dict = all_groups_dict(all_groups)


#traj to SMILES example
from TrajToGraph import *
from GraphToSmiles import *

remove_lower_metal('example.traj','example_r.traj',16.5) #remove all lower metal atoms 
s = locate_adsorbate('example_r.traj','Rh') #pbc condition
graph_traj('example_r.traj','graph.traj',s,'Rh',14.5,2.9,2.75,2.3) #remove all metal atoms not in direct chemical interactions with the adsorbates
adsorbate_graph('graph.traj','graph.gpickle', 'Rh', 1.7, 1.9, 1.3, 1.3, 2.9, 2.75, 2.3) #generate graph representation from trajectory

node_list, adjacency_matrix = gpickle_to_matrix('graph.gpickle') #generate node and edge information from graph
smi = MolFromGraphs(node_list, adjacency_matrix) #generate SMILES notation
print(smi)
