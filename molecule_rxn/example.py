
#recursive bond breaking example (methanol)
from recursive_bond_breaking import *

mol_dict = generate_mol_dict('CO','smi')
G_mol=mol_graph(mol_dict)
all_groups=subgroup([G_mol])
for i, group in enumerate(all_groups):
    nx.write_gpickle(group, './%i.gpickle' % i)
all_group_dict=all_groups_dict(all_groups)


#traj to SMILES example
from TrajToGraph import *
from GraphToSmiles import *

remove_lower_metal('example.traj','example_r.traj',16.5)
s=locate_adsorbate('example_r.traj','Rh')
graph_traj('example_r.traj','graph.traj',s,'Rh',14.5,2.9,2.75,2.3)
adsorbate_graph('graph.traj','graph.gpickle', 'Rh', 1.7, 1.9, 1.3, 1.3, 2.9, 2.75, 2.3)

node_list, adjacency_matrix = gpickle_to_matrix('graph.gpickle')
smi = MolFromGraphs(node_list, adjacency_matrix)
print(smi)
