import numpy as np
import networkx as nx
from ase import Atoms
from ase.neighborlist import NeighborList,natural_cutoffs
from MD_python.chemical_environment import process_atoms, process_atoms2, find_H2O, hbond
from itertools import cycle
from pymatgen.io.vasp import Xdatcar
from vasppy.rdf import RadialDistributionFunction
import homcloud.interface as hc

class Anaylsis():
    def __init__(self,atoms,args) :
        self.atoms = atoms  
        self.graph = self.set_graph(args)  
        self.path_len = None
        self.glob_eff = None
        self.degree_hist = None
        self.pd = self.set_pd(args)

    def set_graph(self, args):
        self.atoms.set_pbc(args.pbc)
        positions = self.atoms.get_positions()
        nl = NeighborList(cutoffs=natural_cutoffs(self.atoms, mult=0.8),
                            self_interaction=False, bothways=True)
        nl.update(self.atoms)
        adjacency_matrix = nl.get_connectivity_matrix(sparse=False)
        graph = nx.from_numpy_matrix(adjacency_matrix)
        O_H_bonded,H_O_bonded= hbond(self.atoms) 
        for i in range(len(O_H_bonded)):
            graph.add_edge(O_H_bonded[i],H_O_bonded[i])
            merged_array = np.append(merged_array, positions[H_O_bonded[i]][2])
        for node, atom in zip(graph.nodes(),self.atoms):
            graph.nodes[node]['symbol'] = atom.symbol
            graph.nodes[node]['position'] = atom.position
        return graph
    
    def set_path_and_globeff(self):
        n = len(self.graph)
        denom = n * (n - 1)
        path_len = np.zeros(n,)
        if denom != 0:
            lengths = nx.all_pairs_shortest_path_length(self.graph)
            g_eff = 0
            for source, targets in lengths:
                for target, distance in targets.items():
                    if distance > 0:
                        g_eff += 1 / distance
                    path_len[int(distance)]+=1
            g_eff /= denom
        else:
            g_eff = 0
        path_len = path_len/2
        self.path_len = path_len
        self.glob_eff = g_eff
        return path_len, g_eff
    
    def set_degree_hist(self):
        self.degree_hist = nx.degree_hist(self.graph)
        return self.degree_hist

    def parse_key_value_pairs(items):
        result = {}
        items = items.split()
        for item in items:
            key, value = item.split("=")
            result[key] = value
        return result
        
    def seek_point(pos, point_list):
        result = []
        for index in range(len(pos)):
            points = pos[index][:-1]
            for tagret_point in point_list:
                if  np.linalg.norm(points - tagret_point)<1e-5:
                    result.append(index)
        return result
        
    def set_pd(self, args):
        pos = self.atoms.get_positions()
        weight = np.empty(pos.shape[0])
        InitialRadius_dict = self.parse_key_value_pairs(args.InitialRadius)
        for i in range(len(weight)):
            weight[i] = InitialRadius_dict[self.atoms[i]]*InitialRadius_dict[self.atoms[i]]
        pos = np.hstack((pos, weight.reshape(-1, 1)))
        hc.PDList.from_alpha_filtration(pos,weight=True,save_to=f"pointcloud.pdgm",
                                        save_boundary_map=True)
        pdlist = hc.PDList(f"pointcloud.pdgm")
        pd1 = pdlist.dth_diagram(1)
        pd_data = []
        for birth, death in zip(pd1.births, pd1.deaths):
            pair = pd1.nearest_pair_to(birth, death)
            if pair.lifetime()<1e-5:
                continue
            stable_volume = pair.stable_volume(pair.lifetime()*0.05)
            point_list = stable_volume.boundary_points()
            pd_data+=[birth, death]
            pd_data+=self.seek_point(pos, point_list)
        return pd_data


    
