from operator import le
import pandas as pd

class Connectome:
    def __init__(self, filename, tau=0.5) -> None:
        self.filename = filename
        self.connectome = None
        self.binary_connectome = None
        self.tau = tau

    def read_connectome(self) -> pd.DataFrame:
        self.connectome = pd.read_csv(self.filename, sep=' ', header=None)
        return self.connectome

    def get_binary_connectome(self, tau) -> pd.DataFrame:
        self.binary_connectome = self.connectome.applymap(lambda x: 1 if x > tau else 0)
        return self.binary_connectome

    def get_tau(self):
        return self.tau
    
    def set_tau(self, tau):
        self.tau = tau
        self.get_binary_connectome(tau)

    def check_binary_connectome(self):
        if type(self.binary_connectome) != pd.DataFrame:
            self.get_binary_connectome(self.get_tau())
        
    def get_net_degree_und(self):
        self.check_binary_connectome()
        return self.binary_connectome.sum(axis=1).sum()/2

    def get_connection_density(self):
        num_e = self.get_net_degree_und()
        total = (218**2) / 2 - 218
        return num_e/ (total/(total-1))
    
    def get_network_strength(self):
        return self.connectome.sum().sum()/2

    def get_node_degree_und(self):
        nodes = []
        for i in range(218):
            nodes.append(self.binary_connectome.iloc[i].sum() - self.binary_connectome.iloc[i][i])
        return nodes

    def get_node_degree_mean(self):
        nodes = self.get_node_degree_und()
        return nodes/len(nodes)
    
    def get_node_strength_und(self):
        nodes = []
        for i in range(218):
            nodes.append(self.connectome.iloc[i].sum())
        return nodes
    
    def get_node_strength_mean(self):
        nodes = self.get_node_strength_und()
        return nodes/len(nodes)