import torch
import torch.nn.functional as F
from torch.nn import Sequential, Linear, ReLU, BatchNorm1d
from torch_geometric.nn import global_mean_pool, knn_graph
from torch_geometric.nn import GraphNorm
from torch_geometric.nn import MessagePassing

class EdgeConv(MessagePassing):
    def __init__(self, in_channels, out_channels):
        super().__init__(aggr='mean')
        self.mlp = Sequential(
                Linear(2 * in_channels, out_channels), BatchNorm1d(out_channels), ReLU(),
                Linear(out_channels, out_channels), BatchNorm1d(out_channels), ReLU(),
                Linear(out_channels, out_channels), BatchNorm1d(out_channels), ReLU())

    def forward(self, x, edge_index, batch=None):
        return self.propagate(edge_index, x=x, batch=batch)

    def message(self, x_i, x_j):
        tmp = torch.cat([x_i, x_j - x_i], dim=1)
        return self.mlp(tmp)


class DynamicEdgeConv(EdgeConv):
    def __init__(self, in_channels, out_channels, k=4):
        super().__init__(in_channels, out_channels)
        self.shortcut = Sequential(Linear(in_channels, out_channels), BatchNorm1d(out_channels))
        self.k = k

    def forward(self, x, edge_index=None, batch=None):
        if edge_index is None:
            edge_index = knn_graph(
                x, self.k, batch, loop=False, flow=self.flow)
        out = super().forward(x, edge_index, batch=batch)
        out += self.shortcut(x)
        return F.relu(out)

class ParticleNet(torch.nn.Module):
    def __init__(self, num_features, num_classes):
        super(ParticleNet, self).__init__()
        self.gn0 = GraphNorm(num_features)
        self.conv1 = DynamicEdgeConv(num_features, 64)
        self.gn1 = GraphNorm(64)
        self.conv2 = DynamicEdgeConv(64, 128)
        self.gn2 = GraphNorm(128)
        self.conv3 = DynamicEdgeConv(128, 128)
        self.gn3 = GraphNorm(128)
        self.dense1 = Linear(128, 64)
        self.dense2 = Linear(64, 64)
        self.output = Linear(64, num_classes)

    def forward(self, x, edge_index, batch=None):
        # Convolution layers
        x = self.gn0(x, batch=batch)
        x = self.conv1(x, edge_index, batch=batch)
        x = self.gn1(x, batch=batch)
        x = self.conv2(x, batch=batch)
        x = self.gn2(x, batch=batch)
        x = self.conv3(x, batch=batch)
        x = self.gn3(x, batch=batch)
        # readout layers
        x = global_mean_pool(x, batch=batch)

        # dense layers
        x = F.relu(self.dense1(x))
        x = F.dropout(x, p=0.2)
        x = F.relu(self.dense2(x))
        x = F.dropout(x, p=0.2)
        x = self.output(x)

        return F.softmax(x, dim=1)


class ParticleNetLite(torch.nn.Module):
    def __init__(self, num_features, num_classes):
        super(ParticleNetLite, self).__init__()
        self.gn0 = GraphNorm(num_features)
        self.conv1 = DynamicEdgeConv(num_features, 32)
        self.gn1 = GraphNorm(32)
        self.conv2 = DynamicEdgeConv(32, 64)
        self.gn2 = GraphNorm(64)
        self.conv3 = DynamicEdgeConv(64, 64)
        self.gn3 = GraphNorm(64)
        self.dense1 = Linear(64, 32)
        self.dense2 = Linear(32, 32)
        self.output = Linear(32, num_classes)

    def forward(self, x, edge_index, batch=None):
        # Convolution layers
        x = self.gn0(x, batch=batch)
        x = self.conv1(x, edge_index, batch=batch)
        x = self.gn1(x, batch=batch)
        x = self.conv2(x, batch=batch)
        x = self.gn2(x, batch=batch)
        x = self.conv3(x, batch=batch)
        x = self.gn3(x, batch=batch)
        # readout layers
        x = global_mean_pool(x, batch=batch)

        # dense layers                                                                  
        x = F.relu(self.dense1(x))
        x = F.dropout(x, p=0.2)     
        x = F.relu(self.dense2(x))
        x = F.dropout(x, p=0.2)    
        x = self.output(x)
   
        return F.softmax(x, dim=1)
