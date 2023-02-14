import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Sequential, Linear, ReLU, Dropout, BatchNorm1d
from torch_geometric.nn import global_mean_pool, global_max_pool, knn_graph
from torch_geometric.nn import GraphNorm
from torch_geometric.nn import MessagePassing


class SNN(nn.Module):
    def __init__(self, nFeatures, nClasses):
        super(SNN, self).__init__()
        # Lecun init is default for pytorch
        self.bn = nn.BatchNorm1d(nFeatures)
        self.dense1 = nn.Linear(nFeatures, 64, bias=True)
        self.dense2 = nn.Linear(64, 128, bias=True)
        self.dense3 = nn.Linear(128, 256, bias=True)
        self.dense4 = nn.Linear(256, 128, bias=True)
        self.dense5 = nn.Linear(128, 64, bias=True)
        self.dense6 = nn.Linear(64, nClasses, bias=True)

    def forward(self, x):
        x = F.selu(self.dense1(x))
        x = F.alpha_dropout(x, p=0.5, training=self.training)
        x = F.selu(self.dense2(x))
        x = F.alpha_dropout(x, p=0.5, training=self.training)
        x = F.selu(self.dense3(x))
        x = F.alpha_dropout(x, p=0.5, training=self.training)
        x = F.selu(self.dense4(x))
        x = F.alpha_dropout(x, p=0.5, training=self.training)
        x = F.selu(self.dense5(x))
        x = F.alpha_dropout(x, p=0.5, training=self.training)
        out = F.softmax(self.dense6(x), dim=1)

        return out

class SNNLite(nn.Module):
    def __init__(self, nFeatures, nClasses):
        super(SNNLite, self).__init__()
        # Lecun init is default for pytorch
        self.bn = nn.BatchNorm1d(nFeatures)
        self.dense1 = nn.Linear(nFeatures, 64, bias=True)
        self.dense2 = nn.Linear(64, 128, bias=True)
        self.dense3 = nn.Linear(128, 64, bias=True)
        self.dense4 = nn.Linear(64, nClasses, bias=True)

    def forward(self, x):
        x = F.selu(self.dense1(x))
        x = F.alpha_dropout(x, p=0.5, training=self.training)
        x = F.selu(self.dense2(x))
        x = F.alpha_dropout(x, p=0.5, training=self.training)
        x = F.selu(self.dense3(x))
        x = F.alpha_dropout(x, p=0.5, training=self.training)
        out = F.softmax(self.dense4(x), dim=1)

        return out

class EdgeConv(MessagePassing):
    def __init__(self, in_channels, out_channels, dropout_p):
        super().__init__(aggr='mean')
        self.mlp = Sequential(
                Linear(2*in_channels, out_channels), ReLU(), BatchNorm1d(out_channels), Dropout(dropout_p),
                Linear(out_channels, out_channels), ReLU(), BatchNorm1d(out_channels), Dropout(dropout_p),
                Linear(out_channels, out_channels), ReLU(), BatchNorm1d(out_channels), Dropout(dropout_p)
                )

    def forward(self, x, edge_index, batch=None):
        return self.propagate(edge_index, x=x, batch=batch)

    def message(self, x_i, x_j):
        tmp = torch.cat([x_i, x_j - x_i], dim=1)
        return self.mlp(tmp)


class DynamicEdgeConv(EdgeConv):
    def __init__(self, in_channels, out_channels, dropout_p, k=4):
        super().__init__(in_channels, out_channels, dropout_p)
        self.shortcut = Sequential(Linear(in_channels, out_channels), BatchNorm1d(out_channels))
        self.k = k

    def forward(self, x, edge_index=None, batch=None):
        if edge_index is None:
            edge_index = knn_graph(
                x, self.k, batch, loop=False, flow=self.flow)
        out = super().forward(x, edge_index, batch=batch)
        out += self.shortcut(x)
        return out


class ParticleNet(torch.nn.Module):
    def __init__(self, num_features, num_classes, dropout_p, readout):
        super(ParticleNet, self).__init__()
        self.gn0 = GraphNorm(num_features)
        self.conv1 = DynamicEdgeConv(num_features, 64, dropout_p)
        self.conv2 = DynamicEdgeConv(64, 128, dropout_p)
        self.conv3 = DynamicEdgeConv(128, 128, dropout_p)
        self.dense1 = Linear(128, 64)
        self.bn1 = BatchNorm1d(64)
        self.dense2 = Linear(64, 64)
        self.bn2 = BatchNorm1d(64)
        self.output = Linear(64, num_classes)
        self.dropout_p = dropout_p
        self.readout = readout

    def forward(self, x, edge_index, batch=None):
        # Convolution layers
        x = self.gn0(x, batch=batch)
        x = self.conv1(x, edge_index, batch=batch)
        x = self.conv2(x, batch=batch)
        x = self.conv3(x, batch=batch)
        # readout layers
        if self.readout == "mean":
            x = global_mean_pool(x, batch=batch)
        elif self.readout == "max":
            x = global_max_pool(x, batch=batch)
        else:
            print(f"Wrong readout {readout}")
            exit(1)

        # dense layers
        x = F.relu(self.dense1(x))
        x = self.bn1(x)
        x = F.dropout1d(x, p=self.dropout_p, training=self.training)
        x = F.relu(self.dense2(x))
        x = self.bn2(x)
        x = F.dropout1d(x, p=self.dropout_p, training=self.training)
        x = self.output(x)

        return F.softmax(x, dim=1)


class ParticleNetLite(torch.nn.Module):
    def __init__(self, num_features, num_classes, dropout_p, readout):
        super(ParticleNetLite, self).__init__()
        self.gn0 = GraphNorm(num_features)
        self.conv1 = DynamicEdgeConv(num_features, 32, dropout_p)
        self.conv2 = DynamicEdgeConv(32, 64, dropout_p)
        self.conv3 = DynamicEdgeConv(64, 64, dropout_p)
        self.dense1 = Linear(64, 32)
        self.bn1 = BatchNorm1d(32)
        self.dense2 = Linear(32, 32)
        self.bn2 = BatchNorm1d(32)
        self.output = Linear(32, num_classes)
        self.dropout_p = dropout_p
        self.readout = readout

    def forward(self, x, edge_index, batch=None):
        # Convolution layers
        x = self.gn0(x, batch=batch)
        x = self.conv1(x, edge_index, batch=batch)
        x = self.conv2(x, batch=batch)
        x = self.conv3(x, batch=batch)
        # readout layers
        if self.readout == "mean":
            x = global_mean_pool(x, batch=batch)
        elif self.readout == "max":
            x = global_max_pool(x, batch=batch)
        else:
            print(f"Wrong readout {readout}")
            exit(1)

        # dense layers
        x = F.relu(self.dense1(x))
        x = self.bn1(x)
        x = F.dropout1d(x, p=self.dropout_p, training=self.training)
        x = F.relu(self.dense2(x))
        x = self.bn2(x)
        x = F.dropout1d(x, p=self.dropout_p, training=self.training)
        x = self.output(x)

        return F.softmax(x, dim=1)
