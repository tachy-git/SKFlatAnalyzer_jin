import numpy as np
import torch
from torch_geometric.data import Data
from ROOT import TLorentzVector


def getEdgeIndices(nodeList, k):
    edgeIndex = []
    edgeAttribute = []
    for i, node in enumerate(nodeList):
        distances = {}
        for j, other in enumerate(nodeList):
            if node is other:    # same node
                continue

            thisNode = TLorentzVector()
            otherNode = TLorentzVector()
            thisNode.SetPxPyPzE(node[1], node[2], node[3], node[0])
            otherNode.SetPxPyPzE(other[1], other[2], other[3], other[0])
            distances[j] = thisNode.DeltaR(otherNode)
        distances = dict(sorted(distances.items(), key=lambda item: item[1]))
        for idx in list(distances.keys())[:k]:  # k-nearest node indices
            edgeIndex.append([i, idx])
            edgeAttribute.append([distances[idx]])

    return torch.tensor(edgeIndex, dtype=torch.long), torch.tensor(edgeAttribute, dtype=torch.float)

def evtToGraph(nodeList, y, k=4):
    x = torch.tensor(nodeList, dtype=torch.float)
    edgeIndex, edgeAttribute = getEdgeIndices(nodeList, k=k)
    graph = Data(x=x, y=y,
                 edge_index=edgeIndex.t().contiguous(),
                 edge_attribute=edgeAttribute)
    return graph

def predictProba(model, x, edgeIndex):
    model.eval()
    with torch.no_grad():
        out = model(x, edgeIndex)
        proba = out.numpy()[0][1]
    
    return proba
