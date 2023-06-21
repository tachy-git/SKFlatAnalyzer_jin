import os
import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data
from ROOT import TLorentzVector
from itertools import product, combinations
from MLTools.models import SNN, ParticleNetV2
from MLTools.formats import NodeParticle

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

def loadModels(network, channel, signals, backgrounds):
    models = {}

    for sig, bkg in product(signals, backgrounds):
        csv = pd.read_csv(f"{os.environ['DATA_DIR']}/Classifiers/{network}/{channel}/summary/info-{sig}_vs_{bkg}.txt",
                              sep=",\s",
                              engine="python",
                              header=None).transpose()
        modelPath = f"{os.environ['DATA_DIR']}/Classifiers/{network}/{channel}/models/{sig}_vs_{bkg}.pt" 
        # note: dropout will not work in model.eval()
        # so in principal, it should have no effect in inference time
        if network == "DenseNeuralNet":
            num_nodes = int(csv[0][3])
            if channel == "Skim1E2Mu": model = SNN(41, 2, num_nodes=num_nodes, dropout_p=0.4)
            elif channel == "Skim3Mu": model = SNN(47, 2, num_nodes=num_nodes, dropout_p=0.4)
            else:
                print(f"Wrong channel {channel}")
                raise(ValueError)
        elif network == "GraphNeuralNet":
            num_hidden = int(csv[0][3])
            if channel == "Skim1E2Mu": model = ParticleNetV2(9, 4, 2, num_hidden=num_hidden, dropout_p=0.4)
            elif channel == "Skim3Mu": model = ParticleNetV2(9, 6, 2, num_hidden=num_hidden, dropout_p=0.4)
            else:
                print(f"Wrong channel {channel}")
                raise(ValueError)
        else:
            print(f"Wrong network {network}")
            raise(ValueError)
        model.load_state_dict(torch.load(modelPath, map_location=torch.device("cpu")))
        model.eval()
        models[f"{sig}_vs_{bkg}"] = model
    return models

def getDenseInput(muons, electrons, jets, bjets, METv):
    inputs = []
    dRjets = []
    for idx1, idx2 in combinations(range(jets.size()), 2):
        dRjets.append(jets[idx1].DeltaR(jets[idx2]))
    dRjets = sum(dRjets) / len(dRjets)

    if electrons.size() == 1 and muons.size() == 2:
        mu1 = muons.at(0)
        mu2 = muons.at(1)
        ele = electrons.at(0)
        j1 = jets.at(0)
        j2 = jets.at(1)
        MT = (ele+METv).Mt()

        inputs.append(mu1.E())
        inputs.append(mu1.Px())
        inputs.append(mu1.Py())
        inputs.append(mu1.Pz())
        inputs.append(mu1.Charge())
        inputs.append(mu2.E())
        inputs.append(mu2.Px())
        inputs.append(mu2.Py())
        inputs.append(mu2.Pz())
        inputs.append(mu2.Charge())
        inputs.append(ele.E())
        inputs.append(ele.Px())
        inputs.append(ele.Py())
        inputs.append(ele.Pz())
        inputs.append(ele.Charge())
        inputs.append(j1.E())
        inputs.append(j1.Px())
        inputs.append(j1.Py())
        inputs.append(j1.Pz())
        inputs.append(j1.Charge())
        inputs.append(j1.GetTaggerResult(3))
        inputs.append(j2.E())
        inputs.append(j2.Px())
        inputs.append(j2.Py())
        inputs.append(j2.Pz())
        inputs.append(j2.Charge())
        inputs.append(j2.GetTaggerResult(3))
        inputs.append(mu1.DeltaR(mu2))
        inputs.append(mu1.DeltaR(ele))
        inputs.append(mu2.DeltaR(ele))
        inputs.append(j1.DeltaR(j2))
        inputs.append(sum([j.Pt() for j in jets]))
        inputs.append(sum([l.Pt() for l in electrons+muons]))
        inputs.append(MT)
        inputs.append(METv.Pt())
        inputs.append(jets.size())
        inputs.append(bjets.size())
        inputs.append(dRjets)
        inputs.append(sum([j.GetTaggerResult(3) for j in jets]) / jets.size())

    elif muons.size() == 3 and electrons.size() == 0:
        mu1 = muons.at(0)
        mu2 = muons.at(1)
        mu3 = muons.at(2)
        j1 = jets.at(0)
        j2 = jets.at(1)
        MT1 = (mu1+METv).Mt()
        MT2 = (mu2+METv).Mt()
        MT3 = (mu3+METv).Mt()

        inputs.append(mu1.E())
        inputs.append(mu1.Px())
        inputs.append(mu1.Py())
        inputs.append(mu1.Pz())
        inputs.append(mu1.Charge())
        inputs.append(mu2.E())
        inputs.append(mu2.Px())
        inputs.append(mu2.Py())
        inputs.append(mu2.Pz())
        inputs.append(mu2.Charge())
        inputs.append(mu3.E())
        inputs.append(mu3.Px())
        inputs.append(mu3.Py())
        inputs.append(mu3.Pz())
        inputs.append(mu3.Charge())
        inputs.append(j1.E())
        inputs.append(j1.Px())
        inputs.append(j1.Py())
        inputs.append(j1.Pz())
        inputs.append(j1.Charge())
        inputs.append(j1.GetTaggerResult(3))
        inputs.append(j2.E())
        inputs.append(j2.Px())
        inputs.append(j2.Py())
        inputs.append(j2.Pz())
        inputs.append(j2.Charge())
        inputs.append(j2.GetTaggerResult(3))
        inputs.append(mu1.DeltaR(mu2))
        inputs.append(mu1.DeltaR(mu3))
        inputs.append(mu2.DeltaR(mu3))
        inputs.append(j1.DeltaR(mu1))
        inputs.append(j1.DeltaR(mu2))
        inputs.append(j1.DeltaR(mu3))
        inputs.append(j2.DeltaR(mu1))
        inputs.append(j2.DeltaR(mu2))
        inputs.append(j2.DeltaR(mu3))
        inputs.append(j1.DeltaR(j2))
        inputs.append(sum([j.Pt() for j in jets]))
        inputs.append(sum([l.Pt() for l in muons]))
        inputs.append(MT1)
        inputs.append(MT2)
        inputs.append(MT3)
        inputs.append(METv.Pt())
        inputs.append(jets.size())
        inputs.append(bjets.size())
        inputs.append(dRjets)
        inputs.append(sum([j.GetTaggerResult(3) for j in jets])/jets.size())
    else:
        print(f"Wrong number of leptons ({electrons.size()}, {muons.size()})")
        raise(ValueError)
    
    return torch.FloatTensor([inputs])

def getDenseScore(self, model, data):
    with torch.no_grad():
        out = self.model(data)
    return out.numpy()[0][1]

def getGraphInput(muons, electrons, jets, bjets, METv):
    particles = []
    for muon in muons:
        node = NodeParticle()
        node.isMuon = True
        node.SetPtEtaPhiM(muon.Pt(), muon.Eta(), muon.Phi(), muon.M())
        node.charge = muon.Charge()
        particles.append(node)
    for ele in electrons:
        node = NodeParticle()
        node.isElectron = True
        node.SetPtEtaPhiM(ele.Pt(), ele.Eta(), ele.Phi(), ele.M())
        node.charge = ele.Charge()
        particles.append(node)
    for jet in jets:
        node = NodeParticle()
        node.isJet = True
        node.SetPtEtaPhiM(jet.Pt(), jet.Eta(), jet.Phi(), jet.M())
        node.charge = jet.Charge()
        node.btagScore = jet.GetTaggerResult(3)
        particles.append(node)
    missing = NodeParticle()
    missing.SetPtEtaPhiM(METv.Pt(), 0., METv.Phi(), 0.)
    particles.append(missing)
    nodeList = []
    for particle in particles:
        nodeList.append([particle.E(), particle.Px(), particle.Py(), particle.Pz(),
                         particle.Charge(), particle.BtagScore(),
                         particle.IsMuon(), particle.IsElectron(), particle.IsJet()])
    data = evtToGraph(nodeList, y=None, k=4)
    
    if muons.size() == 3:
        MT1 = (muons.at(0)+METv).Mt()
        MT2 = (muons.at(1)+METv).Mt()
        MT3 = (muons.at(2)+METv).Mt()
        data.graph_input = torch.tensor([[jets.size(), bjets.size(), METv.Pt(), MT1, MT2, MT3]], dtype=torch.float)
    elif electrons.size() == 1 and muons.size() == 2:
        MT = (electrons.at(0)+METv).Mt()
        data.graph_input = torch.tensor([[jets.size(), bjets.size(), METv.Pt(), MT]], dtype=torch.float)
    else:
        print(f"Wrong size of muons {muons.size()} and electrons {electrons.size()}")
        raise(ValueError)
    return data

def getGraphScore(model, data):
    with torch.no_grad():
        out = model(data.x, data.edge_index, data.graph_input)
    return out.numpy()[0][1]

