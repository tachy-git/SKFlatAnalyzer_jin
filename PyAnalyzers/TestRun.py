from ROOT import gSystem
from ROOT import TLorentzVector
from ROOT import TriLeptonBase
gSystem.Load("/cvmfs/cms.cern.ch/slc7_amd64_gcc900/external/lhapdf/6.2.3/lib/libLHAPDF.so")

from math import sin
import torch
from torch_geometric.data import Data
from models import ParticleNet

def getEdgeIndices(nodeList, k):
    edgeIndex = []
    edgeAttribute = []
    for i, node in enumerate(nodeList):
        distances = {}
        for j, neighbor in enumerate(nodeList):
            # avoid same node
            if node is neighbor:
                continue
            thisPart = TLorentzVector()
            neighborPart = TLorentzVector()
            thisPart.SetPxPyPzE(node[1], node[2], node[3], node[0])
            neighborPart.SetPxPyPzE(neighbor[1], neighbor[2], neighbor[3], neighbor[0])
            distances[j] = thisPart.DeltaR(neighborPart)
        distances = dict(sorted(distances.items(), key=lambda item: item[1]))
        for n in list(distances.keys())[:k]:
            edgeIndex.append([i, n])
            edgeAttribute.append([distances[n]])

    return (torch.tensor(edgeIndex, dtype=torch.long), torch.tensor(edgeAttribute, dtype=torch.float))

def evtToGraph(nodeList, y, k=4):
    x = torch.tensor(nodeList, dtype=torch.float)
    edgeIndex, edgeAttribute = getEdgeIndices(nodeList, k=k)
    data = Data(x=x,
                y=y,
                edge_index=edgeIndex.t().contiguous(),
                edge_attribute=edgeAttribute)
    return data


class NodeParticle(TLorentzVector):
    def __init__(self):
        TLorentzVector.__init__(self)
        self.charge = 0
        self.is_muon = False
        self.is_electron = False
        self.is_jet = False
        self.btagScore = 0.

    def IsMuon(self):
        return self.is_muon

    def IsElectron(self):
        return self.is_electron

    def IsJet(self):
        return self.is_jet

    def Charge(self):
        return self.charge

    def BtagScore(self):
        return self.btagScore


class TestRun(TriLeptonBase):
    def __init__(self):
        super().__init__()
        self.model = ParticleNet(num_features=9, num_classes=2)
        path = "/home/choij/workspace/ChargedHiggsAnalysis/triLepRegion/full/Combine__/MHc-130_MA-90_vs_ttX/models/ParticleNet_Adam_initLR-0p02_CyclicLR.pt"
        self.model.load_state_dict(torch.load(path, map_location=torch.device("cpu")))

    def predictProba(self, x, edgeIndex):
        self.model.eval()
        with torch.no_grad():
            out = self.model(x, edgeIndex)
            proba = out.numpy()[0][1]
        return proba

    # define executeEvent function
    def executeEvent(self):
        if not super().PassMETFilter(): return None
        ev = super().GetEvent()
        truthColl = super().GetGens()
        rawMuons = super().GetAllMuons()
        rawElectrons = super().GetAllElectrons()
        rawJets = super().GetAllJets()
        METv = ev.GetMETVector()

        #### Object definition
        vetoMuons = super().SelectMuons(rawMuons, super().MuonIDs[2], 10., 2.4)
        looseMuons = super().SelectMuons(vetoMuons, super().MuonIDs[1], 10., 2.4)
        tightMuons = super().SelectMuons(looseMuons, super().MuonIDs[0], 10., 2.4)
        vetoElectrons = super().SelectElectrons(rawElectrons, super().ElectronIDs[2], 10., 2.5)
        looseElectrons = super().SelectElectrons(vetoElectrons, super().ElectronIDs[1], 10., 2.5)
        tightElectrons = super().SelectElectrons(looseElectrons, super().ElectronIDs[0], 10., 2.5)
        jets = super().SelectJets(rawJets, "tight", 20., 2.4)
        jets = super().JetsVetoLeptonInside(jets, vetoElectrons, vetoMuons, 0.4)
        bjets = []
        for jet in jets:
            score = jet.GetTaggerResult(3)                      # DeepJet
            wp = super().mcCorr.GetJetTaggingCutValue(3, 1)     # DeepJet Medium
            if score > wp:
                bjets.append(jet)

        # sort objects
        vetoMuons = list(sorted(vetoMuons, key=lambda x: x.Pt(), reverse=True))
        looseMuons = list(sorted(looseMuons, key=lambda x: x.Pt(), reverse=True))
        tightMuons = list(sorted(tightMuons, key=lambda x: x.Pt(), reverse=True))
        vetoElectrons = list(sorted(vetoElectrons, key=lambda x: x.Pt(), reverse=True))
        looseElectrons = list(sorted(looseElectrons, key=lambda x: x.Pt(), reverse=True))
        tightElectrons = list(sorted(tightElectrons, key=lambda x: x.Pt(), reverse=True))
        jets = list(sorted(jets, key=lambda x: x.Pt(), reverse=True))
        bjets = list(sorted(bjets, key = lambda x: x.Pt(), reverse=True))
        
        #### event selection
        is3Mu = (len(looseMuons) == 3 and len(vetoMuons) == 3 and \
                len(looseElectrons) == 0 and len(vetoElectrons) == 0)
        is1E2Mu = len(looseMuons) == 2 and len(vetoMuons) == 2 and \
                  len(looseElectrons) == 1 and len(vetoElectrons) == 1
        if not (is3Mu or is1E2Mu): return None
        ## for fake estimation in data
        if not super().IsDATA:
            if not len(tightMuons) == len(looseMuons): return None
            if not len(tightElectrons) == len(looseElectrons): return None
    
        channel = ""
        ## 1E2Mu baseline
        ## 1. pass EMuTriggers
        ## 2. Exact 2 tight muons and 1 tight electron, no additional lepton
        ## 3. Exists OS muon pair with mass > 12 GeV
        if is1E2Mu:
            if not ev.PassTrigger(super().EMuTriggers): return None
            
        ## 3Mu baseline
        ## 1. pass DblMuTriggers
        ## 2. Exact 3 tight muons, no additional leptons
        ## 3. Exist OS muon pair,
        ## 4. All OS muon pair mass > 12 GeV
        else:
            if not ev.PassTrigger(super().DblMuTriggers): return None
            mu1, mu2, mu3  = tuple(looseMuons)
            if not mu1.Pt() > 20.: return None
            if not mu2.Pt() > 10.: return None
            if not mu3.Pt() > 10.: return None
            if not abs(mu1.Charge()+mu2.Charge()+mu3.Charge()) == 1: return None
            if mu1.Charge() == mu2.Charge():
                pair1 = mu1 + mu3
                pair2 = mu2 + mu3
            elif mu1.Charge() == mu3.Charge():
                pair1 = mu1 + mu2
                pair2 = mu2 + mu3
            else:   # mu2.Charge() == mu3.Charge()
                pair1 = mu1 + mu2
                pair2 = mu1 + mu3
            if not pair1.M() > 12.: return None
            if not pair2.M() > 12.: return None

            # orthogonality of SR and CR done by bjet multiplicity
            if len(bjets) >= 1:
                if len(jets) >= 2: channel = "SR3Mu"
                else: return None
            else:
                mZ = 91.2
                isOnZ = abs(pair1.M() - mZ) < 10. or abs(pair2.M() - mZ) < 10.
                if isOnZ: channel = "ZFake3Mu"
                else:
                    if abs((mu1+mu2+mu3).M() - mZ) < 10.: channel = "ZGamma3Mu"
                    else: return None

        if "3Mu" not in channel: return None
        ## event selection done
        ## make a graph
        particles = []
        for muon in tightMuons:
            node = NodeParticle()
            node.is_muon = True
            node.SetPtEtaPhiM(muon.Pt(), muon.Eta(), muon.Phi(), muon.M())
            node.charge = muon.Charge()
            particles.append(node)
        for ele in tightElectrons:
            node = NodeParticle()
            node.is_electron = True
            node.SetPtEtaPhiM(ele.Pt(), ele.Eta(), ele.Phi(), ele.M())
            node.charge = ele.Charge()
            particles.append(node)
        for jet in jets:
            node = NodeParticle()
            node.is_jet = True
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
                             particle.Charge(), particle.IsMuon(), particle.IsElectron(),
                             particle.IsJet(), particle.BtagScore()])
        data = evtToGraph(nodeList, y=None, k=4)
        score = self.predictProba(data.x, data.edge_index)
        super().FillHist(f"{channel}/score", score, 1., 100, 0., 1.)




if __name__ == "__main__":
    m = TestRun()
    m.SetTreeName("recoTree/SKFlat")
    m.IsDATA = False
    m.MCSample = "TTToHcToWAToMuMu_MHc-130_MA-90"
    m.xsec = 0.015
    m.sumSign = 599702.0
    m.sumW = 3270.46
    m.IsFastSim = False
    m.SetEra("2017")
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_14.root"): exit(1)
    if not m.AddFile("/home/choij/workspace/DATA/SKFlat/Run2UltraLegacy_v3/2017/TTToHcToWAToMuMu_MHc-130_MA-90_MultiLepFilter_TuneCP5_13TeV-madgraph-pythia8/SKFlat_Run2UltraLegacy_v3/220714_084244/0000/SKFlatNtuple_2017_MC_5.root"): exit(1)
    m.SetOutfilePath("hists.root")
    m.Init()
    m.initializeAnalyzer()
    m.initializeAnalyzerTools()
    m.SwitchToTempDir()
    m.Loop()
    m.WriteHist()
