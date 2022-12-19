from ROOT import TLorentzVector

class NodeParticle(TLorentzVector):
    def __init__(self):
        TLorentzVector.__init__(self)
        self.charge = 0
        self.isMuon = False
        self.isElectron = False
        self.isJet = False
        self.btagScore = 0.

    def IsMuon(self):
        return self.isMuon

    def IsElectron(self):
        return self.isElectron

    def IsJet(self):
        return self.isJet

    def Charge(self):
        return self.charge

    def BtagScore(self):
        return self.btagScore
