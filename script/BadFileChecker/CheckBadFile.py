import ROOT,os

def CheckBadFile(filepath):
  f = ROOT.TFile(filepath)

  IsZombie = f.IsZombie()
  IsRecovered = f.TestBit(1024)

  return (IsZombie or IsRecovered)

def CheckHasEntries(filepath):
  

  _localfile = open(filepath, "r")
  nevent_tree1=0
  f = ROOT.TFile(filepath)
  _tree  = f.Get("recoTree/SKFlat")
  if _tree:
    nevent_tree1 = nevent_tree1 + _tree.GetEntries()
  f.Close()
  _localfile.close()


  if nevent_tree1 == 0 :
    return False
  else:
    return True
