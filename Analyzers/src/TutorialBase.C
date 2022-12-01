#include "TutorialBase.h"

//==== Constructor and Destructor
TutorialBase::TutorialBase() {}
TutorialBase::~TutorialBase() {}

//==== Initialize variables
// NOTE: Every local varaible declared in executeEvent() will be initialized for every event, 
// but some variables do not have to be reinitialized for every event.
// for example, we will use the sample trigger and trigger safe pt cut throughout all events.
// For these varibles(which are called global variables) are declared in the TutorialBase.h
// and initialized in initializeAnalyzer step
void TutorialBase::initializeAnalyzer(){

	//==== Dimuon Z-peak events with two muon IDs
	// One can define customized Muon IDs in DataFormat/src/Muon.C
	// For this tutorial, let's use POG based IDs
	// https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonSelection
	MuonIDs = {
		"POGMedium", 
		"POGTight"
	};

	//==== Trigger settings
	// In this tutorial, we will use HLT_IsoMu27_v (HighLevelTrigger_IsolatedMuon ptcut 27)
	// which is the recommended singlemuon tirgger for year 2017
	// https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2017
	// NOTE: for each dataset, you should use corresponding triggers
	// for example, SingleMuon -> HLT_IsoMu27_v
	//              DoubleMuon -> HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v
	//  (HighLevelTrigger_MuonPtCut17GeV_TrackIsolationVeryVeryLoose_
	//  MuonPtCut8GeV_TrakIsolationVeryVeryLoose_DeltaZ_MassGreaterThan3.8GeV)
	if (DataYear == 2017) {
		IsoMuTriggerName = "HLT_IsoMu27_v";
		TriggerSafePtCut = 29.;
	}
	else {
		cerr << "[TutorialBase::InitializeAnalyzer] This tutorial is for year 2017" << endl;
		exit(EXIT_FAILURE);
	}

	cout << "[TutorialBase::initializeAnalyzer] IsoMuTriggerName = " << IsoMuTriggerName << endl;
	cout << "[TutorialBase::initializeAnalyzer] TriggerSafePtCut = " << TriggerSafePtCut << endl;

}

void TutorialBase::executeEvent(){

	//==== *IMPORTANT TO SAVE CPU TIME*
	//==== Every GetMuon() function first collect ALL MINIAOD muons with GetAllMuons() 
	//==== and then check ID booleans
	AllMuons = GetAllMuons();
	AllJets = GetAllJets();

	//==== Select good quality muons using POG muon ID
	//==== Loop over muon IDs
	for (const auto& MuonID: MuonIDs) {

		//==== No cut
		//==== FillHist fuction is defined in Analyzers/src/AnalyzerCore.C
		//==== void FillHist(histkey, value, weight, nbins, left, right)
		FillHist(MuonID+"/Cutflow", 0., 1., 10, 0., 10.);

		//==== MET (Missing transverse energy) Filter
		if (!PassMETFilter()) return;
		FillHist(MuonID+"/Cutflow", 1., 1., 10, 0., 10.);

		Event ev = GetEvent();
		const Particle METv = ev.GetMETVector();

		//==== Trigger
		//==== see Dataformat/src/Event.C
		if (! (ev.PassTrigger(IsoMuTriggerName))) return;
		FillHist(MuonID+"/Cutflow", 2., 1., 10, 0., 10.);

		//==== muons and jets for this event are one the device memory
		//==== select muons and jets for current ID
		//==== using muons with POG Medium/Tight ID, Jets with tight ID
		//==== select muons and jets inside the detector acceptance
		vector<Muon> muons = SelectMuons(AllMuons, MuonID, 20., 2.4);
		vector<Jet> jets = SelectJets(AllJets, "tight", 30., 2.4);

		//==== sort in Pt order
		std::sort(muons.begin(), muons.end(), PtComparing);
		std::sort(jets.begin(), jets.end(), PtComparing);

		//==== b tagging
		//==== in this example, we will use DeepJet algorithm
		//==== see https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation
		vector<Jet> bjets, jets_nonb;
		for (const auto& jet: jets) {
			const double this_discr = jet.GetTaggerResult(JetTagging::DeepJet);
			if (this_discr > mcCorr->GetJetTaggingCutValue(JetTagging::DeepJet, JetTagging::Medium))
				bjets.emplace_back(jet);
			else
				jets_nonb.emplace_back(jet);
		}
		
		//==== Event selection
		//==== Problem 1. Complete the cutflow!
		//==== OS dimuon
		if (muons.size() != 2) return;
		if (muons.at(0).Charge() + muons.at(1).Charge() != 0) return;

		//==== leading muon should pass trigger-safe pt cut
		if (! (muons.at(0).Pt() > TriggerSafePtCut)) return;

		//==== On-Z
		//==== Reconstruct Z for two muons
		//==== see Dataformat/src/Particle.C
		//==== the mass of dimuon do not have to be exact the mass of Z
		//==== due to Heisenberg's uncertainty principle
		//==== therefore, we give a mass window within 15 GeV
		Particle ZCand = muons.at(0) + muons.at(1);
		const double mZ = 91.2; // GeV
		if (! (fabs(mZ - ZCand.M()) < 15.)) return;

		//==== For MC, the events should be normalized to the xsec * Lumi
		double weight = 1.;
		if (! IsDATA) { // MC
			//==== weight_norm_1invpb is set to be event weight normalized to 1/pb
			//==== to get the full luminosity, just pass "Full"
            //==== MCweight is normalized to 1 /pb.
            const double gen_weight = MCweight();
			const double lumi = ev.GetTriggerLumi("Full");
			weight *= gen_weight*lumi;
		}

		//==== Now fill histograms
		//==== see Analyzers/src/AnalyzerCore.C
		FillHist(MuonID+"/ZCand/mass", ZCand.M(), weight, 40, 70., 110.);
		FillHist(MuonID+"/ZCand/pt", ZCand.Pt(), weight, 300, 0., 300.);
		FillHist(MuonID+"/ZCand/eta", ZCand.Eta(), weight, 100, -5., 5.);
		FillHist(MuonID+"/ZCand/phi", ZCand.Phi(), weight, 64, -3.2, 3.2);
		FillHist(MuonID+"/MET", METv.Pt(), weight, 300, 0., 300.);
		FillHist(MuonID+"/bjets/size", bjets.size(), weight, 20, 0., 20.);
		//==== Problem 2. Let's see the kinematic distributions of decayed products. 
		//==== Fill histograms for the leading & subleading muons.
		
		//==== Problem 3. After you compared the data and MCs, 
		//==== you can check many kinematic distributions of DY(Drell-Yan) and other sources.
		//==== Some variables can be used to suppress non-DY backgrounds (such as ttbar process)
		//==== implement additional cuts to suppress backgrounds further.
	}
}	

