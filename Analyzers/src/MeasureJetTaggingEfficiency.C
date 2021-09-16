#include "MeasureJetTaggingEfficiency.h"

MeasureJetTaggingEfficiency::MeasureJetTaggingEfficiency(){
}

void MeasureJetTaggingEfficiency::initializeAnalyzer(){

  TString datapath = getenv("DATA_DIR");
  TString btagpath = datapath+"/"+DataEra+"/BTag/";

  Taggers.clear();
  WPs.clear();
  CutValues.clear();

  ifstream in_tagger(btagpath+"/CutValues.txt");
  string btaggerline;
  while(getline(in_tagger,btaggerline)){
    std::istringstream is_tag( btaggerline );
    TString tstring_taggerline = btaggerline;
    if(tstring_taggerline.Contains("#")) continue;
    TString a;
    string b,c;
    float d;

    is_tag >> a; // ERA
    is_tag >> b; // TAGGER
    is_tag >> c; // WP
    is_tag >> d; // cut value

    if(a == DataEra){
      Taggers.push_back(b);
      WPs.push_back(c);
      CutValues.push_back(d);
    }
  }// end of taggermap loop

  cout << "[MeasureJetTaggingEfficiency::initializeAnalyzer()] What to measure :" << endl;
  cout << "[MeasureJetTaggingEfficiency::initializeAnalyzer()] Tagger\tWP\tCutValue" << endl;
  for(unsigned i_m=0; i_m<Taggers.size(); i_m++){

    string Tagger = Taggers.at(i_m);
    string WP = WPs.at(i_m);
    double CutValue = CutValues.at(i_m);

    cout << "[MeasureJetTaggingEfficiency::initializeAnalyzer()] " << Tagger << "\t" << WP << "\t" << CutValue << endl;
  }
}

MeasureJetTaggingEfficiency::~MeasureJetTaggingEfficiency(){
}

void MeasureJetTaggingEfficiency::executeEvent(){

  //========================
  //==== MET Filter
  //========================

  if(!PassMETFilter()) return;

  Event ev = GetEvent();
  Particle METv = ev.GetMETVector();

  vector<Jet> jets = GetJets("tightLepVeto", 20., 2.5);

  vector<double> vec_etabins = {-2.5, -2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5};
  vector<double> vec_ptbins = {20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 75., 80, 85., 90., 95., 100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 250., 300., 350., 400., 500., 600., 1000.};
  double PtMax = vec_ptbins.at( vec_ptbins.size()-1 );

  const int NEtaBin = vec_etabins.size()-1;
  const int NPtBin = vec_ptbins.size()-1;

  double etabins[NEtaBin+1];
  for(int i=0; i<NEtaBin+1; i++) etabins[i] = vec_etabins.at(i);
  double ptbins[NPtBin+1];
  for(int i=0; i<NPtBin+1; i++) ptbins[i] = vec_ptbins.at(i);

  //==== code to measure btag efficiencies in TT MC
  //==== Reference : https://github.com/rappoccio/usercode/blob/Dev_53x/EDSHyFT/plugins/BTaggingEffAnalyzer.cc
  for(unsigned int ij = 0 ; ij < jets.size(); ij++){

    TString flav= "B";
    if(fabs(jets.at(ij).hadronFlavour()) == 4) flav= "C";
    if(fabs(jets.at(ij).hadronFlavour()) == 0) flav= "Light";

    double this_Eta = jets.at(ij).Eta();
    double this_Pt = jets.at(ij).Pt()<PtMax ? jets.at(ij).Pt() : PtMax-1; // put overflows in the last bin

    //==== First, fill the denominator
    FillHist("Jet_"+DataEra+"_eff_"+flav+"_denom", this_Eta, this_Pt, ev.MCweight(), NEtaBin, etabins, NPtBin, ptbins);

    //==== Now looping over (tagger,working point)
    for(unsigned i_m=0; i_m<Taggers.size(); i_m++){ 

      string Tagger = Taggers.at(i_m);
      string WP = WPs.at(i_m);
      double CutValue = CutValues.at(i_m);

      double this_taggerresult = jets.at(ij).GetTaggerResult( JetTagging::StringToTagger(Tagger) );

      if(this_taggerresult>CutValue){
        FillHist("Jet_"+DataEra+"_"+Tagger+"_"+WP+"_eff_"+flav+"_num", this_Eta, this_Pt, ev.MCweight(), NEtaBin, etabins, NPtBin, ptbins);
      }
    } // END Loop (tagger,working point)
  } // END Loop jet
}
