#ifndef JetTaggingParameters_h
#define JetTaggingParameters_h

#include <iostream>

using namespace std;

namespace JetTagging{

  enum Tagger{
    CSVv2,
    DeepCSV,DeepJet, DeepCvsL, DeepCvsB,
    CvsL, CvsB,
    DeepFlavour_b, DeepFlavour_bb, DeepFlavour_lepb, DeepFlavour_c, DeepFlavour_uds, DeepFlavour_g,
  };
  inline string TaggerToString(Tagger t){
    if(t == CSVv2) return "CSVv2";
    if(t == DeepCSV) return "DeepCSV";
    if(t == DeepJet) return "DeepJet";
    if(t == DeepCvsL) return "DeepCvsL";
    if(t == CvsL) return "CvsL";
    if(t == CvsB) return "CvsB";
    if(t == DeepFlavour_b) return "DeepFlavour_b";
    if(t == DeepFlavour_bb) return "DeepFlavour_bb";
    if(t == DeepFlavour_lepb) return "DeepFlavour_lepb";
    if(t == DeepFlavour_c) return "DeepFlavour_c";
    if(t == DeepFlavour_uds) return "DeepFlavour_uds";
    if(t == DeepFlavour_g) return "DeepFlavour_g";

    cerr << "[JetTaggingParameters::TaggerToString] Wrong tagger enum : " << t << endl;
    exit(ENODATA);

    return "None";
  }
  inline Tagger StringToTagger(string s){
    if(s == "CSVv2") return CSVv2;
    if(s == "DeepCSV") return DeepCSV;
    if(s == "DeepJet") return DeepJet;
    if(s == "DeepCvsL") return DeepCvsL;
    if(s == "CvsL") return CvsL;
    if(s == "CvsB") return CvsB;
    if(s == "DeepFlavour_b") return DeepFlavour_b;
    if(s == "DeepFlavour_bb") return DeepFlavour_bb;
    if(s == "DeepFlavour_lepb") return DeepFlavour_lepb;
    if(s == "DeepFlavour_c") return DeepFlavour_c;
    if(s == "DeepFlavour_uds") return DeepFlavour_uds;
    if(s == "DeepFlavour_g") return DeepFlavour_g;

    cerr << "[JetTaggingParameters:: StringToTagger] Wrong tagger string : " << s << endl;
    exit(ENODATA);

    return CSVv2;
  }

  enum WP{
    Loose, Medium, Tight
  };

  inline string WPToString(WP wp){
    if(wp == Loose) return "Loose";
    if(wp == Medium) return "Medium";
    if(wp == Tight) return "Tight";

    cerr << "[JetTaggingParameters::WPToString] Wrong WP enum : " << wp << endl;
    exit(ENODATA);

    return "None";
  }

  enum MeasurmentType{
    incl, mujets, comb, iterativefit
  };

  inline string MeasurmentTypeToString(MeasurmentType mt){
    if(mt == incl) return "incl";
    if(mt == mujets) return "mujets";
    if(mt == comb) return "comb";
    if(mt == iterativefit) return "iterativefit";

    cerr << "[JetTaggingParameters::MeasurmentTypeToString] Wrong MeasurmentType enum : " << mt << endl;
    exit(ENODATA);

    return "None";

  }

  class Parameters{
  public:
    inline Parameters(){

    }
    inline Parameters(Tagger t, WP wp, MeasurmentType mt_l, MeasurmentType mt_h){
      j_Tagger = t;
      j_WP = wp;
      j_MeasurmentType_Light = mt_l;
      j_MeasurmentType_Heavy = mt_h;
    }

    Tagger j_Tagger;
    WP j_WP;
    MeasurmentType j_MeasurmentType_Light, j_MeasurmentType_Heavy;

    inline void Print(){
      printf("(%s,%s,%s,%s)\n",TaggerToString(j_Tagger).c_str(),WPToString(j_WP).c_str(),MeasurmentTypeToString(j_MeasurmentType_Light).c_str(),MeasurmentTypeToString(j_MeasurmentType_Heavy).c_str());
    }

  };

} // namespace

#endif
