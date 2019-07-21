#include "SampleAnalyzer/User/Analyzer/cms_sus_18_002.h"
using namespace MA5;
using namespace std;


double DeltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}

double Mini_Iso_DeltaR(const RecLeptonFormat* Lepton)
{ double deltaR = 0.0;
  double pt = Lepton->momentum().Pt();
  if(pt<50.0) { deltaR = 0.2;}
  if(pt >= 50.0 and pt <= 200.0) { deltaR = 10.0/pt;}
  if(pt> 200.0) { deltaR = 0.05;}
  return deltaR;
}

// Overlap Removal
template<typename T1, typename T2> std::vector<const T1*> Removal(std::vector<const T1*> &v1,
								  std::vector<const T2*> &v2, const double &drmin)
{
  // Determining with objects should be removed
  std::vector<bool> mask(v1.size(),false);
  for (unsigned int j=0;j<v1.size();j++)
    for (unsigned int i=0;i<v2.size();i++)
      if (v2[i]->dr(v1[j]) < drmin)
	{
	  mask[j]=true;
	  break;
	}

  // Building the cleaned container
  std::vector<const T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
    if (!mask[i]) cleaned_v1.push_back(v1[i]);

  return cleaned_v1;
}

template<typename T1> std::vector<const T1*> Removal(std::vector<const T1*> &v1,
						     const double &drmin)
{
  // Determining with objects should be removed (objects are sorted)
  std::vector<bool> mask(v1.size(),false);
  for (unsigned int j=0;j<v1.size();j++)
    for (unsigned int i=j+1;i<v1.size();i++)
      {
	if (mask[i]) continue;
	if (v1[i]->dr(v1[j]) < drmin)
	  {
	    mask[i]=true;
	    continue;
	  }
      }

  // Building the cleaned container
  std::vector<const T1*> cleaned_v1;
  for (unsigned int i=0;i<v1.size();i++)
    if (!mask[i]) cleaned_v1.push_back(v1[i]);

  return cleaned_v1;
}




// -----------------------------------------------------------------------------
// Initialize
// function called one time at the beginning of the analysis
// -----------------------------------------------------------------------------
bool cms_sus_18_002::Initialize(const MA5::Configuration& cfg, const std::map<std::string,std::string>& parameters)
{
  cout << "BEGIN Initialization" << endl;
  // displaying the information


  INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>"<< endmsg;
  INFO << "        <>    Analysis: cms_sus_18_002                                <>"<< endmsg;
  INFO << "        <>    Photon,MET,Multijet,b-jet = 13 TeV, 35.9 fb^-1 luminosity  <>" << endmsg;
  INFO << "        <>    Recast by:                                              <>" << endmsg;
  INFO << "        <>    Contact:                                                <>" << endmsg;
  INFO << "        <>    Based on MadAnalysis 5 v1.6 and above                   <>" << endmsg;
  INFO << "        <>    For more information, see                               <>" << endmsg;
  INFO << "        <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>" << endmsg;

  // initializing selection region

  Manager()->AddRegionSelection("SR1_NJ2-4_NB0_MISSPT270");
  Manager()->AddRegionSelection("SR2_NJ2-4_NB0_MISSPT350");
  Manager()->AddRegionSelection("SR3_NJ2-4_NB0_MISSPT450");
  Manager()->AddRegionSelection("SR4_NJ2-4_NB0_MISSPT750");
  Manager()->AddRegionSelection("SR5_NJ2-4_NB0_MISSPT>750");
  Manager()->AddRegionSelection("SR6_NJ5-6_NB0_MISSPT270");
  Manager()->AddRegionSelection("SR7_NJ5-6_NB0_MISSPT350");
  Manager()->AddRegionSelection("SR8_NJ5-6_NB0_MISSPT450");
  Manager()->AddRegionSelection("SR9_NJ5-6_NB0_MISSPT>450");
  Manager()->AddRegionSelection("SR10_NJ>7_NB0_MISSPT270");
  Manager()->AddRegionSelection("SR11_NJ>7_NB0_MISSPT350");
  Manager()->AddRegionSelection("SR12_NJ>7_NB0_MISSPT450");
  Manager()->AddRegionSelection("SR13_NJ>7_NB0_MISSPT>450");
  Manager()->AddRegionSelection("SR14_NJ2-4_NB1_MISSPT270");
  Manager()->AddRegionSelection("SR15_NJ2-4_NB1_MISSPT350");
  Manager()->AddRegionSelection("SR16_NJ2-4_NB1_MISSPT450");
  Manager()->AddRegionSelection("SR17_NJ2-4_NB1_MISSPT>450");
  Manager()->AddRegionSelection("SR18_NJ5-6_NB1_MISSPT270");
  Manager()->AddRegionSelection("SR19_NJ5-6_NB1_MISSPT350");
  Manager()->AddRegionSelection("SR20_NJ5-6_NB1_MISSPT450");
  Manager()->AddRegionSelection("SR21_NJ5-6_NB1_MISSPT>450");
  Manager()->AddRegionSelection("SR22_NJ>7_NB1_MISSPT270");
  Manager()->AddRegionSelection("SR23_NJ>7_NB1_MISSPT350");
  Manager()->AddRegionSelection("SR24_NJ>7_NB1_MISSPT450");
  Manager()->AddRegionSelection("SR25_NJ>7_NB1_MISSPT>450");
  Manager()->AddRegionSelection("bjets>=1");
  Manager()->AddRegionSelection("bjets=0");


  // initializing the preselection cuts

  Manager()->AddCut("ptgamma>100");
  Manager()->AddCut("vetoelectron&muons");
  Manager()->AddCut("IsoTrackVeto");
  Manager()->AddCut("rec_met>200");
  Manager()->AddCut("njets>=2");
  Manager()->AddCut("STcut");
  Manager()->AddCut("dphicut");
  Manager()->AddCut("bjetcut>1","bjets>=1");
  Manager()->AddCut("bjetcut=0","bjets=0");


  // initializing the selection region cut

  Manager()->AddCut("SR1_NJ2-4_NB0_MISSPT270","SR1_NJ2-4_NB0_MISSPT270");
  Manager()->AddCut("SR2_NJ2-4_NB0_MISSPT350","SR2_NJ2-4_NB0_MISSPT350");
  Manager()->AddCut("SR3_NJ2-4_NB0_MISSPT450","SR3_NJ2-4_NB0_MISSPT450");
  Manager()->AddCut("SR4_NJ2-4_NB0_MISSPT750","SR4_NJ2-4_NB0_MISSPT750");
  Manager()->AddCut("SR5_NJ2-4_NB0_MISSPT>750","SR5_NJ2-4_NB0_MISSPT>750");
  Manager()->AddCut("SR6_NJ5-6_NB0_MISSPT270","SR6_NJ5-6_NB0_MISSPT270");
  Manager()->AddCut("SR7_NJ5-6_NB0_MISSPT350","SR7_NJ5-6_NB0_MISSPT350");
  Manager()->AddCut("SR8_NJ5-6_NB0_MISSPT450","SR8_NJ5-6_NB0_MISSPT450");
  Manager()->AddCut("SR9_NJ5-6_NB0_MISSPT>450","SR9_NJ5-6_NB0_MISSPT>450");
  Manager()->AddCut("SR10_NJ>7_NB0_MISSPT270","SR10_NJ>7_NB0_MISSPT270");
  Manager()->AddCut("SR11_NJ>7_NB0_MISSPT350","SR11_NJ>7_NB0_MISSPT350");
  Manager()->AddCut("SR12_NJ>7_NB0_MISSPT450","SR12_NJ>7_NB0_MISSPT450");
  Manager()->AddCut("SR13_NJ>7_NB0_MISSPT>450","SR13_NJ>7_NB0_MISSPT>450");
  Manager()->AddCut("SR14_NJ2-4_NB1_MISSPT270","SR14_NJ2-4_NB1_MISSPT270");
  Manager()->AddCut("SR15_NJ2-4_NB1_MISSPT350","SR15_NJ2-4_NB1_MISSPT350");
  Manager()->AddCut("SR16_NJ2-4_NB1_MISSPT450","SR16_NJ2-4_NB1_MISSPT450");
  Manager()->AddCut("SR17_NJ2-4_NB1_MISSPT>450","SR17_NJ2-4_NB1_MISSPT>450");
  Manager()->AddCut("SR18_NJ5-6_NB1_MISSPT270","SR18_NJ5-6_NB1_MISSPT270");
  Manager()->AddCut("SR19_NJ5-6_NB1_MISSPT350","SR19_NJ5-6_NB1_MISSPT350");
  Manager()->AddCut("SR20_NJ5-6_NB1_MISSPT450","SR20_NJ5-6_NB1_MISSPT450");
  Manager()->AddCut("SR21_NJ5-6_NB1_MISSPT>450","SR21_NJ5-6_NB1_MISSPT>450");
  Manager()->AddCut("SR22_NJ>7_NB1_MISSPT270","SR22_NJ>7_NB1_MISSPT270");
  Manager()->AddCut("SR23_NJ>7_NB1_MISSPT350","SR23_NJ>7_NB1_MISSPT350");
  Manager()->AddCut("SR24_NJ>7_NB1_MISSPT450","SR24_NJ>7_NB1_MISSPT450");
  Manager()->AddCut("SR25_NJ>7_NB1_MISSPT>450","SR25_NJ>7_NB1_MISSPT>450");

  
  cout << "END   Initialization" << endl;
  return true;
}

// -----------------------------------------------------------------------------
// Finalize
// function called one time at the end of the analysis
// -----------------------------------------------------------------------------
void cms_sus_18_002::Finalize(const SampleFormat& summary, const std::vector<SampleFormat>& files)
{
  cout << "BEGIN Finalization" << endl;
  // saving histos
  cout << "END   Finalization" << endl;
}

// -----------------------------------------------------------------------------
// Execute
// function called each time one event is read
// -----------------------------------------------------------------------------
bool cms_sus_18_002::Execute(SampleFormat& sample, const EventFormat& event)
{
  // ***************************************************************************
  // Example of analysis with generated particles
  // Concerned samples : LHE/STDHEP/HEPMC
  // ***************************************************************************

  //Event weight no need to change this

  double myEventWeight=1.0;
  if(Configuration().IsNoEventWeight()) myEventWeight=1.0 ;
  else if(event.mc()->weight()!=0.) myEventWeight=event.mc()->weight();
  else
    {
      WARNING << "Found one event with a zero weight. Skipping..." << endmsg;
      return false;
    }
  Manager()->InitializeForNewEvent(myEventWeight);
  
  if(event.rec()!=0)
    {
      vector<const RecLeptonFormat*> rec_electrons, rec_muons, isoelec, isomuon;
      vector<const RecPhotonFormat*> isophoton,recphoton;

      // Looking through the reconstructed photon collection
	
      for (MAuint32 i=0;i<event.rec()->photons().size();i++)
	{
	  const RecPhotonFormat *photon = &(event.rec()->photons()[i]);
	  if((photon->momentum().Pt()>100)&&(fabs(photon->momentum().Eta())<2.4))
	    {
	      if((fabs(photon->momentum().Eta())<1.4442)||(fabs(photon->momentum().Eta())>1.556))
		{		  
		  recphoton.push_back(photon);
	      	  double delrmax = 0.3;
		  double allcomp = PHYSICS->Isol->eflow->sumIsolation(photon,event.rec(),delrmax,0.5,IsolationEFlow::ALL_COMPONENTS);
		  double isolation = allcomp/photon->momentum().Pt();
		  if(isolation <= 0.01 ) { isophoton.push_back(photon);}
		}
	    }
	}

      SORTER->sort(isophoton,PTordering);
	
      // Looking through the reconstructed electron collection
      for (MAuint32 i=0;i<event.rec()->electrons().size();i++)
	{
	  const RecLeptonFormat *elec = &(event.rec()->electrons()[i]);
	  
	  if((elec->momentum().Pt()>10)&&(fabs(elec->momentum().Eta())<2.5)) 
	    {
	      rec_electrons.push_back(elec);
	      double isocone = Mini_Iso_DeltaR(elec);
	      double all = PHYSICS->Isol->eflow->sumIsolation(elec,event.rec(),isocone,0.,IsolationEFlow::ALL_COMPONENTS);
	      double isolation = all/elec->momentum().Pt();
	      if(isolation < 0.1) { isoelec.push_back(elec);}    
	    }
	}
	
      // Looking through the reconstructed muon collection
      for (MAuint32 i=0;i<event.rec()->muons().size();i++)
	{
	  const RecLeptonFormat *mu = &(event.rec()->muons()[i]);
	  
	  if((mu->momentum().Pt()>10)&&(fabs(mu->momentum().Eta())<2.4))
	    {
	      rec_muons.push_back(mu);
	      double isocone = Mini_Iso_DeltaR(mu);
	      double all = PHYSICS->Isol->eflow->sumIsolation(mu,event.rec(),isocone,0.0,IsolationEFlow::ALL_COMPONENTS);
	      double isolation = (all)/mu->momentum().Pt();
	      if(isolation < 0.2) { isomuon.push_back(mu);}
	    }
	}
	
      // Iso track info

      vector<const RecTrackFormat*>  electronIsoTracks, muonIsoTracks, isoTracks ;
      
      for(MAuint32 i=0; i<event.rec()->tracks().size(); i++)
	{  const RecTrackFormat *isotrack = &(event.rec()->tracks()[i]);
	  double pt = isotrack->momentum().Pt(), eta =  isotrack->momentum().Eta(), mT = isotrack->mt_met(event.rec()->MET().momentum());
	  if(!(fabs(eta) < 2.4)) continue;  
	    
	  int particleId    =  isotrack->pdgid() ; // here I identify if the track is an electron or a muon, since this gives different iso requirements    
	  bool iselectron   = ( particleId== 11 or particleId == -11) ;
	  bool ismuon       = ( particleId== 13 or particleId == -13) ;                                               
	  double IsoCone    = 0.3 ;
	  double ChargedSum = PHYSICS->Isol->eflow->sumIsolation(isotrack,event.rec() , IsoCone ,0.,IsolationEFlow::TRACK_COMPONENT) / pt;
	  if((iselectron)&&(pt > 5)&&(mT < 100)&&(ChargedSum < 0.2))           
	    { electronIsoTracks.push_back(isotrack); }
	  if((ismuon)&&(pt > 5)&&(mT < 100)&&(ChargedSum < 0.2))
	    { muonIsoTracks.push_back(isotrack); }   
	  if((!(ismuon or iselectron))&&(pt > 10)&&(mT < 100)&&(ChargedSum < 0.1))
	    {isoTracks.push_back(isotrack); }                                                  

	}
    
      // initializing MET variable
      
      MALorentzVector misspt = event.rec()->MET().momentum();
      double rec_met = misspt.Pt();

      vector<const RecJetFormat*> rec_jets, isr_jets;
      double dr = 100.00;
      MAuint32 nbjets=0;
      double rec_ht=0;

      for (MAuint32 i=0;i<event.rec()->jets().size();i++)
	{ dr=100.0;
	  const RecJetFormat *jet = &(event.rec()->jets()[i]);
	  if(isophoton.size()!=0)
	    {
	      dr = DeltaR(isophoton[0]->momentum().Eta(),isophoton[0]->momentum().Phi(),jet->momentum().Eta(),jet->momentum().Phi());
	    }
	   if(jet->momentum().Pt()>30 and fabs(jet->momentum().Eta())<2.4 and dr>0.3)
	    { rec_ht+=jet->momentum().Pt();
	      rec_jets.push_back(jet);
	      if(jet->btag())
		{ nbjets++;}
	    }
	  else if(jet->momentum().Pt()>30 and fabs(jet->momentum().Eta())<2.4 and dr<0.3)
	    { 
	      rec_ht+=isophoton[0]->momentum().Pt();
	    }

	}

      SORTER->sort(rec_jets,PTordering);
      SORTER->sort(isomuon,PTordering);
      SORTER->sort(isoelec,PTordering);
      rec_jets = Removal(rec_jets,isomuon,0.2);
      rec_jets = Removal(rec_jets,isoelec,0.2);
      rec_jets = Removal(rec_jets,isophoton,0.3);
      MAuint32 rec_njets = rec_jets.size();

      double dphi1 =100.0,dphi2 =100.0;

      if(rec_jets.size()!=0)
	{
	  dphi1 = fabs(DeltaPhi(misspt.Phi(),rec_jets[0]->momentum().Phi()));
	}
      if(rec_jets.size()>1)
	{
	  dphi2 = fabs(DeltaPhi(misspt.Phi(),rec_jets[1]->momentum().Phi()));
	}

      if(isophoton.size()!=0)
	{double gammapt = isophoton[0]->momentum().Pt();

	  // ************************************************
	  // Applying pre selection cuts 
	  // ************************************************

	  // 1) gamma pt > 100
	
	  if(!Manager()->ApplyCut((gammapt > 100),"ptgamma>100"))
	     return true;
	 	  
	  // 2) vetoelectrons and vetomuons
	  
	  if(!Manager()->ApplyCut((isoelec.size() == 0 and isomuon.size() == 0),"vetoelectron&muons")) return true;
	  
	  // 3) no isoelectrontracks,iso muontracks and isopiontracks

	  if(!Manager()->ApplyCut((electronIsoTracks.size() == 0 and muonIsoTracks.size()==0 and  isoTracks.size()==0),"IsoTrackVeto")) return true;
	  
	  // 4) rec_met >200

	  if(!Manager()->ApplyCut((rec_met > 200),"rec_met>200")) return true;

	  // 5) rec_njets >= 2

	  if(!Manager()->ApplyCut((rec_njets >= 2),"njets>=2")) return true;
	  
	  // 6) (photonpt > 190 and rec_ht>500) or (photonpt>100 and rec_ht>800)

	  if(!Manager()->ApplyCut(((gammapt > 190 && rec_ht > 500) || (gammapt > 100 && rec_ht > 800)),"STcut")) return true;

	  // 7) |dphi(met,rec_jet[0])|>0.3    and |dphi(met,rec_jet[1])|>0.3

	  if(!Manager()->ApplyCut((dphi1 > 0.3 && dphi2 > 0.3),"dphicut")) return true;

	  // 8) btag>=1                    
	  
	  if(!Manager()->ApplyCut((nbjets >= 1),"bjetcut>1")) return true; 
	  
	  // 9) or btag == 0

	  if(!Manager()->ApplyCut((nbjets == 0),"bjetcut=0")) return true;
	  
	  // *********************************
	  // Applying the signal regions cuts
	  // *********************************
	  
	  if(!Manager()->ApplyCut((rec_njets>=2 and rec_njets<=4 and nbjets==0 and rec_met>=200 and rec_met<270),"SR1_NJ2-4_NB0_MISSPT270")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=2 and rec_njets<=4 and nbjets==0 and rec_met>=270 and rec_met<350),"SR2_NJ2-4_NB0_MISSPT350")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=2 and rec_njets<=4 and nbjets==0 and rec_met>=350 and rec_met<450),"SR3_NJ2-4_NB0_MISSPT450")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=2 and rec_njets<=4 and nbjets==0 and rec_met>=450 and rec_met<750),"SR4_NJ2-4_NB0_MISSPT750")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=2 and rec_njets<=4 and nbjets==0 and rec_met>=750),"SR5_NJ2-4_NB0_MISSPT>750")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=5 and rec_njets<=6 and nbjets==0 and rec_met>=200 and rec_met<270),"SR6_NJ5-6_NB0_MISSPT270")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=5 and rec_njets<=6 and nbjets==0 and rec_met>=270 and rec_met<350),"SR7_NJ5-6_NB0_MISSPT350")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=5 and rec_njets<=6 and nbjets==0 and rec_met>=350 and rec_met<450),"SR8_NJ5-6_NB0_MISSPT450")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=5 and rec_njets<=6 and nbjets==0 and rec_met>=450),"SR9_NJ5-6_NB0_MISSPT>450")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=7 and nbjets==0 and rec_met>=200 and rec_met<270),"SR10_NJ>7_NB0_MISSPT270")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=7 and nbjets==0 and rec_met>=270 and rec_met<350),"SR11_NJ>7_NB0_MISSPT350")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=7 and nbjets==0 and rec_met>=200 and rec_met<270),"SR12_NJ>7_NB0_MISSPT450")) return true;
	  
	  if(!Manager()->ApplyCut((rec_njets>=7 and nbjets==0 and rec_met>=450),"SR13_NJ>7_NB0_MISSPT>450")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=2 and rec_njets<=4 and nbjets>=1 and rec_met>=200 and rec_met<270),"SR14_NJ2-4_NB1_MISSPT270")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=2 and rec_njets<=4 and nbjets>=1 and rec_met>=270 and rec_met<350),"SR15_NJ2-4_NB1_MISSPT350")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=2 and rec_njets<=4 and nbjets>=1 and rec_met>=350 and rec_met<450),"SR16_NJ2-4_NB1_MISSPT450")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=2 and rec_njets<=4 and nbjets>=1 and rec_met>=450),"SR17_NJ2-4_NB1_MISSPT>450")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=5 and rec_njets<=6 and nbjets>=1 and rec_met>=200 and rec_met<270),"SR18_NJ5-6_NB1_MISSPT270")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=5 and rec_njets<=6 and nbjets>=1 and rec_met>=270 and rec_met<350),"SR19_NJ5-6_NB1_MISSPT350")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=5 and rec_njets<=6 and nbjets>=1 and rec_met>=350 and rec_met<450),"SR20_NJ5-6_NB1_MISSPT450")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=5 and rec_njets<=6 and nbjets>=1 and rec_met>=450),"SR21_NJ5-6_NB1_MISSPT>450")) return true;
	  
	  if(!Manager()->ApplyCut((rec_njets>=7 and nbjets>=1 and rec_met>=200 and rec_met<270),"SR22_NJ>7_NB1_MISSPT270")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=7 and nbjets>=1 and rec_met>=270 and rec_met<350),"SR23_NJ>7_NB1_MISSPT350")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=7 and nbjets>=1 and rec_met>=350 and rec_met<450),"SR24_NJ>7_NB1_MISSPT450")) return true;

	  if(!Manager()->ApplyCut((rec_njets>=7 and nbjets>=1 and rec_met>=450),"SR25_NJ>7_NB1_MISSPT>450")) return true;

     
	}
     
    }

  
  return true;
}

