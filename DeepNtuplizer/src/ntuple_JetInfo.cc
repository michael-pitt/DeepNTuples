/*
 * ntuple_JetInfo.cc
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */




#include "../interface/ntuple_JetInfo.h"
#include "../interface/helpers.h"
//#include "../interface/leptonsInJet.h"
#include <vector>
#include <algorithm>
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

void ntuple_JetInfo::getInput(const edm::ParameterSet& iConfig){

    gluonReduction_=(iConfig.getParameter<double>("gluonReduction"));
    jetPtMin_=(iConfig.getParameter<double>("jetPtMin"));
    jetPtMax_=(iConfig.getParameter<double>("jetPtMax"));
    jetAbsEtaMin_=(iConfig.getParameter<double>("jetAbsEtaMin"));
    jetAbsEtaMax_=(iConfig.getParameter<double>("jetAbsEtaMax"));

    vector<string> disc_names = iConfig.getParameter<vector<string> >("bDiscriminators");
    for(auto& name : disc_names) {
        discriminators_[name] = 0.;
    }
}
void ntuple_JetInfo::initBranches(TTree* tree){

    //more general event info, here applied per jet
    addBranch(tree,"npv"    ,&npv_    ,"npv/f"    );
    addBranch(tree,"rho", &rho_, "rho/f");
    addBranch(tree,"ntrueInt",&ntrueInt_,"ntrueInt/f");
    addBranch(tree,"event_no"    ,&event_no_    ,"event_no/i"    );
    addBranch(tree,"jet_no"    ,&jet_no_    ,"jet_no/i"    );


    // truth labels
    addBranch(tree,"gen_pt"    ,&gen_pt_    ,"gen_pt_/f"    );
    addBranch(tree,"gen_parton_pdgid"    ,&gen_parton_pdgid_    ,"gen_parton_pdgid_/I"    );
    addBranch(tree,"gen_hadron_pdgid"    ,&gen_hadron_pdgid_    ,"gen_hadron_pdgid_/I"    );
    addBranch(tree,"gen_hadron_pt"    ,&gen_hadron_pt_    ,"gen_hadron_pt_/f"    );
    addBranch(tree,"Delta_gen_pt"    ,&Delta_gen_pt_,"Delta_gen_pt_/f"    );
    addBranch(tree,"isB",&isB_, "isB_/i");
    addBranch(tree,"isGBB",&isGBB_, "isGBB_/i");
    addBranch(tree,"isBB",&isBB_, "isBB_/i");
    addBranch(tree,"isLeptonicB",&isLeptonicB_, "isLeptonicB_/i");
    addBranch(tree,"isLeptonicB_C",&isLeptonicB_C_, "isLeptonicB_C_/i");
    addBranch(tree,"isC",&isC_, "isC_/i");
    addBranch(tree,"isGCC",&isGCC_, "isGCC_/i");
    addBranch(tree,"isCC",&isCC_, "isCC_/i");
//    addBranch(tree,"isTau",&isTau_, "isTau_/i");
    addBranch(tree,"isUD",&isUD_, "isUD_/i");
    addBranch(tree,"isS",&isS_, "isS_/i");
    addBranch(tree,"isG",&isG_, "isG_/i");
    addBranch(tree,"isUndefined",&isUndefined_, "isUndefined_/i");
    addBranch(tree,"genDecay",&genDecay_, "genDecay_/f"); //dxy corresponds to the distance the Bhadron traveled

    //truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    addBranch(tree,"isPhysB",&isPhysB_, "isPhysB_/i");
    addBranch(tree,"isPhysGBB",&isPhysGBB_, "isPhysGBB_/i");
    addBranch(tree,"isPhysBB",&isPhysBB_, "isPhysBB_/i");
    addBranch(tree,"isPhysLeptonicB",&isPhysLeptonicB_, "isPhysLeptonicB_/i");
    addBranch(tree,"isPhysLeptonicB_C",&isPhysLeptonicB_C_, "isPhysLeptonicB_C_/i");
    addBranch(tree,"isPhysC",&isPhysC_, "isPhysC_/i");
    addBranch(tree,"isPhysGCC",&isPhysGCC_, "isPhysGCC_/i");
    addBranch(tree,"isPhysCC",&isPhysCC_, "isPhysCC_/i");
//    addBranch(tree,"isPhysTau",&isPhysTau_, "isPhysTau_/i");
    addBranch(tree,"isPhysUD",&isPhysUD_, "isPhysUD_/i");
    addBranch(tree,"isPhysS",&isPhysS_, "isPhysS_/i");
    addBranch(tree,"isPhysG",&isPhysG_, "isPhysG_/i");
    addBranch(tree,"isPhysUndefined",&isPhysUndefined_, "isPhysUndefined_/i");

    // jet variables
    //b=tree->Branch("jet_pt", &jet_pt_);
    addBranch(tree,"jet_pt", &jet_pt_);
    addBranch(tree,"jet_corr_pt", &jet_corr_pt_);
    addBranch(tree,"jet_eta", &jet_eta_);
    addBranch(tree,"jet_phi", &jet_phi_);
    addBranch(tree,"jet_mass", &jet_mass_);
    addBranch(tree,"jet_energy", &jet_energy_);


    //jet id
    addBranch(tree,"jet_looseId", &jet_looseId_);

    // quark gluon
    addBranch(tree,"jet_qgl",   &jet_qgl_);  // qg tagger from jmar
    addBranch(tree,"QG_ptD",   &QG_ptD_);   // momentum fraction per jet constituent
    addBranch(tree,"QG_axis2", &QG_axis2_); // jet shape i.e. gluon are wider than quarks
    addBranch(tree,"QG_mult",  &QG_mult_);  // multiplicity i.e. total num of PFcands reconstructed

    // yutas quark-gluon info
    addBranch(tree,"y_multiplicity"    ,&y_multiplicity_,"y_multiplicity_/f"    );
    addBranch(tree,"y_charged_multiplicity"    ,&y_charged_multiplicity_,"y_charged_multiplicity_/f"    );
    addBranch(tree,"y_neutral_multiplicity"    ,&y_neutral_multiplicity_,"y_neutral_multiplicity_/f"    );
    addBranch(tree,"y_ptD"    ,&y_ptD_,"y_ptD_/f"    );
    addBranch(tree,"y_axis1"    ,&y_axis1_,"y_axis1_/f"    );
    addBranch(tree,"y_axis2"    ,&y_axis2_,"y_axis2_/f"    );
    addBranch(tree,"y_pt_dr_log"    ,&y_pt_dr_log_,"y_pt_dr_log_/f"    );


    // in the jet

    addBranch(tree,"muons_number", &muons_number_, "muons_number_/f");
    addBranch(tree,"electrons_number", &electrons_number_, "electrons_number_/f");
    addBranch(tree,"mn", &mn_, "mn_/i");
    addBranch(tree,"en", &en_, "en_/i");

    addBranch(tree,"muons_isLooseMuon", &muons_isLooseMuon_, "muons_isLooseMuon_[mn_]/f");
    addBranch(tree,"muons_isTightMuon", &muons_isTightMuon_, "muons_isTightMuon_[mn_]/f");
    addBranch(tree,"muons_isSoftMuon", &muons_isSoftMuon_, "muons_isSoftMuon_[mn_]/f");
    addBranch(tree,"muons_isHighPtMuon", &muons_isHighPtMuon_, "muons_isHighPtMuon_[mn_]/f");
    addBranch(tree,"muons_pt", &muons_pt_, "muons_pt_[mn_]/f");
    addBranch(tree,"muons_relEta", &muons_relEta_, "muons_relEta_[mn_]/f");
    addBranch(tree,"muons_relPhi", &muons_relPhi_, "muons_relPhi_[mn_]/f");
    addBranch(tree,"muons_energy", &muons_energy_, "muons_energy_[mn_]/f");
    addBranch(tree,"muons_charge", &muons_charge_, "muons_charge_[mn_]/f");
    addBranch(tree,"electrons_pt", &electrons_pt_, "electrons_pt_[en_]/f");
    addBranch(tree,"electrons_relEta", &electrons_relEta_, "electrons_relEta_[en_]/f");
    addBranch(tree,"electrons_relPhi", &electrons_relPhi_, "electrons_relPhi_[en_]/f");
    addBranch(tree,"electrons_energy", &electrons_energy_, "electrons_energy_[en_]/f");
    addBranch(tree,"electrons_charge", &electrons_charge_, "electrons_charge_[en_]/f");

	addBranch(tree,"muons_isGlobal", &muons_isGlobal_,"muons_isGlobal_[mn_]/f");
	addBranch(tree,"muons_isStandAlone", &muons_isStandAlone_,"muons_isStandAlone_[mn_]/f");
	addBranch(tree,"muons_jetDeltaR", &muons_jetDeltaR_,"muons_jetDeltaR_[mn_]/f");
	addBranch(tree,"muons_numberOfMatchedStations", &muons_numberOfMatchedStations_,"muons_numberOfMatchedStations_[mn_]/f");
	addBranch(tree,"muons_2dIP", &muons_2dIP_,"muons_2dIP_[mn_]/f");
	addBranch(tree,"muons_2dIPSig", &muons_2dIPSig_,"muons_2dIPSig_[mn_]/f");
	addBranch(tree,"muons_3dIP", &muons_3dIP_,"muons_3dIP_[mn_]/f");
	addBranch(tree,"muons_3dIPSig", &muons_3dIPSig_,"muons_3dIPSig_[mn_]/f");
	addBranch(tree,"muons_dxy", &muons_dxy_,"muons_dxy_[mn_]/f");
	addBranch(tree,"muons_dxyError", &muons_dxyError_,"muons_dxyError_[mn_]/f");
	addBranch(tree,"muons_dxySig", &muons_dxySig_,"muons_dxySig_[mn_]/f");
	addBranch(tree,"muons_dz", &muons_dz_,"muons_dz_[mn_]/f");
	addBranch(tree,"muons_dzError", &muons_dzError_,"muons_dzError_[mn_]/f");
	addBranch(tree,"muons_numberOfValidPixelHits", &muons_numberOfValidPixelHits_,"muons_numberOfValidPixelHits_[mn_]/f");
	addBranch(tree,"muons_numberOfpixelLayersWithMeasurement", &muons_numberOfpixelLayersWithMeasurement_,"muons_numberOfpixelLayersWithMeasurement_[mn_]/f");
	addBranch(tree,"muons_numberOfstripLayersWithMeasurement", &muons_numberOfstripLayersWithMeasurement_,"muons_numberOfstripLayersWithMeasurement_[mn_]/f");
	addBranch(tree,"muons_chi2", &muons_chi2_,"muons_chi2_[mn_]/f");
	addBranch(tree,"muons_ndof", &muons_ndof_,"muons_ndof_[mn_]/f");
	addBranch(tree,"muons_caloIso", &muons_caloIso_,"muons_caloIso_[mn_]/f");
	addBranch(tree,"muons_ecalIso", &muons_ecalIso_,"muons_ecalIso_[mn_]/f");
	addBranch(tree,"muons_hcalIso", &muons_hcalIso_,"muons_hcalIso_[mn_]/f");
	addBranch(tree,"muons_sumPfChHadronPt", &muons_sumPfChHadronPt_,"muons_sumPfChHadronPt_[mn_]/f");
	addBranch(tree,"muons_sumPfNeuHadronEt", &muons_sumPfNeuHadronEt_,"muons_sumPfNeuHadronEt_[mn_]/f");
	addBranch(tree,"muons_Pfpileup", &muons_Pfpileup_,"muons_Pfpileup_[mn_]/f");
	addBranch(tree,"muons_sumPfPhotonEt", &muons_sumPfPhotonEt_,"muons_sumPfPhotonEt_[mn_]/f");
	addBranch(tree,"muons_sumPfChHadronPt03", &muons_sumPfChHadronPt03_,"muons_sumPfChHadronPt03_[mn_]/f");
	addBranch(tree,"muons_sumPfNeuHadronEt03", &muons_sumPfNeuHadronEt03_,"muons_sumPfNeuHadronEt03_[mn_]/f");
	addBranch(tree,"muons_Pfpileup03", &muons_Pfpileup03_,"muons_Pfpileup03_[mn_]/f");
	addBranch(tree,"muons_sumPfPhotonEt03", &muons_sumPfPhotonEt03_,"muons_sumPfPhotonEt03_[mn_]/f");
	addBranch(tree,"electrons_jetDeltaR", &electrons_jetDeltaR_,"electrons_jetDeltaR_[en_]/f");
	addBranch(tree,"electrons_EtFromCaloEn", &electrons_EtFromCaloEn_,"electrons_EtFromCaloEn_[en_]/f");
	addBranch(tree,"electrons_ecalDrivenSeed", &electrons_ecalDrivenSeed_,"electrons_ecalDrivenSeed_[en_]/f");
	addBranch(tree,"electrons_isEB", &electrons_isEB_,"electrons_isEB_[en_]/f");
	addBranch(tree,"electrons_isEE", &electrons_isEE_,"electrons_isEE_[en_]/f");
	addBranch(tree,"electrons_ecalEnergy", &electrons_ecalEnergy_,"electrons_ecalEnergy_[en_]/f");
	addBranch(tree,"electrons_isPassConversionVeto", &electrons_isPassConversionVeto_,"electrons_isPassConversionVeto_[en_]/f");
	addBranch(tree,"electrons_convDist", &electrons_convDist_,"electrons_convDist_[en_]/f");
	addBranch(tree,"electrons_convFlags", &electrons_convFlags_,"electrons_convFlags_[en_]/f");
	addBranch(tree,"electrons_convRadius", &electrons_convRadius_,"electrons_convRadius_[en_]/f");
	addBranch(tree,"electrons_3dIP", &electrons_3dIP_,"electrons_3dIP_[en_]/f");
	addBranch(tree,"electrons_3dIPSig", &electrons_3dIPSig_,"electrons_3dIPSig_[en_]/f");
	addBranch(tree,"electrons_2dIP", &electrons_2dIP_,"electrons_2dIP_[en_]/f");
	addBranch(tree,"electrons_2dIPSig", &electrons_2dIPSig_,"electrons_2dIPSig_[en_]/f");
	addBranch(tree,"electrons_sCseedEta", &electrons_sCseedEta_,"electrons_sCseedEta_[en_]/f");
	addBranch(tree,"electrons_eSeedClusterOverP", &electrons_eSeedClusterOverP_,"electrons_eSeedClusterOverP_[en_]/f");
	addBranch(tree,"electrons_eSeedClusterOverPout", &electrons_eSeedClusterOverPout_,"electrons_eSeedClusterOverPout_[en_]/f");
	addBranch(tree,"electrons_eSuperClusterOverP", &electrons_eSuperClusterOverP_,"electrons_eSuperClusterOverP_[en_]/f");
	addBranch(tree,"electrons_hadronicOverEm", &electrons_hadronicOverEm_,"electrons_hadronicOverEm_[en_]/f");
	addBranch(tree,"electrons_deltaEtaEleClusterTrackAtCalo", &electrons_deltaEtaEleClusterTrackAtCalo_,"electrons_deltaEtaEleClusterTrackAtCalo_[en_]/f");
	addBranch(tree,"electrons_deltaPhiEleClusterTrackAtCalo", &electrons_deltaPhiEleClusterTrackAtCalo_,"electrons_deltaPhiEleClusterTrackAtCalo_[en_]/f");
	addBranch(tree,"electrons_deltaEtaSeedClusterTrackAtCalo", &electrons_deltaEtaSeedClusterTrackAtCalo_,"electrons_deltaEtaSeedClusterTrackAtCalo_[en_]/f")    ;
	addBranch(tree,"electrons_deltaPhiSeedClusterTrackAtCalo", &electrons_deltaPhiSeedClusterTrackAtCalo_,"electrons_deltaPhiSeedClusterTrackAtCalo_[en_]/f")    ;
	addBranch(tree,"electrons_deltaEtaSeedClusterTrackAtVtx", &electrons_deltaEtaSeedClusterTrackAtVtx_,"electrons_deltaEtaSeedClusterTrackAtVtx_[en_]/f");
	addBranch(tree,"electrons_deltaEtaSuperClusterTrackAtVtx", &electrons_deltaEtaSuperClusterTrackAtVtx_,"electrons_deltaEtaSuperClusterTrackAtVtx_[en_]/f")    ;
	addBranch(tree,"electrons_deltaPhiSuperClusterTrackAtVtx", &electrons_deltaPhiSuperClusterTrackAtVtx_,"electrons_deltaPhiSuperClusterTrackAtVtx_[en_]/f")    ;
	addBranch(tree,"electrons_dxy", &electrons_dxy_,"electrons_dxy_[en_]/f");
	addBranch(tree,"electrons_dz", &electrons_dz_,"electrons_dz_[en_]/f");
	addBranch(tree,"electrons_nbOfMissingHits", &electrons_nbOfMissingHits_,"electrons_nbOfMissingHits_[en_]/f");
	addBranch(tree,"electrons_gsfCharge", &electrons_gsfCharge_,"electrons_gsfCharge_[en_]/f");
	addBranch(tree,"electrons_ndof", &electrons_ndof_,"electrons_ndof_[en_]/f");
	addBranch(tree,"electrons_chi2", &electrons_chi2_,"electrons_chi2_[en_]/f");
	addBranch(tree,"electrons_SC", &electrons_SC_energy_,"electrons_SC_energy_[en_]/f");
	addBranch(tree,"electrons_SC", &electrons_SC_deta_,"electrons_SC_deta_[en_]/f");
	addBranch(tree,"electrons_SC", &electrons_SC_dphi_,"electrons_SC_dphi_[en_]/f");
	addBranch(tree,"electrons_SC", &electrons_SC_et_,"electrons_SC_et_[en_]/f");
	addBranch(tree,"electrons_scPixCharge", &electrons_scPixCharge_,"electrons_scPixCharge_[en_]/f");
	addBranch(tree,"electrons_numberOfBrems", &electrons_numberOfBrems_,"electrons_numberOfBrems_[en_]/f");
	addBranch(tree,"electrons_fbrem", &electrons_fbrem_,"electrons_fbrem_[en_]/f");
	addBranch(tree,"electrons_sigmaEtaEta", &electrons_sigmaEtaEta_,"electrons_sigmaEtaEta_[en_]/f");
	addBranch(tree,"electrons_sigmaIetaIeta", &electrons_sigmaIetaIeta_,"electrons_sigmaIetaIeta_[en_]/f");
	addBranch(tree,"electrons_sigmaIphiIphi", &electrons_sigmaIphiIphi_,"electrons_sigmaIphiIphi_[en_]/f");
	addBranch(tree,"electrons_r9", &electrons_r9_,"electrons_r9_[en_]/f");
	addBranch(tree,"electrons_superClusterFbrem", &electrons_superClusterFbrem_,"electrons_superClusterFbrem_[en_]/f");
	addBranch(tree,"electrons_e5x5", &electrons_e5x5_,"electrons_e5x5_[en_]/f");
	addBranch(tree,"electrons_e5x5Rel", &electrons_e5x5Rel_,"electrons_e5x5Rel_[en_]/f");
	addBranch(tree,"electrons_e1x5Overe5x5", &electrons_e1x5Overe5x5_,"electrons_e1x5Overe5x5_[en_]/f");
	addBranch(tree,"electrons_e2x5MaxOvere5x5", &electrons_e2x5MaxOvere5x5_,"electrons_e2x5MaxOvere5x5_[en_]/f");
	addBranch(tree,"electrons_hcalOverEcal", &electrons_hcalOverEcal_,"electrons_hcalOverEcal_[en_]/f");
	addBranch(tree,"electrons_SC", &electrons_SC_eSuperClusterOverP_,"electrons_SC_eSuperClusterOverP_[en_]/f");
	addBranch(tree,"electrons_neutralHadronIso", &electrons_neutralHadronIso_,"electrons_neutralHadronIso_[en_]/f");
	addBranch(tree,"electrons_photonIso", &electrons_photonIso_,"electrons_photonIso_[en_]/f");
	addBranch(tree,"electrons_puChargedHadronIso", &electrons_puChargedHadronIso_,"electrons_puChargedHadronIso_[en_]/f");
	addBranch(tree,"electrons_trackIso", &electrons_trackIso_,"electrons_trackIso_[en_]/f");
	addBranch(tree,"electrons_hcalDepth1OverEcal", &electrons_hcalDepth1OverEcal_,"electrons_hcalDepth1OverEcal_[en_]/f");
	addBranch(tree,"electrons_hcalDepth2OverEcal", &electrons_hcalDepth2OverEcal_,"electrons_hcalDepth2OverEcal_[en_]/f");
	addBranch(tree,"electrons_ecalPFClusterIso", &electrons_ecalPFClusterIso_,"electrons_ecalPFClusterIso_[en_]/f");
	addBranch(tree,"electrons_hcalPFClusterIso", &electrons_hcalPFClusterIso_,"electrons_hcalPFClusterIso_[en_]/f");
	addBranch(tree,"electrons_pfSumPhotonEt", &electrons_pfSumPhotonEt_,"electrons_pfSumPhotonEt_[en_]/f");
	addBranch(tree,"electrons_pfSumChargedHadronPt", &electrons_pfSumChargedHadronPt_,"electrons_pfSumChargedHadronPt_[en_]/f");
	addBranch(tree,"electrons_pfSumNeutralHadronEt", &electrons_pfSumNeutralHadronEt_,"electrons_pfSumNeutralHadronEt_[en_]/f");
	addBranch(tree,"electrons_pfSumPUPt", &electrons_pfSumPUPt_,"electrons_pfSumPUPt_[en_]/f");
	addBranch(tree,"electrons_dr04TkSumPt", &electrons_dr04TkSumPt_,"electrons_dr04TkSumPt_[en_]/f");
	addBranch(tree,"electrons_dr04EcalRecHitSumEt", &electrons_dr04EcalRecHitSumEt_,"electrons_dr04EcalRecHitSumEt_[en_]/f");
	addBranch(tree,"electrons_dr04HcalDepth1TowerSumEt", &electrons_dr04HcalDepth1TowerSumEt_,"electrons_dr04HcalDepth1TowerSumEt_[en_]/f");
	addBranch(tree,"electrons_dr04HcalDepth1TowerSumEtBc", &electrons_dr04HcalDepth1TowerSumEtBc_,"electrons_dr04HcalDepth1TowerSumEtBc_[en_]/f");
	addBranch(tree,"electrons_dr04HcalDepth2TowerSumEt", &electrons_dr04HcalDepth2TowerSumEt_,"electrons_dr04HcalDepth2TowerSumEt_[en_]/f");
	addBranch(tree,"electrons_dr04HcalDepth2TowerSumEtBc", &electrons_dr04HcalDepth2TowerSumEtBc_,"electrons_dr04HcalDepth2TowerSumEtBc_[en_]/f");
	addBranch(tree,"electrons_dr04HcalTowerSumEt", &electrons_dr04HcalTowerSumEt_,"electrons_dr04HcalTowerSumEt_[en_]/f");
	addBranch(tree,"electrons_dr04HcalTowerSumEtBc", &electrons_dr04HcalTowerSumEtBc_,"electrons_dr04HcalTowerSumEtBc_[en_]/f");

    addBranch(tree,"gen_pt_Recluster"    ,&gen_pt_Recluster_    ,"gen_pt_Recluster_/f"    );
    addBranch(tree,"gen_pt_WithNu"    ,&gen_pt_WithNu_    ,"gen_pt_WithNu_/f"    );
    addBranch(tree,"Delta_gen_pt_Recluster"    ,&Delta_gen_pt_Recluster_    ,"Delta_gen_pt_Recluster_/f"    );
    addBranch(tree,"Delta_gen_pt_WithNu"    ,&Delta_gen_pt_WithNu_    ,"Delta_gen_pt_WithNu_/f"    );

    if(1) // discriminators might need to be filled differently. FIXME
        for(auto& entry : discriminators_) {
            string better_name(entry.first);
            std::replace(better_name.begin(), better_name.end(), ':', '_');
            addBranch(tree,better_name.c_str(), &entry.second, (better_name+"/F").c_str());
        }
}


void ntuple_JetInfo::readEvent(const edm::Event& iEvent){

    iEvent.getByToken(qglToken_, qglHandle);
    iEvent.getByToken(ptDToken_, ptDHandle);
    iEvent.getByToken(axis2Token_, axis2Handle);
    iEvent.getByToken(multToken_, multHandle);

    iEvent.getByToken(genJetMatchReclusterToken_, genJetMatchRecluster);
    iEvent.getByToken(genJetMatchWithNuToken_, genJetMatchWithNu);

    iEvent.getByToken(genParticlesToken_, genParticlesHandle);
    iEvent.getByToken(jetFlavourInfosToken_, jetFlavourInfos);


    iEvent.getByToken(muonsToken_, muonsHandle);
    iEvent.getByToken(electronsToken_, electronsHandle);

    event_no_=iEvent.id().event();

    //presumably this whole part can be removed!


    neutrinosLepB.clear();
    neutrinosLepB_C.clear();
    gToBB.clear();
    gToCC.clear();
    alltaus_.clear();
    Bhadron_.clear();
    Bhadron_daughter_.clear();

    //std::cout << " start search for a b in this event "<<std::endl;
 for (const reco::Candidate &genC : *genParticlesHandle)
   {
     const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);
     
     if((abs(gen.pdgId())>500&&abs(gen.pdgId())<600)||(abs(gen.pdgId())>5000&&abs(gen.pdgId())<6000)) {

       //std::cout<<gen.end_vertex()<<endl;

       Bhadron_.push_back(gen);
       if(gen.numberOfDaughters()>0){
     
	 if( (abs(gen.daughter(0)->pdgId())>500&&abs(gen.daughter(0)->pdgId())<600)||(abs(gen.daughter(0)->pdgId())>5000&&abs(gen.daughter(0)->pdgId())<6000))
	   {
	     if(gen.daughter(0)->numberOfDaughters()>0)
	       {
		
		 const reco::GenParticle &daughter_ = static_cast< const reco::GenParticle &>(*(gen.daughter(0)->daughter(0)));
		 
		 if(daughter_.vx()!=gen.vx())
		   { 
		     Bhadron_daughter_.push_back(daughter_);
		   }
                 else Bhadron_daughter_.push_back(gen);
		 //  std::cout << "only b daughters " << endl;
		 // }
	       }
	     else  Bhadron_daughter_.push_back(gen);
	     
	   }
	 else{
	   //  std::cout<<gen.daughter(0)->vx()<< " oh  " <<gen.vx()<<" "<<gen.pt() <<" "<<  gen.daughter(0)->pdgId() <<std::endl; 
	  
	   const reco::GenParticle &daughter_ = static_cast< const reco::GenParticle &>(*gen.daughter(0));
	   Bhadron_daughter_.push_back(daughter_);
	 }

       }// if daughter is there
       else {
	 
	 //std::cout << " lonly B hadron, has NO daughter??? "<<std::endl;
	 Bhadron_daughter_.push_back(gen);
       }
     }
   }

 for (const reco::Candidate &genC : *genParticlesHandle) {
        const reco::GenParticle &gen = static_cast< const reco::GenParticle &>(genC);
        if(abs(gen.pdgId())==12||abs(gen.pdgId())==14||abs(gen.pdgId())==16) {
            const reco::GenParticle* mother =  static_cast< const reco::GenParticle*> (gen.mother());
            if(mother!=NULL) {
                if((abs(mother->pdgId())>500&&abs(mother->pdgId())<600)||(abs(mother->pdgId())>5000&&abs(mother->pdgId())<6000)) {
                    neutrinosLepB.emplace_back(gen);
                }
                if((abs(mother->pdgId())>400&&abs(mother->pdgId())<500)||(abs(mother->pdgId())>4000&&abs(mother->pdgId())<5000)) {
                    neutrinosLepB_C.emplace_back(gen);
                }
            }
            else {
                std::cout << "No mother" << std::endl;
            }
        }

        int id(std::abs(gen.pdgId())); 
        int status(gen.status());

        if (id == 21 && status >= 21 && status <= 59) { //// Pythia8 hard scatter, ISR, or FSR
            if ( gen.numberOfDaughters() == 2 ) {
                const reco::Candidate* d0 = gen.daughter(0);
                const reco::Candidate* d1 = gen.daughter(1);
                if ( std::abs(d0->pdgId()) == 5 && std::abs(d1->pdgId()) == 5
                        && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4) gToBB.push_back(gen) ;
                if ( std::abs(d0->pdgId()) == 4 && std::abs(d1->pdgId()) == 4
                        && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR(*d0, *d1) < 0.4) gToCC.push_back(gen) ;
            }
        }

        if(id == 15 && false){
            alltaus_.push_back(gen);
        }

    }
    //technically a branch fill but per event, therefore here
}

//use either of these functions

bool ntuple_JetInfo::fillBranches(const pat::Jet & jet, const size_t& jetidx, const edm::View<pat::Jet> * coll){
    if(!coll)
        throw std::runtime_error("ntuple_JetInfo::fillBranches: no jet collection");

    /// cuts ///
    bool returnval=true;

    // some cuts to contrin training region
    if ( jet.pt() < jetPtMin_ ||  jet.pt() > jetPtMax_ ) returnval=false;                  // apply jet pT cut
    if ( fabs(jet.eta()) < jetAbsEtaMin_ || fabs(jet.eta()) > jetAbsEtaMax_ ) returnval=false; // apply jet eta cut


    // often we have way to many gluons that we do not need. This randomply reduces the gluons
    if (gluonReduction_>0 && jet.partonFlavour()==21)
        if(TRandom_.Uniform()>gluonReduction_) returnval=false;

    if(jet.genJet()==NULL)returnval=false;


    //branch fills
    for(auto& entry : discriminators_) {
        entry.second = catchInfs(jet.bDiscriminator(entry.first),-0.1);
    }

    npv_ = vertices()->size();

    for (auto const& v : *pupInfo()) {
        int bx = v.getBunchCrossing();
        if (bx == 0) {
            ntrueInt_ = v.getTrueNumInteractions();
        }
    }
    rho_ = rhoInfo()[0];


    jet_no_=jetidx;

    const auto jetRef = reco::CandidatePtr(coll->ptrs().at( jetidx));

    jet_qgl_ = (*qglHandle)[jetRef];
    QG_ptD_ = (*ptDHandle)[jetRef];
    QG_axis2_ = (*axis2Handle)[jetRef];
    QG_mult_ = (*multHandle)[jetRef];


    //std::vector<Ptr<pat::Jet> > p= coll->ptrs();

    isB_=0; isGBB_=0; isBB_=0; isC_=0; isGCC_=0; isCC_=0; isUD_=0;isTau_=0;
    isS_=0; isG_=0, isLeptonicB_=0, isLeptonicB_C_=0, isUndefined_=0;
    auto muIds = deep_ntuples::jet_muonsIds(jet,*muonsHandle);
    auto elecIds = deep_ntuples::jet_electronsIds(jet,*electronsHandle);

    mn_ = 0;
    en_ = 0;

    float etasign = 1.;
    if (jet.eta()<0) etasign = -1.;

    for(std::size_t i=0; i<max_num_lept; i++) {
        if (i < muIds.size()) {
            const auto & muon = (*muonsHandle).at(muIds.at(i));
            muons_isLooseMuon_[i] = float(muon.isLooseMuon());
            muons_isTightMuon_[i] = float(muon.isTightMuon(vertices()->at(0)));
            muons_isSoftMuon_[i] = float(muon.isSoftMuon(vertices()->at(0)));
            muons_isHighPtMuon_[i] = float(muon.isHighPtMuon(vertices()->at(0)));
            muons_pt_[i] = muon.pt();
            muons_relEta_[i] = etasign*(muon.eta()-jet.eta());
            muons_relPhi_[i] = reco::deltaPhi(muon.phi(),jet.phi());
            muons_energy_[i] = muon.energy()/jet.energy();
			muons_charge_[i] = muon.charge();
            // Extra vars
            muons_isGlobal_[i] = float(muon.isLooseMuon());
            muons_isStandAlone_[i] = float(muon.isStandAloneMuon());
            muons_jetDeltaR_[i] = reco::deltaR(muon, jet);
            muons_numberOfMatchedStations_[i] = muon.numberOfMatchedStations();
            muons_2dIP_[i] = muon.dB();
            muons_2dIPSig_[i] = muon.dB()/muon.edB();
            muons_3dIP_[i] = muon.dB(pat::Muon::PV3D);
            muons_3dIPSig_[i] = muon.dB(pat::Muon::PV3D)/muon.edB(pat::Muon::PV3D);
            muons_dxy_[i] = muon.bestTrack()->dxy(vertices()->at(0).position());
            muons_dxyError_[i] = muon.bestTrack()->dxyError();
            muons_dxySig_[i] = muon.bestTrack()->dxy(vertices()->at(0).position())/muon.bestTrack()->dxyError(); 
            muons_dz_[i] = muon.bestTrack()->dz(vertices()->at(0).position());
            muons_dzError_[i] = muon.bestTrack()->dzError();
            muons_numberOfValidPixelHits_[i] = muon.bestTrack()->hitPattern().numberOfValidPixelHits();
            muons_numberOfpixelLayersWithMeasurement_[i] = muon.bestTrack()->hitPattern().pixelLayersWithMeasurement();
            muons_numberOfstripLayersWithMeasurement_[i] = muon.bestTrack()->hitPattern().stripLayersWithMeasurement();
            muons_chi2_[i] = muon.bestTrack()->chi2();
            muons_ndof_[i] = muon.bestTrack()->ndof();
            muons_caloIso_[i] =  muon.caloIso()/muon.pt();
            muons_ecalIso_[i] =  muon.ecalIso()/muon.pt(); 
            muons_hcalIso_[i] =  muon.hcalIso()/muon.pt();     
            muons_sumPfChHadronPt_[i]  = muon.pfIsolationR04().sumChargedHadronPt/muon.pt();
            muons_sumPfNeuHadronEt_[i]  = muon.pfIsolationR04().sumNeutralHadronEt/muon.pt();
            muons_Pfpileup_[i]  = muon.pfIsolationR04().sumPUPt/muon.pt();
            muons_sumPfPhotonEt_[i] = muon.pfIsolationR04().sumPhotonEt/muon.pt();
            muons_sumPfChHadronPt03_[i]  = muon.pfIsolationR03().sumChargedHadronPt/muon.pt();
            muons_sumPfNeuHadronEt03_[i]  = muon.pfIsolationR03().sumNeutralHadronEt/muon.pt();
            muons_Pfpileup03_[i]  = muon.pfIsolationR03().sumPUPt/muon.pt();
            muons_sumPfPhotonEt03_[i] = muon.pfIsolationR03().sumPhotonEt/muon.pt();       
			mn_++;
        }
        if (i < elecIds.size()) {
            const auto & electron = (*electronsHandle).at(elecIds.at(i));
            electrons_pt_[i] = electron.pt();
            electrons_relEta_[i] = etasign*(electron.eta()-jet.eta());
            electrons_relPhi_[i] = reco::deltaPhi(electron.phi(),jet.phi());
            electrons_energy_[i] = electron.energy()/jet.energy();
			electrons_charge_[i] = electron.charge();
            // Extra vars
			electrons_jetDeltaR_[i] = reco::deltaR(electron,jet); 
			electrons_EtFromCaloEn_[i] = electron.caloEnergy() * sin(electron.p4().theta())/ electron.pt();
			electrons_ecalDrivenSeed_[i] = electron.ecalDrivenSeed();

			electrons_isEB_[i] = electron.isEB();
			electrons_isEE_[i] = electron.isEE();
			electrons_ecalEnergy_[i] = electron.ecalEnergy()/electron.pt();
			electrons_isPassConversionVeto_[i] = electron.passConversionVeto();
			if(electron.convDist() >= 0.){
				electrons_convDist_[i] = electron.convDist(); 
				electrons_convFlags_[i] = electron.convFlags(); 
				electrons_convRadius_[i] = electron.convRadius();
			}
			else{
				electrons_convDist_[i] = -1.; 
				electrons_convFlags_[i] = -1.; 
				electrons_convRadius_[i] = -1.;
			}


			electrons_3dIP_[i] = electron.dB(pat::Electron::PV3D); 
			electrons_3dIPSig_[i] = electron.dB(pat::Electron::PV3D); 
			electrons_2dIP_[i] = electron.dB();
			electrons_2dIPSig_[i] = electron.dB()/electron.edB();

			if (std::isnan(electrons_2dIPSig_[i]) || std::isnan(electrons_3dIPSig_[i]))
			{
				electrons_2dIPSig_[i] = 0.;
				electrons_3dIPSig_[i] = 0.;
			}

			electrons_sCseedEta_[i] = electron.superCluster()->seed()->eta();


			electrons_eSeedClusterOverP_[i] = electron.eSeedClusterOverP();
			electrons_eSeedClusterOverPout_[i] = electron.eSeedClusterOverPout();
			electrons_eSuperClusterOverP_[i] = electron.eSuperClusterOverP();
			electrons_hadronicOverEm_[i] = electron.hadronicOverEm();


			electrons_deltaEtaEleClusterTrackAtCalo_[i] = electron.deltaEtaEleClusterTrackAtCalo();
			electrons_deltaPhiEleClusterTrackAtCalo_[i] = electron.deltaPhiEleClusterTrackAtCalo();

			electrons_deltaEtaSeedClusterTrackAtCalo_[i] = electron.deltaEtaSeedClusterTrackAtCalo(); 
			electrons_deltaPhiSeedClusterTrackAtCalo_[i] = electron.deltaPhiSeedClusterTrackAtCalo(); 

			electrons_deltaEtaSeedClusterTrackAtVtx_[i] = electron.deltaEtaSeedClusterTrackAtVtx(); 
			electrons_deltaEtaSuperClusterTrackAtVtx_[i] = electron.deltaEtaSuperClusterTrackAtVtx();  
			electrons_deltaPhiSuperClusterTrackAtVtx_[i] = electron.deltaPhiSuperClusterTrackAtVtx();

			electrons_dxy_[i] = electron.gsfTrack()->dxy(vertices()->at(0).position());
			electrons_dz_[i] = electron.gsfTrack()->dz(vertices()->at(0).position());
			electrons_nbOfMissingHits_[i] = electron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
			electrons_gsfCharge_[i] = electron.gsfTrack()->charge();
			electrons_ndof_[i] = electron.gsfTrack()->ndof(); 
			electrons_chi2_[i] = electron.gsfTrack()->chi2();
			electrons_SC_energy_[i] = electron.superCluster()->energy()/electron.pt(); 
			electrons_SC_deta_[i] = electron.superCluster()->eta()-electron.gsfTrack()->eta();
			electrons_SC_dphi_[i] = reco::deltaPhi(electron.superCluster()->phi(),electron.gsfTrack()->phi());
			electrons_SC_et_[i] = electron.superCluster()->energy() * sin(electron.p4().theta())/electron.pt();
			electrons_scPixCharge_[i] = electron.scPixCharge();
			electrons_numberOfBrems_[i] = electron.numberOfBrems();
			if(electron.pt() >= 5.){
				electrons_fbrem_[i] = electron.fbrem();
				electrons_sigmaEtaEta_[i] = electron.sigmaEtaEta();
				electrons_sigmaIetaIeta_[i] = electron.sigmaIetaIeta();
				electrons_sigmaIphiIphi_[i] = electron.sigmaIphiIphi();
				electrons_r9_[i] = electron.r9();
				electrons_superClusterFbrem_[i] = electron.superClusterFbrem();
			}
			else 
			{
				electrons_fbrem_[i] = -1.;
				electrons_sigmaEtaEta_[i] = -1.;
				electrons_sigmaIetaIeta_[i] = -1.;
				electrons_sigmaIphiIphi_[i] = -1;
				electrons_superClusterFbrem_[i] = -1.;
			}
			electrons_e5x5_[i] = electron.e5x5();
			electrons_e5x5Rel_[i] = electron.e5x5()/jet.pt();
			electrons_e1x5Overe5x5_[i] = electron.e1x5()/electron.e5x5();
			electrons_e2x5MaxOvere5x5_[i] = electron.e2x5Max()/electron.e5x5();
			if (electron.e5x5() == 0){
				electrons_e1x5Overe5x5_[i] = -1.;
				electrons_e2x5MaxOvere5x5_[i] = -1.;
			}
			electrons_hcalOverEcal_[i] = electron.hcalOverEcal();
			electrons_SC_eSuperClusterOverP_[i] = electron.eSuperClusterOverP();
			electrons_neutralHadronIso_[i] = electron.neutralHadronIso()/electron.pt();
			electrons_photonIso_[i] = electron.photonIso()/electron.pt(); 
			electrons_puChargedHadronIso_[i] = electron.puChargedHadronIso()/electron.pt(); 
			electrons_trackIso_[i] = electron.trackIso()/electron.pt();
			electrons_hcalDepth1OverEcal_[i] = electron.hcalDepth1OverEcal(); 
			electrons_hcalDepth2OverEcal_[i] = electron.hcalDepth2OverEcal();  
			electrons_ecalPFClusterIso_[i] = electron.ecalPFClusterIso()/electron.pt(); 
			electrons_hcalPFClusterIso_[i] = electron.hcalPFClusterIso()/electron.pt(); 
			electrons_pfSumPhotonEt_[i] = electron.pfIsolationVariables().sumPhotonEt/electron.pt(); 
			electrons_pfSumChargedHadronPt_[i] = electron.pfIsolationVariables().sumChargedHadronPt/electron.pt(); 
			electrons_pfSumNeutralHadronEt_[i] = electron.pfIsolationVariables().sumNeutralHadronEt/electron.pt(); 
			electrons_pfSumPUPt_[i] = electron.pfIsolationVariables().sumPUPt/electron.pt(); 
			electrons_dr04TkSumPt_[i] = electron.dr04TkSumPt()/electron.pt();
			electrons_dr04EcalRecHitSumEt_[i] = electron.dr04EcalRecHitSumEt()/electron.pt(); 
			electrons_dr04HcalDepth1TowerSumEt_[i] = electron.dr04HcalDepth1TowerSumEt()/electron.pt(); 
			electrons_dr04HcalDepth1TowerSumEtBc_[i] = electron.dr04HcalDepth1TowerSumEtBc()/electron.pt(); 
			electrons_dr04HcalDepth2TowerSumEt_[i] = electron.dr04HcalDepth2TowerSumEt()/electron.pt(); 
			electrons_dr04HcalDepth2TowerSumEtBc_[i] = electron.dr04HcalDepth2TowerSumEtBc()/electron.pt();
			electrons_dr04HcalTowerSumEt_[i] = electron.dr04HcalTowerSumEt()/electron.pt();
			electrons_dr04HcalTowerSumEtBc_[i] = electron.dr04HcalTowerSumEtBc()/electron.pt();
			en_++;
        }
    }
    muons_number_ = float(mn_);
    electrons_number_ = float(en_);

    //// Note that jets with gluon->bb (cc) and x->bb (cc) are in the same categories
    if(jet.genJet()!=NULL){
        switch(deep_ntuples::jet_flavour(jet, gToBB, gToCC, neutrinosLepB, neutrinosLepB_C, alltaus_)) {
        case deep_ntuples::JetFlavor::B:  isB_=1; break;
        case deep_ntuples::JetFlavor::LeptonicB: isLeptonicB_=1; break;
        case deep_ntuples::JetFlavor::LeptonicB_C: isLeptonicB_C_=1; break;
        case deep_ntuples::JetFlavor::GBB: isGBB_=1; break;
        case deep_ntuples::JetFlavor::BB: isBB_=1; break;
        case deep_ntuples::JetFlavor::C:  isC_=1; break;
        case deep_ntuples::JetFlavor::GCC: isGCC_=1; break;
        case deep_ntuples::JetFlavor::CC: isCC_=1; break;
        case deep_ntuples::JetFlavor::TAU: isTau_=1;break;
        case deep_ntuples::JetFlavor::G:  isG_=1; break;
        case deep_ntuples::JetFlavor::UD: isUD_=1; break;
        case deep_ntuples::JetFlavor::S:  isS_=1; break;
        default : isUndefined_=1; break;
        }
    }

    //truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    //// Note that jets with gluon->bb (cc) and x->bb (cc) are in the same categories
    isPhysB_=0; isPhysBB_=0; isPhysGBB_=0; isPhysC_=0; isPhysCC_=0;
    isPhysGCC_=0; isPhysUD_=0; isPhysS_=0; isPhysG_=0, isPhysLeptonicB_=0, isPhysLeptonicB_C_=0, isPhysUndefined_=0;
    isPhysTau_=0;
    if(jet.genJet()!=NULL){
        switch(deep_ntuples::jet_flavour(jet, gToBB, gToCC, neutrinosLepB, neutrinosLepB_C, alltaus_,true)) {
        case deep_ntuples::JetFlavor::UD: isPhysUD_=1; break;
        case deep_ntuples::JetFlavor::S:  isPhysS_=1; break;
        case deep_ntuples::JetFlavor::B:  isPhysB_=1; break;
        case deep_ntuples::JetFlavor::BB: isPhysBB_=1; break;
        case deep_ntuples::JetFlavor::GBB: isPhysGBB_=1; break;
        case deep_ntuples::JetFlavor::C:  isPhysC_=1; break;
        case deep_ntuples::JetFlavor::CC: isPhysCC_=1; break;
        case deep_ntuples::JetFlavor::GCC: isPhysGCC_=1; break;
        case deep_ntuples::JetFlavor::TAU: isPhysTau_=1;break;
        case deep_ntuples::JetFlavor::G:  isPhysG_=1; break;
        case deep_ntuples::JetFlavor::LeptonicB: isPhysLeptonicB_=1; break;
        case deep_ntuples::JetFlavor::LeptonicB_C: isPhysLeptonicB_C_=1; break;
        default : isPhysUndefined_=1; break;
        }
    }

    if(!jet.genJet()){//for data
        isUndefined_=1;isPhysUndefined_=1;
    }

    if(isUndefined_ && isPhysUndefined_) returnval=false; //skip event, if neither standard flavor definition nor physics definition fallback define a "proper flavor"

    pat::JetCollection h;

    jet_pt_ = jet.correctedJet("Uncorrected").pt();
    jet_eta_ = jet.eta();
    jet_phi_ = jet.phi();
    jet_corr_pt_ = jet.pt();
    jet_mass_ = jet.mass();
    jet_energy_ = jet.energy();

    genDecay_ = -1.;
    // jet pdgId implementation (using slimmedGenJetsFlavourInfos):
    gen_parton_pdgid_ = 0; gen_hadron_pdgid_ = 0; gen_hadron_pt_ = -1;
	for (const reco::JetFlavourInfoMatching & jetFlavourInfoMatching : *jetFlavourInfos) {
		if (deltaR(jet.p4(), jetFlavourInfoMatching.first->p4()) > 0.4) continue;
		gen_parton_pdgid_ = jetFlavourInfoMatching.second.getPartonFlavour();
		const reco::GenParticleRefVector & bHadrons = jetFlavourInfoMatching.second.getbHadrons();
		if (bHadrons.size()==0) continue;
		gen_hadron_pdgid_= bHadrons.at(0)->pdgId();
		gen_hadron_pt_= bHadrons.at(0)->pt();
		break;
	}
	
	// MP 30/09/2020: skip non b-jets:
	if(!isPhysBB_ && !isPhysB_ && !isPhysGBB_ && !isPhysLeptonicB_ && !isPhysLeptonicB_C_)
		returnval=false;
			

    try {
        reco::GenParticleRefVector Bhadrons_in_jet = jet.jetFlavourInfo().getbHadrons();

        if (Bhadrons_in_jet.size() > 0){ 

            for (unsigned int idx=0; idx<Bhadron_.size(); ++idx){

                reco::GenParticle bhad = Bhadron_[idx];

                bool bhad_is_in_jet = false;

                for (reco::GenParticleRefVector::const_iterator bhad_in_jet = Bhadrons_in_jet.begin(); bhad_in_jet!=Bhadrons_in_jet.end(); ++bhad_in_jet) {

                    //check if bhad is identical to bhad_in_jet
                    if ( (*bhad_in_jet)->pt() == bhad.pt() && (*bhad_in_jet)->eta() == bhad.eta()
                            && (*bhad_in_jet)->phi() == bhad.phi() && (*bhad_in_jet)->pdgId() == bhad.pdgId())              
                        bhad_is_in_jet = true;
                }
                if (bhad_is_in_jet){

                    if (Bhadron_daughter_[idx].vx()!=bhad.vx()){

                        float vx = Bhadron_daughter_[idx].vx() - bhad.vx();
                        float vy = Bhadron_daughter_[idx].vy() - bhad.vy();

                        float dxy = sqrt(vx*vx+vy*vy);
                        if (dxy > genDecay_)
                            genDecay_= dxy;
                    }
                    else if (genDecay_ < 0) 
                        genDecay_ = -0.1;
                }
            }
        }
    }
    catch (const cms::Exception &e){
        genDecay_ = -1.;
    }



    //https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2016
    try{
        float NHF  = jet.neutralHadronEnergyFraction();
        float NEMF = jet.neutralEmEnergyFraction();
        float CHF  = jet.chargedHadronEnergyFraction();
        //float MUF  = jet.muonEnergyFraction();
        float CEMF = jet.chargedEmEnergyFraction();
        float NumConst = jet.chargedMultiplicity()+jet.neutralMultiplicity();
        float NumNeutralParticles =jet.neutralMultiplicity();
        float CHM      = jet.chargedMultiplicity();

        jet_looseId_ = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(jet_eta_)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jet_eta_)>2.4) && abs(jet_eta_)<=2.7) ||
                (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 && abs(jet_eta_)>2.7 && abs(jet_eta_)<=3.0 ) ||
                (NEMF<0.90 && NumNeutralParticles>10 && abs(jet_eta_)>3.0 );
    }catch(const cms::Exception &e){
        jet_looseId_ = 1;
    }


    gen_pt_ =  0;
    Delta_gen_pt_ =  0;
    gen_pt_Recluster_=0;
    gen_pt_WithNu_=0;
    Delta_gen_pt_Recluster_=0;
    Delta_gen_pt_WithNu_=0;

    if(jet.genJet()){
        gen_pt_ =  jet.genJet()->pt();
        Delta_gen_pt_ =  jet.genJet()->pt()- jet_pt_;

        const edm::RefToBase<pat::Jet> patJetRef = coll->refAt(jetidx);
        reco::GenJetRef genjetRecluster = (*genJetMatchRecluster)[patJetRef];

        gen_pt_Recluster_ = 0.;
        if (genjetRecluster.isNonnull() && genjetRecluster.isAvailable()) {
            gen_pt_Recluster_ = genjetRecluster->pt();
        }
        reco::GenJetRef genjetWithNu = (*genJetMatchWithNu)[patJetRef];

        gen_pt_WithNu_ = 0.;
        if (genjetWithNu.isNonnull() && genjetWithNu.isAvailable()) {
            gen_pt_WithNu_ = genjetWithNu->pt();
        }

        Delta_gen_pt_Recluster_=gen_pt_Recluster_-jet.pt();
        Delta_gen_pt_WithNu_=gen_pt_WithNu_-jet.pt();
    }


    auto qgtuple=yuta::calcVariables(&jet);
    //(multiplicity, charged_multiplicity, neutral_multiplicity, ptD, axis1, axis2, pt_dr_log);

    y_multiplicity_=std::get<0>(qgtuple);
    y_charged_multiplicity_=std::get<1>(qgtuple);
    y_neutral_multiplicity_=std::get<2>(qgtuple);
    y_ptD_    =  std::get<3>(qgtuple);
    y_axis1_  =  std::get<4>(qgtuple);
    y_axis2_  =  std::get<5>(qgtuple);
    y_pt_dr_log_=std::get<6>(qgtuple);

   



    return returnval;
}
