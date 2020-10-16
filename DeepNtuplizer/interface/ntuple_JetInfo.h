/*
 * ntuple_JetInfo.h
 *
 *  Created on: 13 Feb 2017
 *      Author: jkiesele
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_

#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"

#include "ntuple_content.h"
#include "TRandom3.h"
#include <map>
#include <string>

/*
 * For global jet info such as eta, pt, gen info
 */
class ntuple_JetInfo: public ntuple_content{
public:
    ntuple_JetInfo():ntuple_content(),
    gluonReduction_(0),
    useherwcompat_matching_(false),
    isherwig_(false)
{}

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);

    //use either of these functions

    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);

    void setAxis2Token(edm::EDGetTokenT<edm::ValueMap<float> > axis2Token) {
        axis2Token_ = axis2Token;
    }

    void setMultToken(edm::EDGetTokenT<edm::ValueMap<int> > multToken) {
        multToken_ = multToken;
    }

    void setPtDToken(edm::EDGetTokenT<edm::ValueMap<float> > ptDToken) {
        ptDToken_ = ptDToken;
    }

    void setQglToken(edm::EDGetTokenT<edm::ValueMap<float> > qglToken) {
        qglToken_ = qglToken;
    }

    void setGenJetMatchReclusterToken(
            edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchReclusterToken) {
        genJetMatchReclusterToken_ = genJetMatchReclusterToken;
    }

    void setGenJetMatchWithNuToken(
            edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchWithNuToken) {
        genJetMatchWithNuToken_ = genJetMatchWithNuToken;
    }

    void setGenParticlesToken(edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken) {
        genParticlesToken_ = genParticlesToken;
    }

    void setJetInfoToken(edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken) {
        jetFlavourInfosToken_ = jetFlavourInfosToken;
    }
    
    void setMuonsToken(edm::EDGetTokenT<pat::MuonCollection> muonsToken) {
        muonsToken_ = muonsToken;
    }

    void setElectronsToken(edm::EDGetTokenT<pat::ElectronCollection> electronsToken) {
        electronsToken_ = electronsToken;
    }

    void setUseHerwigCompatibleMatching(const bool use){
        useherwcompat_matching_=use;
    }
    void setIsHerwig(const bool use){
        isherwig_=use;
    }

    //private:

    double                    jetPtMin_;
    double                    jetPtMax_;
    double                    jetAbsEtaMin_;
    double                    jetAbsEtaMax_;

    //Quark gluon likelihood
    edm::EDGetTokenT<edm::ValueMap<float>>   qglToken_;
    edm::EDGetTokenT<edm::ValueMap<float>>   ptDToken_;
    edm::EDGetTokenT<edm::ValueMap<float>>   axis2Token_;
    edm::EDGetTokenT<edm::ValueMap<int>>     multToken_;

    edm::Handle<edm::ValueMap<float>> qglHandle;
    edm::Handle<edm::ValueMap<float>> ptDHandle;
    edm::Handle<edm::ValueMap<float>> axis2Handle;
    edm::Handle<edm::ValueMap<int>> multHandle;


    edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchReclusterToken_;
    edm::EDGetTokenT<edm::Association<reco::GenJetCollection> > genJetMatchWithNuToken_;

    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

    edm::EDGetTokenT<pat::MuonCollection> muonsToken_;       
    edm::EDGetTokenT<pat::ElectronCollection> electronsToken_;

    edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatchRecluster;
    edm::Handle<edm::Association<reco::GenJetCollection> > genJetMatchWithNu;

    edm::Handle<reco::GenParticleCollection> genParticlesHandle;

    edm::Handle<pat::MuonCollection> muonsHandle;
    edm::Handle<pat::ElectronCollection> electronsHandle;
    
    edm::EDGetTokenT<reco::JetFlavourInfoMatchingCollection> jetFlavourInfosToken_;
    edm::Handle<reco::JetFlavourInfoMatchingCollection> jetFlavourInfos;


    TRandom3 TRandom_;
    float gluonReduction_;

    std::vector <reco::GenParticle> neutrinosLepB;
    std::vector <reco::GenParticle> neutrinosLepB_C;

    std::vector<reco::GenParticle> gToBB;
    std::vector<reco::GenParticle> gToCC;
    std::vector<reco::GenParticle> alltaus_;

    std::vector<reco::GenParticle> Bhadron_;
    std::vector<reco::GenParticle> Bhadron_daughter_;



    bool useherwcompat_matching_;
    bool isherwig_;

    /////////branches

    // labels (MC truth)
    // regressions pt, Deta, Dphi
    float gen_pt_;
    int gen_parton_pdgid_;
    int gen_hadron_pdgid_;
    float gen_hadron_pt_;
    float Delta_gen_pt_;
    //classification
    int isB_;
    int isGBB_;
    int isBB_;
    int isC_;
    int isGCC_;
    int isCC_;
    int isUD_;
    int isS_;
    int isG_;
    int isUndefined_;
    float genDecay_;
    int isLeptonicB_;
    int isLeptonicB_C_;
    int isTau_;

    //truth labeling with fallback to physics definition for light/gluon/undefined of standard flavor definition
    int isPhysB_;
    int isPhysGBB_;
    int isPhysBB_;
    int isPhysC_;
    int isPhysGCC_;
    int isPhysCC_;
    int isPhysUD_;
    int isPhysS_;
    int isPhysG_;
    int isPhysUndefined_;
    int isPhysLeptonicB_;
    int isPhysLeptonicB_C_;
    int isPhysTau_;

    // global variables
    float npv_;
    float ntrueInt_;
    float rho_;
    unsigned int event_no_;
    unsigned int jet_no_;

    // jet variables
    float jet_pt_;
    float jet_corr_pt_;
    float  jet_eta_;
    float  jet_phi_;
    float  jet_mass_;
    float  jet_energy_;

    float jet_looseId_;

    // quark/gluon
    float jet_qgl_;
    float QG_ptD_;
    float QG_axis2_;
    float QG_mult_;


    float y_multiplicity_;
    float y_charged_multiplicity_;
    float y_neutral_multiplicity_;
    float y_ptD_;
    float y_axis1_;
    float y_axis2_;
    float y_pt_dr_log_;

    static constexpr std::size_t max_num_lept = 5;
    float muons_isLooseMuon_[max_num_lept];
    float muons_isTightMuon_[max_num_lept];
    float muons_isSoftMuon_[max_num_lept];
    float muons_isHighPtMuon_[max_num_lept]; 
    float muons_pt_[max_num_lept]; 
    float muons_relEta_[max_num_lept]; 
    float muons_relPhi_[max_num_lept]; 
    float muons_energy_[max_num_lept]; 
    float muons_charge_[max_num_lept];
    float electrons_pt_[max_num_lept]; 
    float electrons_relEta_[max_num_lept]; 
    float electrons_relPhi_[max_num_lept]; 
    float electrons_energy_[max_num_lept];
    float electrons_charge_[max_num_lept];

    // lepton extra vars
    float muons_isGlobal_[max_num_lept];
    float muons_isStandAlone_[max_num_lept];
    float muons_jetDeltaR_[max_num_lept];
    float muons_numberOfMatchedStations_[max_num_lept];
    float muons_2dIP_[max_num_lept];
    float muons_2dIPSig_[max_num_lept];
    float muons_3dIP_[max_num_lept];
    float muons_3dIPSig_[max_num_lept];
    float muons_dxy_[max_num_lept];
    float muons_dxyError_[max_num_lept];
    float muons_dxySig_[max_num_lept];
    float muons_dz_[max_num_lept];
    float muons_dzError_[max_num_lept];
    float muons_numberOfValidPixelHits_[max_num_lept];
    float muons_numberOfpixelLayersWithMeasurement_[max_num_lept];
    float muons_numberOfstripLayersWithMeasurement_[max_num_lept];
    float muons_chi2_[max_num_lept];
    float muons_ndof_[max_num_lept];
    float muons_caloIso_[max_num_lept];
    float muons_ecalIso_[max_num_lept];
    float muons_hcalIso_[max_num_lept];
    float muons_sumPfChHadronPt_[max_num_lept];
    float muons_sumPfNeuHadronEt_[max_num_lept];
    float muons_Pfpileup_[max_num_lept];
    float muons_sumPfPhotonEt_[max_num_lept];
    float muons_sumPfChHadronPt03_[max_num_lept];
    float muons_sumPfNeuHadronEt03_[max_num_lept];
    float muons_Pfpileup03_[max_num_lept];
    float muons_sumPfPhotonEt03_[max_num_lept];
    float electrons_jetDeltaR_[max_num_lept];
    float electrons_EtFromCaloEn_[max_num_lept];
    float electrons_ecalDrivenSeed_[max_num_lept];
    float electrons_isEB_[max_num_lept];
    float electrons_isEE_[max_num_lept];
    float electrons_ecalEnergy_[max_num_lept];
    float electrons_isPassConversionVeto_[max_num_lept];
    float electrons_convDist_[max_num_lept];
    float electrons_convFlags_[max_num_lept];
    float electrons_convRadius_[max_num_lept];
    float electrons_3dIP_[max_num_lept];
    float electrons_3dIPSig_[max_num_lept];
    float electrons_2dIP_[max_num_lept];
    float electrons_2dIPSig_[max_num_lept];
    float electrons_sCseedEta_[max_num_lept];
    float electrons_eSeedClusterOverP_[max_num_lept];
    float electrons_eSeedClusterOverPout_[max_num_lept];
    float electrons_eSuperClusterOverP_[max_num_lept];
    float electrons_hadronicOverEm_[max_num_lept];
    float electrons_deltaEtaEleClusterTrackAtCalo_[max_num_lept];
    float electrons_deltaPhiEleClusterTrackAtCalo_[max_num_lept];
    float electrons_deltaEtaSeedClusterTrackAtCalo_[max_num_lept];
    float electrons_deltaPhiSeedClusterTrackAtCalo_[max_num_lept];
    float electrons_deltaEtaSeedClusterTrackAtVtx_[max_num_lept];
    float electrons_deltaEtaSuperClusterTrackAtVtx_[max_num_lept];
    float electrons_deltaPhiSuperClusterTrackAtVtx_[max_num_lept];
    float electrons_dxy_[max_num_lept];
    float electrons_dz_[max_num_lept];
    float electrons_nbOfMissingHits_[max_num_lept];
    float electrons_gsfCharge_[max_num_lept];
    float electrons_ndof_[max_num_lept];
    float electrons_chi2_[max_num_lept];
    float electrons_SC_energy_[max_num_lept];
    float electrons_SC_deta_[max_num_lept];
    float electrons_SC_dphi_[max_num_lept];
    float electrons_SC_et_[max_num_lept];
    float electrons_scPixCharge_[max_num_lept];
    float electrons_numberOfBrems_[max_num_lept];
    float electrons_fbrem_[max_num_lept];
    float electrons_sigmaEtaEta_[max_num_lept];
    float electrons_sigmaIetaIeta_[max_num_lept];
    float electrons_sigmaIphiIphi_[max_num_lept];
    float electrons_r9_[max_num_lept];
    float electrons_superClusterFbrem_[max_num_lept];
    float electrons_e5x5_[max_num_lept];
    float electrons_e5x5Rel_[max_num_lept];
    float electrons_e1x5Overe5x5_[max_num_lept];
    float electrons_e2x5MaxOvere5x5_[max_num_lept];
    float electrons_hcalOverEcal_[max_num_lept];
    float electrons_SC_eSuperClusterOverP_[max_num_lept];
    float electrons_neutralHadronIso_[max_num_lept];
    float electrons_photonIso_[max_num_lept];
    float electrons_puChargedHadronIso_[max_num_lept];
    float electrons_trackIso_[max_num_lept];
    float electrons_hcalDepth1OverEcal_[max_num_lept];
    float electrons_hcalDepth2OverEcal_[max_num_lept];
    float electrons_ecalPFClusterIso_[max_num_lept];
    float electrons_hcalPFClusterIso_[max_num_lept];
    float electrons_pfSumPhotonEt_[max_num_lept];
    float electrons_pfSumChargedHadronPt_[max_num_lept];
    float electrons_pfSumNeutralHadronEt_[max_num_lept];
    float electrons_pfSumPUPt_[max_num_lept];
    float electrons_dr04TkSumPt_[max_num_lept];
    float electrons_dr04EcalRecHitSumEt_[max_num_lept];
    float electrons_dr04HcalDepth1TowerSumEt_[max_num_lept];
    float electrons_dr04HcalDepth1TowerSumEtBc_[max_num_lept];
    float electrons_dr04HcalDepth2TowerSumEt_[max_num_lept];
    float electrons_dr04HcalDepth2TowerSumEtBc_[max_num_lept];
    float electrons_dr04HcalTowerSumEt_[max_num_lept];
    float electrons_dr04HcalTowerSumEtBc_[max_num_lept];
            
    float muons_number_ = 0; int mn_ = 0;
    float electrons_number_ = 0; int en_ = 0;

    float gen_pt_Recluster_;
    float gen_pt_WithNu_;
    float Delta_gen_pt_Recluster_;
    float Delta_gen_pt_WithNu_;
    std::map<std::string, float> discriminators_;
};


#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_JETINFO_H_ */
