#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexPrimitives/interface/ConvertToFromReco.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "RecoBTag/FeatureTools/interface/sorting_modules.h"

#include <string> 

using namespace btagbtvdeep;

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

// To store the gen info to get the truth flavour of the jet
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "../interface/helpers.h"

template<typename T>
class JetConstituentTableProducerCustomTagger : public edm::stream::EDProducer<> {
public:
  explicit JetConstituentTableProducerCustomTagger(const edm::ParameterSet &);
  ~JetConstituentTableProducerCustomTagger() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  
  int jet_flavour(const pat::Jet& jet,
    const std::vector<reco::GenParticle>& gToBB,
    const std::vector<reco::GenParticle>& gToCC,
    const std::vector<reco::GenParticle>& neutrinosLepB,
    const std::vector<reco::GenParticle>& neutrinosLepB_C,
    const std::vector<reco::GenParticle>& alltaus,
    bool usePhysForLightAndUndefined) { 

      int hflav = abs(jet.hadronFlavour());
      int pflav = abs(jet.partonFlavour());
      int physflav = 0;
      if ( !( jet.genJet() ) ){
          if (pflav == 0) return 999;
          else return 1000;
      }
      if(jet.genParton()) physflav=abs(jet.genParton()->pdgId());
      std::size_t nbs = jet.jetFlavourInfo().getbHadrons().size();
      std::size_t ncs = jet.jetFlavourInfo().getcHadrons().size();
  
      unsigned int nbFromGSP(0);
      for (reco::GenParticle p : gToBB) {
          double dr2(reco::deltaR2(jet, p));
          if (dr2 < jetR_ * jetR_) ++nbFromGSP;
      }
  
      unsigned int ncFromGSP(0);
      for (reco::GenParticle p : gToCC) {
          double dr2(reco::deltaR2(jet, p));
          if (dr2 < jetR_ * jetR_) ++ncFromGSP;
      }
  
      //std::cout << " jet pt = " << jet.pt() << " hfl = " << hflav << " pfl = " << pflav << " genpart = " << physflav
              //  << " nbFromGSP = " << nbFromGSP << " ncFromGSP = " << ncFromGSP
      //  << " nBhadrons " << nbs << " nCHadrons " << ncs << std::endl;
      if(hflav == 5) { //B jet
          if(nbs > 1) {
              if (nbFromGSP > 0) return 511;
              else return 510;
          }
          else if(nbs == 1) {
              for (std::vector<reco::GenParticle>::const_iterator it = neutrinosLepB.begin(); it != neutrinosLepB.end(); ++it){
                  if(reco::deltaR(it->eta(),it->phi(),jet.eta(),jet.phi()) < 0.4) {
                      return 520;
                  }
              }
              for (std::vector<reco::GenParticle>::const_iterator it = neutrinosLepB_C.begin(); it != neutrinosLepB_C.end(); ++it){
                  if(reco::deltaR(it->eta(),it->phi(),jet.eta(),jet.phi()) < 0.4) {
                      return 521;
                  }
              }
              return 500;
          }
          else {
              if(usePhysForLightAndUndefined){
                  if(physflav == 21) return 0;
                  else if(physflav == 3) return 2;
                  else if(physflav == 2 || physflav ==1) return 1;
                  else return 1000;
              }
              else return 1000;
          }
      }
      else if(hflav == 4) { //C jet
          if (ncs > 1) {
              if (ncFromGSP > 0) return 411;
              else return 410;
          }
          else return 400;
      }
      else { //not a heavy jet
          if(alltaus.size()>0){ //check for tau in a simplistic way
              bool ishadrtaucontained=true;
              for(const auto& p:alltaus){
                  size_t ndau=p.numberOfDaughters();
                  for(size_t i=0;i<ndau;i++){
                      const reco::Candidate* dau=p.daughter(i);
                      int daupid=std::abs(dau->pdgId());
                      if(daupid == 13 || daupid == 11){
                          ishadrtaucontained=false;
                          break;
                      }
                      if(daupid != 12 && daupid!=14 && daupid!=16 &&
                              reco::deltaR2(*dau,jet) > jetR_ * jetR_){
                          ishadrtaucontained=false;
                          break;
                      }
                  }
              }
              if(ishadrtaucontained) return 600;
          }
          if(std::abs(pflav) == 4 || std::abs(pflav) == 5 || nbs || ncs) {
              if(usePhysForLightAndUndefined){
                  if(physflav == 21) return 0;
                  else if(physflav == 3) return 2;
                  else if(physflav == 2 || physflav ==1) return 1;
                  else return 1000;
              }
              else return 1000;
          }
          else if(usePhysForLightAndUndefined){
              if(physflav == 21) return 0;
              else if(physflav == 3) return 2;
              else if(physflav == 2 || physflav ==1) return 1;
              else return 1000;
          }
          else {
              if(pflav == 21) return 0;
              else if(pflav == 3) return 2;
              else if(pflav == 2 || pflav ==1) return 1;
              else return 1000;
          }
    }
        
    }
  
private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  typedef reco::VertexCollection VertexCollection;
  //=====
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;

  //const std::string name_;
  //const std::string name_;
  //const std::string nameSV_;
  const std::string nameCustomTagger_;
  //const std::string idx_name_;
  //const std::string idx_nameSV_;
  const std::string idx_nameCustomTagger_;
  const std::string storeAK4Truth_;    
  const bool readBtag_;
  const double jet_radius_;
  const bool add_CustomTagger_noclip_;

  const edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<edm::View<T>> jet_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  edm::EDGetTokenT<reco::CandidateView> cand_token_;
  edm::EDGetTokenT<SVCollection> sv_token_;

  edm::Handle<VertexCollection> vtxs_;
  edm::Handle<reco::CandidateView> cands_;
  edm::Handle<SVCollection> svs_;
  edm::ESHandle<TransientTrackBuilder> track_builder_;

  const reco::Vertex *pv_ = nullptr;
    
  const float min_candidate_pt_ = 0.95;
  
  constexpr static unsigned n_cpf_ = 25;
  constexpr static unsigned n_npf_ = 25;
  constexpr static unsigned n_sv_ = 5;
    
  constexpr static double jetR_ = 0.4;    
    
  constexpr static bool usePhysForLightAndUndefined = false;  
  
};

//
// constructors and destructor
//
template< typename T>
JetConstituentTableProducerCustomTagger<T>::JetConstituentTableProducerCustomTagger(const edm::ParameterSet &iConfig)
    : //name_(iConfig.getParameter<std::string>("name")),
      //nameSV_(iConfig.getParameter<std::string>("nameSV")),
      nameCustomTagger_(iConfig.getParameter<std::string>("nameCustomTagger")),
      //idx_name_(iConfig.getParameter<std::string>("idx_name")),
      //idx_nameSV_(iConfig.getParameter<std::string>("idx_nameSV")),
      idx_nameCustomTagger_(iConfig.getParameter<std::string>("idx_nameCustomTagger")),
      storeAK4Truth_(iConfig.getParameter<std::string>("storeAK4Truth")),
      readBtag_(iConfig.getParameter<bool>("readBtag")),
      jet_radius_(iConfig.getParameter<double>("jet_radius")),
      add_CustomTagger_noclip_(iConfig.getParameter<bool>("add_CustomTagger_noclip")),
      genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"))),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))){
  // for some reason, the tokens need to come after the other parameters (i.e. keeping the order from above)
  //produces<nanoaod::FlatTable>(name_);
  //produces<nanoaod::FlatTable>(name_);
  //produces<nanoaod::FlatTable>(nameSV_);
  produces<nanoaod::FlatTable>(nameCustomTagger_);
  //produces<std::vector<reco::CandidatePtr>>();
}

template< typename T>
JetConstituentTableProducerCustomTagger<T>::~JetConstituentTableProducerCustomTagger() {}

template< typename T>
void JetConstituentTableProducerCustomTagger<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  //auto outCands = std::make_unique<std::vector<reco::CandidatePtr>>();
  //auto outSVs = std::make_unique<std::vector<const reco::VertexCompositePtrCandidate *>> ();
  std::vector<int> jetIdx_dj;
  
    
  auto jets = iEvent.getHandle(jet_token_);
  iEvent.getByToken(vtx_token_, vtxs_);
  iEvent.getByToken(cand_token_, cands_);
  iEvent.getByToken(sv_token_, svs_);

  if(readBtag_){
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", track_builder_);
  }
    
  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  iEvent.getByToken(genParticlesToken_, genParticlesHandle);
    
  std::vector<reco::GenParticle> neutrinosLepB;
  std::vector<reco::GenParticle> neutrinosLepB_C;
  
  std::vector<reco::GenParticle> gToBB;
  std::vector<reco::GenParticle> gToCC;
  std::vector<reco::GenParticle> alltaus;
    
    
  unsigned nJets = jets->size();
    
  std::vector<int> jet_N_CPFCands(nJets);
  std::vector<int> jet_N_NPFCands(nJets);
  std::vector<int> jet_N_PVs(nJets); 
  std::vector<int> jet_N_SVs(nJets); 
    
    
  std::vector<unsigned> jet_FlavSplit(nJets); 
    
    
    
    
    
  // should default to 0 if less than 25 cpf with information
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackEtaRel_nCpf(n_cpf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPtRel_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPPar_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackDeltaR_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPParRatio_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip2dVal_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip2dSig_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip3dVal_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip3dSig_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackJetDistVal_nCpf(n_cpf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Cpfcan_ptrel_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_drminsv_nCpf(n_cpf_, std::vector<float>(nJets));
  std::vector<std::vector<int>> Cpfcan_VTX_ass_nCpf(n_cpf_, std::vector<int>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_puppiw_nCpf(n_cpf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Cpfcan_chi2_nCpf(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<int>> Cpfcan_quality_nCpf(n_cpf_, std::vector<int>(nJets));
  // no clip versions
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackEtaRel_nCpf_noclip(n_cpf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPtRel_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPPar_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackDeltaR_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackPParRatio_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip2dVal_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip2dSig_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip3dVal_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackSip3dSig_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_BtagPf_trackJetDistVal_nCpf_noclip(n_cpf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Cpfcan_ptrel_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Cpfcan_drminsv_nCpf_noclip(n_cpf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Cpfcan_chi2_nCpf_noclip(n_cpf_, std::vector<float>(nJets)); 
    
    
  // should default to 0 if less than 25 npf with information
  std::vector<std::vector<float>> Npfcan_ptrel_nNpf(n_npf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Npfcan_deltaR_nNpf(n_npf_, std::vector<float>(nJets)); 
  std::vector<std::vector<int>> Npfcan_isGamma_nNpf(n_npf_, std::vector<int>(nJets)); 
  std::vector<std::vector<float>> Npfcan_HadFrac_nNpf(n_npf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Npfcan_drminsv_nNpf(n_npf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Npfcan_puppiw_nNpf(n_npf_, std::vector<float>(nJets)); 
  // no clip versions
  std::vector<std::vector<float>> Npfcan_ptrel_nNpf_noclip(n_npf_, std::vector<float>(nJets));
  std::vector<std::vector<float>> Npfcan_deltaR_nNpf_noclip(n_npf_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> Npfcan_drminsv_nNpf_noclip(n_npf_, std::vector<float>(nJets)); 
    

  // should default to 0 if less than four SVs with information
  std::vector<std::vector<float>> sv_mass_nSV(n_sv_, std::vector<float>(nJets));
  std::vector<std::vector<float>> sv_pt_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<int>> sv_ntracks_nSV(n_sv_, std::vector<int>(nJets)); 
  std::vector<std::vector<float>> sv_chi2_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_normchi2_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_dxy_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_dxysig_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_d3d_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_d3dsig_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_costhetasvpv_nSV(n_sv_, std::vector<float>(nJets));
  std::vector<std::vector<float>> sv_ptrel_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_phirel_nSV(n_sv_, std::vector<float>(nJets));
  std::vector<std::vector<float>> sv_deltaR_nSV(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_enratio_nSV(n_sv_, std::vector<float>(nJets));
  // no clip versions
  std::vector<std::vector<float>> sv_normchi2_nSV_noclip(n_sv_, std::vector<float>(nJets));
  std::vector<std::vector<float>> sv_deltaR_nSV_noclip(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_dxysig_nSV_noclip(n_sv_, std::vector<float>(nJets)); 
  std::vector<std::vector<float>> sv_d3dsig_nSV_noclip(n_sv_, std::vector<float>(nJets));
  
  if (storeAK4Truth_ == "yes") { 

    neutrinosLepB.clear();
    neutrinosLepB_C.clear();
    gToBB.clear();
    gToCC.clear();
    alltaus.clear();  

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
                      && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR2(*d0, *d1) < jetR_ * jetR_) gToBB.push_back(gen) ;
              if ( std::abs(d0->pdgId()) == 4 && std::abs(d1->pdgId()) == 4
                      && d0->pdgId()*d1->pdgId() < 0 && reco::deltaR2(*d0, *d1) < jetR_ * jetR_) gToCC.push_back(gen) ;
          }
      }

      if(id == 15 && false){
          alltaus.push_back(gen);
      }

    }
      
  }
    
  //std::cout << "Successfully initialized all vectors." << std::endl;    
    
  //std::cout << "In producer: add_CustomTagger_noclip_ ? " << add_CustomTagger_noclip_ << std::endl;  
  
  //std::cout << "Start jet loop." << std::endl;  
    
  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    const auto &jet = jets->at(i_jet);
    math::XYZVector jet_dir = jet.momentum().Unit();
    GlobalVector jet_ref_track_dir(jet.px(), jet.py(), jet.pz());
    VertexDistance3D vdist;

    pv_ = &vtxs_->at(0);
    jet_N_PVs[i_jet] = vtxs_->size();  

    //////////////////////
    // Secondary Vertices
    std::vector<const reco::VertexCompositePtrCandidate *> jetSVs;
    std::vector<const reco::VertexCompositePtrCandidate *> allSVs;  
      
    if (storeAK4Truth_ == "yes") { 
      jet_FlavSplit[i_jet] =  jet_flavour(jet, gToBB, gToCC, neutrinosLepB, neutrinosLepB_C, alltaus, usePhysForLightAndUndefined);  
    }
    
    jetIdx_dj.push_back(i_jet);
      
    jet_N_SVs[i_jet] = 0;
    for (const auto &sv : *svs_) {
      // Factor in cuts in NanoAOD for indexing
      Measurement1D dl= vdist.distance(vtxs_->front(), VertexState(RecoVertex::convertPos(sv.position()),RecoVertex::convertError(sv.error())));
      if(dl.significance() > 3){
        allSVs.push_back(&sv);
      }
      if (reco::deltaR2(sv, jet) < jet_radius_ * jet_radius_) {
        jetSVs.push_back(&sv);
        jet_N_SVs[i_jet]++;
      }
    }
      
      
    // sort by dxy significance
    std::sort(jetSVs.begin(),
              jetSVs.end(),
              [&](const reco::VertexCompositePtrCandidate *sva, const reco::VertexCompositePtrCandidate *svb) {
                return sv_vertex_comparator(*sva, *svb, *pv_);
              });
    
    // counter to get flat info per jet for SVs
    unsigned i_sv_in_jet = 0;  

    
    for (const auto &sv : jetSVs) {

      if (readBtag_ && !vtxs_->empty()) {
          
        if (i_sv_in_jet < n_sv_) {   
          sv_mass_nSV[i_sv_in_jet][i_jet] = sv->mass();
          sv_pt_nSV[i_sv_in_jet][i_jet] = sv->pt();
          sv_ntracks_nSV[i_sv_in_jet][i_jet] = sv->numberOfDaughters();
          sv_chi2_nSV[i_sv_in_jet][i_jet] = sv->vertexChi2();
          sv_normchi2_nSV[i_sv_in_jet][i_jet] = catch_infs_and_bound(sv->vertexChi2() / sv->vertexNdof(), 1000, -1000, 1000);
          const auto& dxy_meas = vertexDxy(*sv, *pv_);
          sv_dxy_nSV[i_sv_in_jet][i_jet] = dxy_meas.value();
          sv_dxysig_nSV[i_sv_in_jet][i_jet] = catch_infs_and_bound(dxy_meas.value() / dxy_meas.error(), 0, -1, 800);
          const auto& d3d_meas = vertexD3d(*sv, *pv_);
          sv_d3d_nSV[i_sv_in_jet][i_jet] = d3d_meas.value();
          sv_d3dsig_nSV[i_sv_in_jet][i_jet] = catch_infs_and_bound(d3d_meas.value() / d3d_meas.error(), 0, -1, 800);
          sv_costhetasvpv_nSV[i_sv_in_jet][i_jet] = vertexDdotP(*sv, *pv_);
          // Jet related
          sv_ptrel_nSV[i_sv_in_jet][i_jet] = sv->pt() / jet.pt();
          sv_phirel_nSV[i_sv_in_jet][i_jet] = reco::deltaPhi(*sv, jet);
          sv_deltaR_nSV[i_sv_in_jet][i_jet] = catch_infs_and_bound(std::fabs(reco::deltaR(*sv, jet_dir)) - 0.5, 0, -2, 0);
          sv_enratio_nSV[i_sv_in_jet][i_jet] = sv->energy() / jet.energy();
            
          // no clip
          if (add_CustomTagger_noclip_) {
            sv_normchi2_nSV_noclip[i_sv_in_jet][i_jet] = sv->vertexChi2() / sv->vertexNdof();
            sv_dxysig_nSV_noclip[i_sv_in_jet][i_jet] = dxy_meas.value() / dxy_meas.error();
            sv_d3dsig_nSV_noclip[i_sv_in_jet][i_jet] = d3d_meas.value() / d3d_meas.error();
            sv_deltaR_nSV_noclip[i_sv_in_jet][i_jet] = std::fabs(reco::deltaR(*sv, jet_dir));
          }
        } else {
                continue;
        }
          
      } 
      i_sv_in_jet++;
    }

    //std::cout << "Successfully filled SV info, start sorting PF cands now." << std::endl;    
      
    // PF Cands    
    std::vector<reco::CandidatePtr> const & daughters = jet.daughterPtrVector();

    const auto& svs_unsorted = *svs_;  
      
      
      
    std::vector<btagbtvdeep::SortingClass<size_t>> c_sorted, n_sorted;
      
      
      
    // first time looping over all pf candidates
    //     to fill sorted indices and get a connection back to the old indices
    
    jet_N_CPFCands[i_jet] = 0;  
    jet_N_NPFCands[i_jet] = 0;  
      
    for (unsigned int i = 0; i < jet.numberOfDaughters(); i++) {
      auto cand = jet.daughter(i);
      //if ( cand.isNull() ) continue;
      //auto const *packedCand = dynamic_cast <pat::PackedCandidate const *>(cand.get());
      //if ( packedCand == nullptr ) continue;
      if (cand) {
/*
        // not used currently as it was not part of the JetConstituentTableProducer at first, only for DeepFlavourTagInfoProducer
        // candidates under 950MeV (configurable) are not considered
        // might change if we use also white-listing
        if (cand->pt() < min_candidate_pt_)
          continue;
*/
        if (cand->charge() != 0) {
          //auto& trackinfo = trackinfos.emplace(i, track_builder).first->second;
          //trackinfo.buildTrackInfo(cand, jet_dir, jet_ref_track_dir, *pv_); // *pv_ is an atlernative to vtxs_->at(0)
          btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
            // similar to https://github.com/cms-sw/cmssw/blob/master/RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h
          trkinfo.buildTrackInfo(cand, jet_dir, jet_ref_track_dir, vtxs_->at(0));
          c_sorted.emplace_back(
              i, trkinfo.getTrackSip2dSig(), -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet.pt());
          jet_N_CPFCands[i_jet]++;
        } else {
          n_sorted.emplace_back(i, -1, -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet.pt());
          jet_N_NPFCands[i_jet]++;
        }
      }
    }
    
      
    // sort collections (open the black-box if you please)
    std::sort(c_sorted.begin(), c_sorted.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);
    std::sort(n_sorted.begin(), n_sorted.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);

    std::vector<size_t> c_sortedindices, n_sortedindices;

    // this puts 0 everywhere and the right position in ind
    c_sortedindices = btagbtvdeep::invertSortingVector(c_sorted);
    n_sortedindices = btagbtvdeep::invertSortingVector(n_sorted);
  
      
    //std::cout << "Start looping over PF cands to fill info." << std::endl;    
      
    int i_pf_in_jet = 0;
    for (const auto &cand : daughters) {
      
      auto candPtrs = cands_->ptrs();
      auto candInNewList = std::find( candPtrs.begin(), candPtrs.end(), cand );
      if ( candInNewList == candPtrs.end() ) {
        //std::cout << "Cannot find candidate : " << cand.id() << ", " << cand.key() << ", pt = " << cand->pt() << std::endl;
        continue;
      }
/*
      outCands->push_back(cand);
      jetIdx_pf.push_back(i_jet);
      pfcandIdx.push_back(candInNewList - candPtrs.begin());
      //cand_pt.push_back(cand->pt());
*/
      if (readBtag_ && !vtxs_->empty()) {
        if ( cand.isNull() ) continue;
        auto const *packedCand = dynamic_cast <pat::PackedCandidate const *>(cand.get());
        if ( packedCand == nullptr ) continue;
        float drminpfcandsv = btagbtvdeep::mindrsvpfcand(svs_unsorted, &(*cand));
        if ( packedCand->charge() != 0 ) {
            
            // is charged candidate
            auto entry_sorted = c_sortedindices.at(i_pf_in_jet);
            //std::cout << "Current candidate is " << i_pf_in_jet << " and entry_sorted = c_sortedindices.at(i_pf_in_jet) = " << entry_sorted << std::endl; 
            // need only the first 25 cpfs for CustomTagger
            if (entry_sorted >= n_cpf_) {
                continue;
            }
            
            Cpfcan_puppiw_nCpf[entry_sorted][i_jet] = packedCand->puppiWeight();
            Cpfcan_VTX_ass_nCpf[entry_sorted][i_jet] = packedCand->pvAssociationQuality();
            Cpfcan_drminsv_nCpf[entry_sorted][i_jet] = catch_infs_and_bound(drminpfcandsv, 0, -1. * jet_radius_, 0, -1. * jet_radius_);
            Cpfcan_ptrel_nCpf[entry_sorted][i_jet] = catch_infs_and_bound(packedCand->pt() / jet.pt(), 0, -1, 0, -1);
            
            // no clip
            if (add_CustomTagger_noclip_) {
              Cpfcan_drminsv_nCpf_noclip[entry_sorted][i_jet] = drminpfcandsv;
              Cpfcan_ptrel_nCpf_noclip[entry_sorted][i_jet] = packedCand->pt() / jet.pt();
            }
            
            if ( packedCand && packedCand->hasTrackDetails()){
              const auto& pseudo_track = packedCand->pseudoTrack();
              Cpfcan_chi2_nCpf[entry_sorted][i_jet] = catch_infs_and_bound(pseudo_track.normalizedChi2(), 300, -1, 300);
                
              // this returns the quality enum not a mask.
              Cpfcan_quality_nCpf[entry_sorted][i_jet] = pseudo_track.qualityMask();  
                
                
              btagbtvdeep::TrackInfoBuilder trkinfo(track_builder_);
                // similar to https://github.com/cms-sw/cmssw/blob/master/RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h
              trkinfo.buildTrackInfo(&(*packedCand), jet_dir, jet_ref_track_dir, vtxs_->at(0));
              
              
              // with clip
              Cpfcan_BtagPf_trackEtaRel_nCpf[entry_sorted][i_jet]     = catch_infs_and_bound(trkinfo.getTrackEtaRel(), 0, -5, 15);
              Cpfcan_BtagPf_trackPtRel_nCpf[entry_sorted][i_jet]      = catch_infs_and_bound(trkinfo.getTrackPtRel(), 0, -1, 4);
              Cpfcan_BtagPf_trackPPar_nCpf[entry_sorted][i_jet]       = catch_infs_and_bound(trkinfo.getTrackPPar(), 0, -1e5, 1e5);
              Cpfcan_BtagPf_trackDeltaR_nCpf[entry_sorted][i_jet]     = catch_infs_and_bound(trkinfo.getTrackDeltaR(), 0, -5, 5);
              Cpfcan_BtagPf_trackPParRatio_nCpf[entry_sorted][i_jet]  = catch_infs_and_bound(trkinfo.getTrackPParRatio(), 0, -10, 100);
              Cpfcan_BtagPf_trackSip2dVal_nCpf[entry_sorted][i_jet]   = catch_infs_and_bound(trkinfo.getTrackSip2dVal(), 0, -1, 70);
              Cpfcan_BtagPf_trackSip2dSig_nCpf[entry_sorted][i_jet]   = catch_infs_and_bound(trkinfo.getTrackSip2dSig(), 0, -1, 4e4);
              Cpfcan_BtagPf_trackSip3dVal_nCpf[entry_sorted][i_jet]   = catch_infs_and_bound(trkinfo.getTrackSip3dVal(), 0, -1, 1e5);
              Cpfcan_BtagPf_trackSip3dSig_nCpf[entry_sorted][i_jet]   = catch_infs_and_bound(trkinfo.getTrackSip3dSig(), 0, -1, 4e4);
              Cpfcan_BtagPf_trackJetDistVal_nCpf[entry_sorted][i_jet] = catch_infs_and_bound(trkinfo.getTrackJetDistVal(), 0, -20, 1);  
                
              // no clip was default inside JetConstituentTableProducer, but CustomTagger seems to use the clipped versions
              if (add_CustomTagger_noclip_) {
                Cpfcan_BtagPf_trackEtaRel_nCpf_noclip[entry_sorted][i_jet]     = trkinfo.getTrackEtaRel();
                Cpfcan_BtagPf_trackPtRel_nCpf_noclip[entry_sorted][i_jet]      = trkinfo.getTrackPtRel();
                Cpfcan_BtagPf_trackPPar_nCpf_noclip[entry_sorted][i_jet]       = trkinfo.getTrackPPar();
                Cpfcan_BtagPf_trackDeltaR_nCpf_noclip[entry_sorted][i_jet]     = trkinfo.getTrackDeltaR();
                Cpfcan_BtagPf_trackPParRatio_nCpf_noclip[entry_sorted][i_jet]  = trkinfo.getTrackPParRatio();
                Cpfcan_BtagPf_trackSip2dVal_nCpf_noclip[entry_sorted][i_jet]   = trkinfo.getTrackSip2dVal();
                Cpfcan_BtagPf_trackSip2dSig_nCpf_noclip[entry_sorted][i_jet]   = trkinfo.getTrackSip2dSig();
                Cpfcan_BtagPf_trackSip3dVal_nCpf_noclip[entry_sorted][i_jet]   = trkinfo.getTrackSip3dVal();
                Cpfcan_BtagPf_trackSip3dSig_nCpf_noclip[entry_sorted][i_jet]   = trkinfo.getTrackSip3dSig();
                Cpfcan_BtagPf_trackJetDistVal_nCpf_noclip[entry_sorted][i_jet] = trkinfo.getTrackJetDistVal();
                
                Cpfcan_chi2_nCpf_noclip[entry_sorted][i_jet] = pseudo_track.normalizedChi2();
              }
                
              //c_pf_features.btagPf_trackPtRatio    = catch_infs_and_bound(track_info.getTrackPtRatio(), 0, -1, 10);
                
              //Cpfcan_ptrel.push_back(catch_infs_and_bound(cand->pt() / jet.pt(), 0, -1, 0, -1));
              //Cpfcan_drminsv.push_back(catch_infs_and_bound(mindrsvpfcand(svs_unsorted, &(*cand), 0.4), 0, -1. * jet_radius_, 0, -1. * jet_radius_));
            } else {
                    // default negative chi2 and loose track if notTrackDetails
                    Cpfcan_chi2_nCpf[entry_sorted][i_jet] = catch_infs_and_bound(-1, 300, -1, 300);
                    Cpfcan_quality_nCpf[entry_sorted][i_jet] = (1 << reco::TrackBase::loose);
                    
                    
                    // no clip
                    if (add_CustomTagger_noclip_) {
                      Cpfcan_chi2_nCpf_noclip[entry_sorted][i_jet] = -1;
                    }
           }
        } else {
            
            
            // is neutral candidate
            auto entry_sorted = n_sortedindices.at(i_pf_in_jet);
            //std::cout << "Current candidate is " << i_pf_in_jet << " and entry_sorted = n_sortedindices.at(i_pf_in_jet) = " << entry_sorted << std::endl; 
            // need only the first 25 cpfs for CustomTagger
            if (entry_sorted >= n_npf_) {
                continue;
            }
            
            Npfcan_puppiw_nNpf[entry_sorted][i_jet] = packedCand->puppiWeight();
            Npfcan_HadFrac_nNpf[entry_sorted][i_jet] = packedCand->hcalFraction();
            if (std::abs(packedCand->pdgId()) == 22)
              Npfcan_isGamma_nNpf[entry_sorted][i_jet] = 1;
            // catch
            Npfcan_deltaR_nNpf[entry_sorted][i_jet] = catch_infs_and_bound(reco::deltaR(*packedCand, jet), 0, -0.6, 0, -0.6);
            Npfcan_drminsv_nNpf[entry_sorted][i_jet] = catch_infs_and_bound(drminpfcandsv, 0, -1. * jet_radius_, 0, -1. * jet_radius_);
            Npfcan_ptrel_nNpf[entry_sorted][i_jet] = catch_infs_and_bound(packedCand->pt() / jet.pt(), 0, -1, 0, -1);
            // for all cases where catch_infs_and_bound appears, also store the raw quantities without any modifications (no clip)
            if (add_CustomTagger_noclip_) {
              Npfcan_deltaR_nNpf_noclip[entry_sorted][i_jet] = reco::deltaR(*packedCand, jet);
              Npfcan_drminsv_nNpf_noclip[entry_sorted][i_jet] = drminpfcandsv;
              Npfcan_ptrel_nNpf_noclip[entry_sorted][i_jet] = packedCand->pt() / jet.pt();
            }
        }
      }
      i_pf_in_jet++;
    }  // end jet loop
  }

    
  //std::cout << "Successfully made it through the jet loop, now starting filling the columns into the table." << std::endl; 
    
  // CustomTaggerInputs table
  auto customTable = std::make_unique<nanoaod::FlatTable>(jetIdx_dj.size(), nameCustomTagger_, false, true);
  customTable->addColumn<int>("CustomTagger_jetIdx", jetIdx_dj, "Index of the parent jet", nanoaod::FlatTable::IntColumn);
    
    
  // customTable->addColumn<int>("CustomTagger_nCpfcand", jet_N_CPFCands, "Number of charged PF candidates in the jet", nanoaod::FlatTable::IntColumn);
  // customTable->addColumn<int>("CustomTagger_nNpfcand", jet_N_NPFCands, "Number of neutral PF candidates in the jet", nanoaod::FlatTable::IntColumn);
  // customTable->addColumn<int>("CustomTagger_nsv", jet_N_SVs, "Number of secondary vertices in the jet", nanoaod::FlatTable::IntColumn);
    
  customTable->addColumn<int>("CustomTagger_nCpfcand",
                          jet_N_CPFCands,
                          "Number of charged PF candidates in the jet",
                          nanoaod::FlatTable::IntColumn);
  customTable->addColumn<int>("CustomTagger_nNpfcand",
                          jet_N_NPFCands,
                          "Number of neutral PF candidates in the jet",
                          nanoaod::FlatTable::IntColumn);
  customTable->addColumn<int>("CustomTagger_nsv",
                          jet_N_SVs,
                          "Number of secondary vertices in the jet",
                          nanoaod::FlatTable::IntColumn);
  customTable->addColumn<int>("CustomTagger_npv",
                          jet_N_PVs,
                          "Number of primary vertices",
                          nanoaod::FlatTable::IntColumn);
    
  if (storeAK4Truth_ == "yes") { 
      //std::cout << "Start filling table with truth info" << std::endl;
      customTable->addColumn<int>("FlavSplit",
                              jet_FlavSplit,
                              "Flavour of the jet, numerical codes: "
                              "isG: 0, "
                              "isUD: 1, "
                              "isS: 2, "
                              "isC: 400, "
                              "isCC: 410, "
                              "isGCC: 411, "
                              "isB: 500, "
                              "isBB: 510, "
                              "isGBB: 511, "
                              "isLeptonicB: 520, "
                              "isLeptonicB_C: 521, "
                              "isTAU: 600, "
                              "isPU: 999,"
                              "isUndefined: 1000. "
                              "May be combined to form coarse labels for tagger training and flavour dependent attacks using the loss surface.",
                              nanoaod::FlatTable::IntColumn);
  }
    
  std::string input_name;
  std::string description;  
  for (unsigned int p = 0; p < n_cpf_; p++) {
      auto s = std::to_string(p);
      
      
      // ============================================================== Cpfs ===================================================================
      input_name = "CustomTagger_Cpfcan_puppiw_" + s;
      description = "charged candidate PUPPI weight of the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_puppiw_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_VTX_ass_" + s;
      description = "integer flag that indicates whether the track was used in the primary vertex fit for the " + s + ". cpf";
      customTable->addColumn<int>(input_name, Cpfcan_VTX_ass_nCpf[p], description, nanoaod::FlatTable::IntColumn, 10);
      input_name = "CustomTagger_Cpfcan_drminsv_" + s;
      description = "track pseudoangular distance from the closest secondary vertex of the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_drminsv_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_ptrel_" + s;
      description = "fraction of the jet momentum carried by the track for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_ptrel_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_quality_" + s;
      description = "integer flag which indicates the quality of the fitted track, based on number of detector hits used for the reconstruction as well as the overall chi2 of the charged track fit for the " + s + ". cpf";
      customTable->addColumn<int>(input_name, Cpfcan_quality_nCpf[p], description, nanoaod::FlatTable::IntColumn, 10);
      input_name = "CustomTagger_Cpfcan_chi2_" + s;
      description = "chi2 of the charged track fit for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_chi2_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      
      input_name = "CustomTagger_Cpfcan_BtagPf_trackDeltaR_" + s;
      description = "track pseudoangular distance from the jet axis for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackDeltaR_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_BtagPf_trackEtaRel_" + s;
      description = "track pseudorapidity, relative to the jet axis for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackEtaRel_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_BtagPf_trackJetDistVal_" + s;
      description = "minimum track approach distance to jet axis for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackJetDistVal_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_BtagPf_trackPPar_" + s;
      description = "dot product of the jet and track momentum for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPPar_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_BtagPf_trackPParRatio_" + s;
      description = "dot product of the jet and track momentum divided by the magnitude of the jet momentum for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPParRatio_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_BtagPf_trackPtRel_" + s;
      description = "track transverse momentum, relative to the jet axis for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPtRel_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_BtagPf_trackSip2dSig_" + s;
      description = "track 2D signed impact parameter significance for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip2dSig_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_BtagPf_trackSip3dSig_" + s;
      description = "track 3D signed impact parameter significance for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip3dSig_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_BtagPf_trackSip2dVal_" + s;
      description = "track 2D signed impact parameter for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip2dVal_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Cpfcan_BtagPf_trackSip3dVal_" + s;
      description = "track 3D signed impact parameter for the " + s + ". cpf";
      customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip3dVal_nCpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      
  }
      
  // ============================================================== Npfs ===================================================================
  for (unsigned int p = 0; p < n_npf_; p++) {
      auto s = std::to_string(p);
      input_name = "CustomTagger_Npfcan_puppiw_" + s;
      description = "neutral candidate PUPPI weight for the " + s + ". npf";
      customTable->addColumn<float>(input_name, Npfcan_puppiw_nNpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Npfcan_deltaR_" + s;
      description = "pseudoangular distance between the neutral candidate and the jet axis for the " + s + ". npf";
      customTable->addColumn<float>(input_name, Npfcan_deltaR_nNpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Npfcan_drminsv_" + s;
      description = "pseudoangular distance between the neutral candidate and the closest secondary vertex for the " + s + ". npf";
      customTable->addColumn<float>(input_name, Npfcan_drminsv_nNpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Npfcan_HadFrac_" + s;
      description = "fraction of the neutral candidate energy deposited in the hadronic calorimeter for the " + s + ". npf";
      customTable->addColumn<float>(input_name, Npfcan_HadFrac_nNpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Npfcan_ptrel_" + s;
      description = "fraction of the jet momentum carried by the neutral candidate for the " + s + ". npf";
      customTable->addColumn<float>(input_name, Npfcan_ptrel_nNpf[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_Npfcan_isGamma_" + s;
      description = "integer flag indicating whether the neutral candidate is a photon for the " + s + ". npf";
      customTable->addColumn<int>(input_name, Npfcan_isGamma_nNpf[p], description, nanoaod::FlatTable::IntColumn, 10);
      
      
      if (add_CustomTagger_noclip_) {
          // cpf
          input_name = "CustomTagger_Cpfcan_drminsv_noclip_" + s;
          description = "track pseudoangular distance from the closest secondary vertex of the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_drminsv_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_ptrel_noclip_" + s;
          description = "fraction of the jet momentum carried by the track for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_ptrel_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_chi2_noclip_" + s;
          description = "chi2 of the charged track fit for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_chi2_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          
          // cpf from track_builder
          input_name = "CustomTagger_Cpfcan_BtagPf_trackDeltaR_noclip_" + s;
          description = "track pseudoangular distance from the jet axis for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackDeltaR_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_BtagPf_trackEtaRel_noclip_" + s;
          description = "track pseudorapidity, relative to the jet axis for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackEtaRel_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_BtagPf_trackJetDistVal_noclip_" + s;
          description = "minimum track approach distance to jet axis for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackJetDistVal_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_BtagPf_trackPPar_noclip_" + s;
          description = "dot product of the jet and track momentum for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPPar_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_BtagPf_trackPParRatio_noclip_" + s;
          description = "dot product of the jet and track momentum divided by the magnitude of the jet momentum for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPParRatio_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_BtagPf_trackPtRel_noclip_" + s;
          description = "track transverse momentum, relative to the jet axis for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackPtRel_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_BtagPf_trackSip2dSig_noclip_" + s;
          description = "track 2D signed impact parameter significance for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip2dSig_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_BtagPf_trackSip3dSig_noclip_" + s;
          description = "track 3D signed impact parameter significance for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip3dSig_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_BtagPf_trackSip2dVal_noclip_" + s;
          description = "track 2D signed impact parameter for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip2dVal_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Cpfcan_BtagPf_trackSip3dVal_noclip_" + s;
          description = "track 3D signed impact parameter for the " + s + ". cpf (no clip)";
          customTable->addColumn<float>(input_name, Cpfcan_BtagPf_trackSip3dVal_nCpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          
          // npf
          input_name = "CustomTagger_Npfcan_deltaR_noclip_" + s;
          description = "pseudoangular distance between the neutral candidate and the jet axis for the " + s + ". npf (no clip)";
          customTable->addColumn<float>(input_name, Npfcan_deltaR_nNpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Npfcan_drminsv_noclip_" + s;
          description = "pseudoangular distance between the neutral candidate and the closest secondary vertex for the " + s + ". npf (no clip)";
          customTable->addColumn<float>(input_name, Npfcan_drminsv_nNpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_Npfcan_ptrel_noclip_" + s;
          description = "fraction of the jet momentum carried by the neutral candidate for the " + s + ". npf (no clip)";
          customTable->addColumn<float>(input_name, Npfcan_ptrel_nNpf_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
      }
      
  }  
    
  // ============================================================== SVs ===================================================================
  for (unsigned int p = 0; p < n_sv_; p++) {
      auto s = std::to_string(p);
      
      input_name = "CustomTagger_sv_mass_" + s;
      description = "SV mass of the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_mass_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_pt_" + s;
      description = "SV pt of the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_pt_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_ntracks_" + s;
      description = "Number of tracks asociated to the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_ntracks_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_chi2_" + s;
      description = "chi2 of the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_chi2_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_normchi2_" + s;
      description = "chi2/dof of the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_normchi2_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_dxy_" + s;
      description = "2D impact parameter (flight distance) value of the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_dxy_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_dxysig_" + s;
      description = "2D impact parameter (flight distance) significance of the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_dxysig_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_d3d_" + s;
      description = "3D impact parameter (flight distance) value of the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_d3d_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_d3dsig_" + s;
      description = "3D impact parameter (flight distance) significance of the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_d3dsig_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_costhetasvpv_" + s;
      description = "cosine of the angle between the " + s + ". SV flight direction and the direction of the " + s + ". SV momentum";
      customTable->addColumn<float>(input_name, sv_costhetasvpv_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTaggerExtra_sv_phirel_" + s;
      description = "DeltaPhi(sv, jet) for the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_phirel_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTaggerExtra_sv_ptrel_" + s;
      description = "pT relative to parent jet for the " + s + ". SV";
      customTable->addColumn<float>(input_name, sv_ptrel_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_deltaR_" + s;
      description = "pseudoangular distance between jet axis and the " + s + ". SV direction";
      customTable->addColumn<float>(input_name, sv_deltaR_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      input_name = "CustomTagger_sv_enratio_" + s;
      description = "ratio of the " + s + ". SV energy ratio to the jet energy";
      customTable->addColumn<float>(input_name, sv_enratio_nSV[p], description, nanoaod::FlatTable::FloatColumn, 10);
      
      if (add_CustomTagger_noclip_) {
          input_name = "CustomTagger_sv_normchi2_noclip_" + s;
          description = "chi2/dof of the " + s + ". SV (no clip)";
          customTable->addColumn<float>(input_name, sv_normchi2_nSV_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_sv_dxysig_noclip_" + s;
          description = "2D impact parameter (flight distance) significance of the " + s + ". SV (no clip)";
          customTable->addColumn<float>(input_name, sv_dxysig_nSV_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_sv_d3dsig_noclip_" + s;
          description = "3D impact parameter (flight distance) significance of the " + s + ". SV (no clip)";
          customTable->addColumn<float>(input_name, sv_d3dsig_nSV_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
          input_name = "CustomTagger_sv_deltaR_noclip_" + s;
          description = "pseudoangular distance between jet axis and the " + s + ". SV direction (no clip)";
          customTable->addColumn<float>(input_name, sv_deltaR_nSV_noclip[p], description, nanoaod::FlatTable::FloatColumn, 10);
      }
  }
  
    
  iEvent.put(std::move(customTable), nameCustomTagger_);  
  
}

template< typename T>
void JetConstituentTableProducerCustomTagger<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("nameCustomTagger", "Jet");
  desc.add<std::string>("idx_nameCustomTagger", "customIdx");
  desc.add<std::string>("storeAK4Truth","no");
  desc.add<double>("jet_radius", 0.4);
  desc.add<bool>("readBtag", true);
  desc.add<edm::InputTag>("genparticles", edm::InputTag("prunedGenParticles"));
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJets"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("slimmedSecondaryVertices"));
  desc.add<bool>("add_CustomTagger_noclip", false);
  descriptions.addWithDefaultLabel(desc);
}

typedef JetConstituentTableProducerCustomTagger<pat::Jet> PatJetConstituentTableProducerCustomTagger;
//typedef JetConstituentTableProducerCustomTagger<reco::GenJet> GenJetConstituentTableProducerCustomTagger;

DEFINE_FWK_MODULE(PatJetConstituentTableProducerCustomTagger);
//DEFINE_FWK_MODULE(GenJetConstituentTableProducerCustomTagger);