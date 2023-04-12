import FWCore.ParameterSet.Config as cms
# from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.NanoAOD.jets_cff import jetTable, fatJetTable, subJetTable
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from PhysicsTools.PatAlgos.tools.helpers import addToProcessAndTask, getPatAlgosToolsTask


# def update_jets_AK4(process):
#     # Based on ``nanoAOD_addDeepInfo``
#     # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
#     # DeepJet flav_names as found in
#     # https://github.com/cms-sw/cmssw/blob/master/RecoBTag/ONNXRuntime/plugins/DeepFlavourONNXJetTagsProducer.cc#L86
#     # and https://twiki.cern.ch/twiki/bin/view/CMS/DeepJet
#     _btagDiscriminators = [
#         'pfJetProbabilityBJetTags',
#        # 'pfDeepCSVJetTags:probb',
#        # 'pfDeepCSVJetTags:probc',
#        # 'pfDeepCSVJetTags:probbb',
#        # 'pfDeepCSVJetTags:probudsg',       
#        # 'pfDeepFlavourJetTags:probb',
#        # 'pfDeepFlavourJetTags:probbb',
#        # 'pfDeepFlavourJetTags:problepb',
#        # 'pfDeepFlavourJetTags:probc',
#        # 'pfDeepFlavourJetTags:probuds',
#        # 'pfDeepFlavourJetTags:probg'
#     ]
#     
#     updateJetCollection(
#         process,
#         jetSource=cms.InputTag('slimmedJets'),
#         jetCorrections=('AK4PFchs',
#                         cms.vstring(
#                             ['L1FastJet', 'L2Relative', 'L3Absolute',
#                              'L2L3Residual']), 'None'),
#         btagDiscriminators=_btagDiscriminators,
#         postfix='WithDeepInfo',
#     )
#     process.load("Configuration.StandardSequences.MagneticField_cff")
#     process.jetCorrFactorsNano.src = "selectedUpdatedPatJetsWithDeepInfo"
#     process.updatedJets.jetSource = "selectedUpdatedPatJetsWithDeepInfo"
# 
#    # process.updatedPatJetsTransientCorrectedWithDeepInfo.tagInfoSources.append(cms.InputTag("pfDeepCSVTagInfosWithDeepInfo"))
#    # process.updatedPatJetsTransientCorrectedWithDeepInfo.tagInfoSources.append(cms.InputTag("pfDeepFlavourTagInfosWithDeepInfo"))
#     
#     process.updatedPatJetsTransientCorrectedWithDeepInfo.addTagInfos = cms.bool(True)
#     
#     return process


# def update_jets_AK8(process):
#     # Based on ``nanoAOD_addDeepInfoAK8``
#     # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
#     # Care needs to be taken to make sure no discriminators from stock Nano are excluded -> would results in unfilled vars
#     _btagDiscriminators = [
#         'pfJetProbabilityBJetTags',
#         'pfDeepCSVJetTags:probb',
#         'pfDeepCSVJetTags:probc',
#         'pfDeepCSVJetTags:probbb',
#         'pfDeepCSVJetTags:probudsg',
#         'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
#         'pfMassIndependentDeepDoubleCvLJetTags:probHcc',
#         'pfMassIndependentDeepDoubleCvBJetTags:probHcc',
#         'pfMassIndependentDeepDoubleBvLV2JetTags:probHbb',
#         'pfMassIndependentDeepDoubleCvLV2JetTags:probHcc',
#         'pfMassIndependentDeepDoubleCvBV2JetTags:probHcc',
#         ]
#     from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
#     _btagDiscriminators += pfParticleNetJetTagsAll
#     updateJetCollection(
#         process,
#         jetSource=cms.InputTag('slimmedJetsAK8'),
#         pvSource=cms.InputTag('offlineSlimmedPrimaryVertices'),
#         svSource=cms.InputTag('slimmedSecondaryVertices'),
#         rParam=0.8,
#         jetCorrections=('AK8PFPuppi',
#                         cms.vstring([
#                             'L1FastJet', 'L2Relative', 'L3Absolute',
#                             'L2L3Residual'
#                         ]), 'None'),
#         btagDiscriminators=_btagDiscriminators,
#         postfix='AK8WithDeepInfo',
#         # this should work but doesn't seem to enable the tag info with addTagInfos
#         # btagInfos=['pfDeepDoubleXTagInfos'],
#         printWarning=False)
#     process.jetCorrFactorsAK8.src = "selectedUpdatedPatJetsAK8WithDeepInfo"
#     process.updatedJetsAK8.jetSource = "selectedUpdatedPatJetsAK8WithDeepInfo"
#     # add DeepDoubleX taginfos
#     process.updatedPatJetsTransientCorrectedAK8WithDeepInfo.tagInfoSources.append(cms.InputTag("pfDeepDoubleXTagInfosAK8WithDeepInfo"))
#     process.updatedPatJetsTransientCorrectedAK8WithDeepInfo.addTagInfos = cms.bool(True)
#     return process


# def update_jets_AK8_subjet(process):
#     # Based on ``nanoAOD_addDeepInfoAK8``
#     # in https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/nano_cff.py
#     # and https://github.com/alefisico/RecoBTag-PerformanceMeasurements/blob/10_2_X_boostedCommissioning/test/runBTagAnalyzer_cfg.py
#     _btagDiscriminators = [
#         'pfJetProbabilityBJetTags',
#         'pfDeepCSVJetTags:probb',
#         'pfDeepCSVJetTags:probc',
#         'pfDeepCSVJetTags:probbb',
#         'pfDeepCSVJetTags:probudsg',
#         ]
#     updateJetCollection(
#         process,
#         labelName='SoftDropSubjetsPF',
#         jetSource=cms.InputTag("slimmedJetsAK8PFPuppiSoftDropPacked", "SubJets"),
#         jetCorrections=('AK4PFPuppi',
#                         ['L2Relative', 'L3Absolute'], 'None'),
#         btagDiscriminators=list(_btagDiscriminators),
#         explicitJTA=True,  # needed for subjet b tagging
#         svClustering=False,  # needed for subjet b tagging (IMPORTANT: Needs to be set to False to disable ghost-association which does not work with slimmed jets)
#         fatJets=cms.InputTag('slimmedJetsAK8'),  # needed for subjet b tagging
#         rParam=0.8,  # needed for subjet b tagging
#         sortByPt=False, # Don't change order (would mess with subJetIdx for FatJets)
#         postfix='AK8SubjetsWithDeepInfo')
# 
#     process.subJetTable.src = 'selectedUpdatedPatJetsSoftDropSubjetsPFAK8SubjetsWithDeepInfo' 
#     
# 
#     return process



def add_CustomTagger(process, runOnMC=False, keepInputs=['CustomTagger'], storeAK4Truth='no'):
    # addAK4 = not onlyAK8
    # addAK8 = not onlyAK4

    # if addAK4:
    #     process = update_jets_AK4(process)
    # if addAK8:
    #     process = update_jets_AK8(process)
    #     process = update_jets_AK8_subjet(process)

    process.customizeJetTask = cms.Task()
    process.schedule.associate(process.customizeJetTask)

    # CommonVars = cms.PSet(
    #     Proba=Var("bDiscriminator('pfJetProbabilityBJetTags')",
    #               float,
    #               doc="Jet Probability (Usage:BTV)",
    #               precision=10),
    #   #  btagDeepB_b=Var("bDiscriminator('pfDeepCSVJetTags:probb')",
    #   #                  float,
    #   #                  doc="DeepCSV b tag discriminator",
    #   #                  precision=10),
    #   #  btagDeepB_bb=Var("bDiscriminator('pfDeepCSVJetTags:probbb')",
    #   #                   float,
    #   #                   doc="DeepCSV bb tag discriminator",
    #   #                   precision=10),
    #   #  #btagDeepC=Var("bDiscriminator('pfDeepCSVJetTags:probc')",
    #   #  #              float,
    #   #  #              doc="DeepCSV c btag discriminator",
    #   #  #              precision=10), # only necessary after 106X
    #   #  btagDeepL=Var("bDiscriminator('pfDeepCSVJetTags:probudsg')",
    #   #                float,
    #   #                doc="DeepCSV light btag discriminator",
    #   #                precision=10),
    # )
    
    # decouple these from CommonVars, not relevant for data
    HadronCountingVars = cms.PSet(
        nBHadrons=Var("jetFlavourInfo().getbHadrons().size()",
                      int,
                      doc="number of b-hadrons"),
        nCHadrons=Var("jetFlavourInfo().getcHadrons().size()",
                      int,
                      doc="number of c-hadrons")
    )
    

    # AK4
    process.customJetExtTable = cms.EDProducer(
        "SimpleCandidateFlatTableProducer",
        src=jetTable.src,
        cut=jetTable.cut,
        name=jetTable.name,
        doc=jetTable.doc,
        singleton=cms.bool(False),  # the number of entries is variable
        extension=cms.bool(True),  # this is the extension table for Jets
        variables=cms.PSet(
          #   CommonVars,
            HadronCountingVars if runOnMC else cms.PSet(), # hadrons from Generator only relevant for MC
          #  get_DeepCSV_vars() if ('DeepCSV' in keepInputs) else cms.PSet(),
          #  get_DeepJet_outputs()  # outputs are added in any case, inputs only if requested
        ))
    
    if ('CustomTagger' in keepInputs):
        if runOnMC == False and storeAK4Truth == "yes":
            storeAK4Truth = "no" # data does not have truth information, avoid crashes in producer.
        process.customAK4ConstituentsForCustomTaggerTable = cms.EDProducer("PatJetConstituentTableProducerCustomTagger",
                                                                      jets = cms.InputTag("finalJets"),
                                                                      storeAK4Truth = cms.string(storeAK4Truth)
                                                                      )
    
    
  #  # AK8
  #  process.customFatJetExtTable = cms.EDProducer(
  #      "SimpleCandidateFlatTableProducer",
  #      src=fatJetTable.src,
  #      cut=fatJetTable.cut,
  #      name=fatJetTable.name,
  #      doc=fatJetTable.doc,
  #      singleton=cms.bool(False),  # the number of entries is variable
  #      extension=cms.bool(True),  # this is the extension table for FatJets
  #      variables=cms.PSet(
  #          CommonVars,
  #          #HadronCountingVars if runOnMC else cms.PSet(), # only necessary before 106x
  #          #cms.PSet(
  #          #    btagDDBvLV2 = Var("bDiscriminator('pfMassIndependentDeepDoubleBvLV2JetTags:probHbb')",float,doc="DeepDoubleX V2 discriminator for H(Z)->bb vs QCD",precision=10),
  #          #    btagDDCvLV2 = Var("bDiscriminator('pfMassIndependentDeepDoubleCvLV2JetTags:probHcc')",float,doc="DeepDoubleX V2 discriminator for H(Z)->cc vs QCD",precision=10),
  #          #    btagDDCvBV2 = Var("bDiscriminator('pfMassIndependentDeepDoubleCvBV2JetTags:probHcc')",float,doc="DeepDoubleX V2 discriminator for H(Z)->cc vs H(Z)->bb",precision=10),
  #          #), # only necessary before 10_6_19
  #          get_DDX_vars() if ('DDX' in keepInputs) else cms.PSet(),
  #          btagDeepC = Var("bDiscriminator('pfDeepCSVJetTags:probc')",
  #                      float,
  #                      doc="DeepCSV charm btag discriminator",
  #                      precision=10),
  #      ))

  #  # Subjets
  #  process.customSubJetExtTable = cms.EDProducer(
  #      "SimpleCandidateFlatTableProducer",
  #      src=subJetTable.src,
  #      cut=subJetTable.cut,
  #      name=subJetTable.name,
  #      doc=subJetTable.doc,
  #      singleton=cms.bool(False),  # the number of entries is variable
  #      extension=cms.bool(True),  # this is the extension table for FatJets
  #      variables=cms.PSet(
  #          CommonVars,
  #          #HadronCountingVars if runOnMC else cms.PSet(), # only necessary before 106x
  #          btagDeepC = Var("bDiscriminator('pfDeepCSVJetTags:probc')",
  #                      float,
  #                      doc="DeepCSV charm btag discriminator",
  #                      precision=10),

  #  ))

  #  process.customSubJetMCExtTable = cms.EDProducer(
  #      "SimpleCandidateFlatTableProducer",
  #      src = subJetTable.src,
  #      cut = subJetTable.cut,
  #      name = subJetTable.name,
  #      doc=subJetTable.doc,
  #      singleton = cms.bool(False),
  #      extension = cms.bool(True),
  #      variables = cms.PSet(
  #      subGenJetAK8Idx = Var("?genJetFwdRef().backRef().isNonnull()?genJetFwdRef().backRef().key():-1",
  #      int,
  #      doc="index of matched gen Sub jet"),
  #     )
  #  )

   #  if addAK4:
    process.customizeJetTask.add(process.customJetExtTable)
    if ('CustomTagger' in keepInputs):
        process.customizeJetTask.add(process.customAK4ConstituentsForCustomTaggerTable)

  #  if addAK8:
  #      process.customizeJetTask.add(process.customFatJetExtTable)
  #      process.customizeJetTask.add(process.customSubJetExtTable)
  #      if runOnMC:
  #          process.customizeJetTask.add(process.customSubJetMCExtTable)

    return process
