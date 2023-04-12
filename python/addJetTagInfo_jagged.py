import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import Var
from PhysicsTools.NanoAOD.jets_cff import jetTable, fatJetTable, subJetTable
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
from PhysicsTools.PatAlgos.tools.helpers import addToProcessAndTask, getPatAlgosToolsTask

def add_CustomTagger(process, runOnMC=False, keepInputs=['CustomTagger'], storeAK4Truth='no'):

    process.customizeJetTask = cms.Task()
    process.schedule.associate(process.customizeJetTask)

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
            HadronCountingVars if runOnMC else cms.PSet(), # hadrons from Generator only relevant for MC
        ))
    
    if ('CustomTagger' in keepInputs):
        if runOnMC == False and storeAK4Truth == "yes":
            storeAK4Truth = "no" # data does not have truth information, avoid crashes in producer.
        process.customAK4ConstituentsForCustomTaggerTable = cms.EDProducer("PatJetConstituentTableProducerCustomTaggerJagged",
                                                                      jets = cms.InputTag("finalJets"),
                                                                      storeAK4Truth = cms.string(storeAK4Truth)
                                                                      )

   #  if addAK4:
    process.customizeJetTask.add(process.customJetExtTable)
    if ('CustomTagger' in keepInputs):
        process.customizeJetTask.add(process.customAK4ConstituentsForCustomTaggerTable)

    return process
