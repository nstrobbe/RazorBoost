#ifndef RZRBTANALYZER_H
#define RZRBTANALYZER_H
//-----------------------------------------------------------------------------
// File:        rzrBTanalyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Sun Aug  4 13:32:09 2013 by mkanalyzer.py
// Author:      Sezen Sekmen
//-----------------------------------------------------------------------------
// -- System

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>

#ifdef PROJECT_NAME
#include "PhysicsTools/TheNtupleMaker/interface/treestream.h"
#include "PhysicsTools/TheNtupleMaker/interface/pdg.h"
#else
#include "treestream.h"
#include "pdg.h"
#endif

// -- Root

#include "TROOT.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1D.h"
#include "TH2D.h"
//-----------------------------------------------------------------------------
// -- Declare variables to be read
//-----------------------------------------------------------------------------
std::vector<double>	calomet1_energy(200,0);
std::vector<double>	calomet1_mEtSig(200,0);
std::vector<double>	calomet1_phi(200,0);
std::vector<double>	calomet1_pt(200,0);
std::vector<double>	calomet1_significance(200,0);
std::vector<double>	calomet1_sumEt(200,0);
std::vector<double>	calomet_energy(200,0);
std::vector<double>	calomet_mEtSig(200,0);
std::vector<double>	calomet_phi(200,0);
std::vector<double>	calomet_pt(200,0);
std::vector<double>	calomet_significance(200,0);
std::vector<double>	calomet_sumEt(200,0);
std::vector<double>	cmgbasemet1_energy(200,0);
std::vector<double>	cmgbasemet1_et(200,0);
std::vector<double>	cmgbasemet1_phi(200,0);
std::vector<double>	cmgbasemet1_pt(200,0);
std::vector<double>	cmgbasemet1_sumEt(200,0);
std::vector<double>	cmgbasemet2_energy(200,0);
std::vector<double>	cmgbasemet2_et(200,0);
std::vector<double>	cmgbasemet2_phi(200,0);
std::vector<double>	cmgbasemet2_pt(200,0);
std::vector<double>	cmgbasemet2_sumEt(200,0);
std::vector<double>	cmgbasemet_energy(200,0);
std::vector<double>	cmgbasemet_et(200,0);
std::vector<double>	cmgbasemet_phi(200,0);
std::vector<double>	cmgbasemet_pt(200,0);
std::vector<double>	cmgbasemet_sumEt(200,0);
std::vector<int>	cmgelectron1_charge(200,0);
std::vector<double>	cmgelectron1_energy(200,0);
std::vector<double>	cmgelectron1_et(200,0);
std::vector<double>	cmgelectron1_eta(200,0);
std::vector<double>	cmgelectron1_phi(200,0);
std::vector<double>	cmgelectron1_pt(200,0);
std::vector<int>	cmgelectron_charge(200,0);
std::vector<double>	cmgelectron_energy(200,0);
std::vector<double>	cmgelectron_et(200,0);
std::vector<double>	cmgelectron_eta(200,0);
std::vector<double>	cmgelectron_phi(200,0);
std::vector<double>	cmgelectron_pt(200,0);
std::vector<int>	cmgmuon1_charge(200,0);
std::vector<double>	cmgmuon1_energy(200,0);
std::vector<double>	cmgmuon1_et(200,0);
std::vector<double>	cmgmuon1_eta(200,0);
std::vector<double>	cmgmuon1_phi(200,0);
std::vector<double>	cmgmuon1_pt(200,0);
std::vector<int>	cmgmuon_charge(200,0);
std::vector<double>	cmgmuon_energy(200,0);
std::vector<double>	cmgmuon_et(200,0);
std::vector<double>	cmgmuon_eta(200,0);
std::vector<double>	cmgmuon_phi(200,0);
std::vector<double>	cmgmuon_pt(200,0);
std::vector<double>	cmgpfjet_combinedSecondaryVertexBJetTags(200,0);
std::vector<double>	cmgpfjet_combinedSecondaryVertexMVABJetTags(200,0);
std::vector<double>	cmgpfjet_component_0_fraction(200,0);
std::vector<double>	cmgpfjet_component_0_number(200,0);
std::vector<double>	cmgpfjet_component_1_fraction(200,0);
std::vector<double>	cmgpfjet_component_1_number(200,0);
std::vector<double>	cmgpfjet_component_2_fraction(200,0);
std::vector<double>	cmgpfjet_component_2_number(200,0);
std::vector<double>	cmgpfjet_component_3_fraction(200,0);
std::vector<double>	cmgpfjet_component_3_number(200,0);
std::vector<double>	cmgpfjet_component_4_fraction(200,0);
std::vector<double>	cmgpfjet_component_4_number(200,0);
std::vector<double>	cmgpfjet_component_5_fraction(200,0);
std::vector<double>	cmgpfjet_component_5_number(200,0);
std::vector<double>	cmgpfjet_component_6_fraction(200,0);
std::vector<double>	cmgpfjet_component_6_number(200,0);
std::vector<double>	cmgpfjet_component_7_fraction(200,0);
std::vector<double>	cmgpfjet_component_7_number(200,0);
std::vector<double>	cmgpfjet_energy(200,0);
std::vector<double>	cmgpfjet_et(200,0);
std::vector<double>	cmgpfjet_eta(200,0);
std::vector<float>	cmgpfjet_jetArea(200,0);
std::vector<double>	cmgpfjet_jetBProbabilityBJetTags(200,0);
std::vector<double>	cmgpfjet_jetProbabilityBJetTags(200,0);
std::vector<double>	cmgpfjet_mass(200,0);
std::vector<int>	cmgpfjet_nConstituents(200,0);
std::vector<int>	cmgpfjet_partonFlavour(200,0);
std::vector<double>	cmgpfjet_phi(200,0);
std::vector<double>	cmgpfjet_pt(200,0);
std::vector<double>	cmgpfjet_rapidity(200,0);
std::vector<double>	cmgpfjet_trackCountingHighEffBJetTag(200,0);
std::vector<double>	cmgpfjet_trackCountingHighPurBJetTags(200,0);
std::vector<float>	cmgtau1_byTightCombinedIsolationDeltaBetaCorr(200,0);
std::vector<int>	cmgtau1_charge(200,0);
std::vector<double>	cmgtau1_energy(200,0);
std::vector<double>	cmgtau1_et(200,0);
std::vector<double>	cmgtau1_eta(200,0);
std::vector<double>	cmgtau1_phi(200,0);
std::vector<double>	cmgtau1_pt(200,0);
std::vector<float>	cmgtau_byLooseCombinedIsolationDeltaBetaCorr(200,0);
std::vector<int>	cmgtau_charge(200,0);
std::vector<double>	cmgtau_energy(200,0);
std::vector<double>	cmgtau_et(200,0);
std::vector<double>	cmgtau_eta(200,0);
std::vector<double>	cmgtau_phi(200,0);
std::vector<double>	cmgtau_pt(200,0);
std::vector<float>	electronhelper_caloIso(200,0);
std::vector<int>	electronhelper_charge(200,0);
std::vector<float>	electronhelper_deltaEtaSuperClusterTrackAtVtx(200,0);
std::vector<float>	electronhelper_deltaPhiSuperClusterTrackAtVtx(200,0);
std::vector<double>	electronhelper_dxywrtPV(200,0);
std::vector<double>	electronhelper_dzwrtPV(200,0);
std::vector<float>	electronhelper_eSuperClusterOverP(200,0);
std::vector<float>	electronhelper_ecalIso(200,0);
std::vector<double>	electronhelper_energy(200,0);
std::vector<double>	electronhelper_et(200,0);
std::vector<double>	electronhelper_eta(200,0);
std::vector<unsigned short>	electronhelper_gsfTrack_trackerExpectedHitsInner_numberOfHits(200,0);
std::vector<float>	electronhelper_hadronicOverEm(200,0);
std::vector<float>	electronhelper_hcalIso(200,0);
std::vector<int>	electronhelper_isEB(200,0);
std::vector<int>	electronhelper_isEE(200,0);
std::vector<double>	electronhelper_phi(200,0);
std::vector<double>	electronhelper_pt(200,0);
std::vector<float>	electronhelper_scSigmaIEtaIEta(200,0);
std::vector<float>	electronhelper_simpleEleId80relIso(200,0);
std::vector<float>	electronhelper_simpleEleId95relIso(200,0);
std::vector<double>	electronhelper_superCluster_energy(200,0);
std::vector<float>	electronhelper_trackIso(200,0);
int	eventhelper_bunchCrossing;
int	eventhelper_event;
int	eventhelper_isRealData;
int	eventhelper_luminosityBlock;
int	eventhelper_orbitNumber;
int	eventhelper_run;
int	eventhelperextra_numberOfPrimaryVertices;
int	eventhelperextra_trackIso;
int	eventhelperextra_trackIsoLep;
double	geneventinfoproduct_alphaQCD;
double	geneventinfoproduct_alphaQED;
int	geneventinfoproduct_hasBinningValues;
int	geneventinfoproduct_hasPDF;
double	geneventinfoproduct_qScale;
unsigned int	geneventinfoproduct_signalProcessID;
double	geneventinfoproduct_weight;
std::vector<double>	genjet_energy(200,0);
std::vector<double>	genjet_et(200,0);
std::vector<double>	genjet_eta(200,0);
std::vector<double>	genjet_mass(200,0);
std::vector<int>	genjet_nConstituents(200,0);
std::vector<double>	genjet_phi(200,0);
std::vector<double>	genjet_pt(200,0);
std::vector<double>	genjet_rapidity(200,0);
std::vector<int>	genparticlehelper_charge(200,0);
std::vector<double>	genparticlehelper_eta(200,0);
std::vector<int>	genparticlehelper_firstDaughter(200,0);
std::vector<int>	genparticlehelper_firstMother(200,0);
std::vector<int>	genparticlehelper_lastDaughter(200,0);
std::vector<int>	genparticlehelper_lastMother(200,0);
std::vector<double>	genparticlehelper_mass(200,0);
std::vector<int>	genparticlehelper_pdgId(200,0);
std::vector<double>	genparticlehelper_phi(200,0);
std::vector<double>	genparticlehelper_pt(200,0);
std::vector<int>	genparticlehelper_status(200,0);
int     geneventinfoproducthelper_id1;
int     geneventinfoproducthelper_id2;
double  geneventinfoproducthelper_q;
double  geneventinfoproducthelper_x1;
double  geneventinfoproducthelper_x2;
double	genruninfoproduct_crossSection;
double	genruninfoproduct_filterEfficiency;
std::vector<float>	jethelper1_chargedEmEnergyFraction(200,0);
std::vector<float>	jethelper1_chargedHadronEnergyFraction(200,0);
std::vector<int>	jethelper1_chargedMultiplicity(200,0);
std::vector<float>	jethelper1_combinedSecondaryVertexBJetTags(200,0);
std::vector<float>	jethelper1_combinedSecondaryVertexMVABJetTags(200,0);
std::vector<double>	jethelper1_energy(200,0);
std::vector<double>	jethelper1_et(200,0);
std::vector<double>	jethelper1_eta(200,0);
std::vector<float>	jethelper1_jetArea(200,0);
std::vector<float>	jethelper1_jetBProbabilityBJetTags(200,0);
std::vector<float>	jethelper1_jetCharge03(200,0);
std::vector<float>	jethelper1_jetCharge05(200,0);
std::vector<float>	jethelper1_jetCharge10(200,0);
std::vector<float>	jethelper1_jetProbabilityBJetTags(200,0);
std::vector<double>	jethelper1_mass(200,0);
std::vector<float>	jethelper1_muonEnergyFraction(200,0);
std::vector<int>	jethelper1_nConstituents(200,0);
std::vector<float>	jethelper1_neutralEmEnergyFraction(200,0);
std::vector<float>	jethelper1_neutralHadronEnergyFraction(200,0);
std::vector<int>	jethelper1_partonFlavour(200,0);
std::vector<double>	jethelper1_phi(200,0);
std::vector<float>	jethelper1_photonEnergyFraction(200,0);
std::vector<double>	jethelper1_pt(200,0);
std::vector<double>	jethelper1_rapidity(200,0);
std::vector<float>	jethelper1_trackCountingHighEffBJetTags(200,0);
std::vector<float>	jethelper1_trackCountingHighPurBJetTags(200,0);
std::vector<double>	jethelper1_uncor_energy(200,0);
std::vector<double>	jethelper1_uncor_et(200,0);
std::vector<double>	jethelper1_uncor_pt(200,0);
std::vector<float>	jethelper2_combinedSecondaryVertexBJetTags(200,0);
std::vector<float>	jethelper2_combinedSecondaryVertexMVABJetTags(200,0);
std::vector<double>	jethelper2_daughter_0_energy(200,0);
std::vector<double>	jethelper2_daughter_0_eta(200,0);
std::vector<float>	jethelper2_daughter_0_jetCharge03(200,0);
std::vector<float>	jethelper2_daughter_0_jetCharge05(200,0);
std::vector<float>	jethelper2_daughter_0_jetCharge10(200,0);
std::vector<double>	jethelper2_daughter_0_mass(200,0);
std::vector<double>	jethelper2_daughter_0_phi(200,0);
std::vector<double>	jethelper2_daughter_0_pt(200,0);
std::vector<double>	jethelper2_daughter_0_rapidity(200,0);
std::vector<double>	jethelper2_daughter_1_energy(200,0);
std::vector<double>	jethelper2_daughter_1_eta(200,0);
std::vector<float>	jethelper2_daughter_1_jetCharge03(200,0);
std::vector<float>	jethelper2_daughter_1_jetCharge05(200,0);
std::vector<float>	jethelper2_daughter_1_jetCharge10(200,0);
std::vector<double>	jethelper2_daughter_1_mass(200,0);
std::vector<double>	jethelper2_daughter_1_phi(200,0);
std::vector<double>	jethelper2_daughter_1_pt(200,0);
std::vector<double>	jethelper2_daughter_1_rapidity(200,0);
std::vector<double>	jethelper2_energy(200,0);
std::vector<double>	jethelper2_et(200,0);
std::vector<double>	jethelper2_eta(200,0);
std::vector<float>	jethelper2_jetArea(200,0);
std::vector<float>	jethelper2_jetBProbabilityBJetTags(200,0);
std::vector<float>	jethelper2_jetCharge03(200,0);
std::vector<float>	jethelper2_jetCharge05(200,0);
std::vector<float>	jethelper2_jetCharge10(200,0);
std::vector<float>	jethelper2_jetProbabilityBJetTags(200,0);
std::vector<double>	jethelper2_mass(200,0);
std::vector<int>	jethelper2_nConstituents(200,0);
std::vector<size_t>	jethelper2_numberOfDaughters(200,0);
std::vector<int>	jethelper2_partonFlavour(200,0);
std::vector<double>	jethelper2_phi(200,0);
std::vector<double>	jethelper2_pt(200,0);
std::vector<double>	jethelper2_rapidity(200,0);
std::vector<float>	jethelper2_trackCountingHighEffBJetTags(200,0);
std::vector<float>	jethelper2_trackCountingHighPurBJetTags(200,0);
std::vector<double>	jethelper2_uncor_energy(200,0);
std::vector<double>	jethelper2_uncor_et(200,0);
std::vector<double>	jethelper2_uncor_pt(200,0);
std::vector<float>	jethelper3_C2beta17(200,0);
std::vector<float>	jethelper3_chargedEmEnergyFraction(200,0);
std::vector<float>	jethelper3_chargedHadronEnergyFraction(200,0);
std::vector<int>	jethelper3_chargedMultiplicity(200,0);
std::vector<float>	jethelper3_combinedSecondaryVertexBJetTags(200,0);
std::vector<float>	jethelper3_combinedSecondaryVertexMVABJetTags(200,0);
std::vector<double>	jethelper3_energy(200,0);
std::vector<double>	jethelper3_et(200,0);
std::vector<double>	jethelper3_eta(200,0);
std::vector<float>	jethelper3_jetArea(200,0);
std::vector<float>	jethelper3_jetBProbabilityBJetTags(200,0);
std::vector<float>	jethelper3_jetCharge03(200,0);
std::vector<float>	jethelper3_jetCharge05(200,0);
std::vector<float>	jethelper3_jetCharge10(200,0);
std::vector<float>	jethelper3_jetProbabilityBJetTags(200,0);
std::vector<double>	jethelper3_mass(200,0);
std::vector<float>	jethelper3_muonEnergyFraction(200,0);
std::vector<int>	jethelper3_nConstituents(200,0);
std::vector<float>	jethelper3_neutralEmEnergyFraction(200,0);
std::vector<float>	jethelper3_neutralHadronEnergyFraction(200,0);
std::vector<int>	jethelper3_partonFlavour(200,0);
std::vector<double>	jethelper3_phi(200,0);
std::vector<float>	jethelper3_photonEnergyFraction(200,0);
std::vector<double>	jethelper3_pt(200,0);
std::vector<float>	jethelper3_qjetsvolatility(200,0);
std::vector<double>	jethelper3_rapidity(200,0);
std::vector<float>	jethelper3_tau1(200,0);
std::vector<float>	jethelper3_tau2(200,0);
std::vector<float>	jethelper3_tau3(200,0);
std::vector<float>	jethelper3_trackCountingHighEffBJetTags(200,0);
std::vector<float>	jethelper3_trackCountingHighPurBJetTags(200,0);
std::vector<double>	jethelper3_uncor_energy(200,0);
std::vector<double>	jethelper3_uncor_et(200,0);
std::vector<double>	jethelper3_uncor_pt(200,0);
std::vector<float>	jethelper4_combinedSecondaryVertexBJetTags(200,0);
std::vector<float>	jethelper4_combinedSecondaryVertexMVABJetTags(200,0);
std::vector<double>	jethelper4_daughter_0_energy(200,0);
std::vector<double>	jethelper4_daughter_0_eta(200,0);
std::vector<float>	jethelper4_daughter_0_jetCharge03(200,0);
std::vector<float>	jethelper4_daughter_0_jetCharge05(200,0);
std::vector<float>	jethelper4_daughter_0_jetCharge10(200,0);
std::vector<double>	jethelper4_daughter_0_mass(200,0);
std::vector<double>	jethelper4_daughter_0_phi(200,0);
std::vector<double>	jethelper4_daughter_0_pt(200,0);
std::vector<double>	jethelper4_daughter_0_rapidity(200,0);
std::vector<double>	jethelper4_daughter_1_energy(200,0);
std::vector<double>	jethelper4_daughter_1_eta(200,0);
std::vector<float>	jethelper4_daughter_1_jetCharge03(200,0);
std::vector<float>	jethelper4_daughter_1_jetCharge05(200,0);
std::vector<float>	jethelper4_daughter_1_jetCharge10(200,0);
std::vector<double>	jethelper4_daughter_1_mass(200,0);
std::vector<double>	jethelper4_daughter_1_phi(200,0);
std::vector<double>	jethelper4_daughter_1_pt(200,0);
std::vector<double>	jethelper4_daughter_1_rapidity(200,0);
std::vector<double>	jethelper4_energy(200,0);
std::vector<double>	jethelper4_et(200,0);
std::vector<double>	jethelper4_eta(200,0);
std::vector<float>	jethelper4_jetArea(200,0);
std::vector<float>	jethelper4_jetBProbabilityBJetTags(200,0);
std::vector<float>	jethelper4_jetCharge03(200,0);
std::vector<float>	jethelper4_jetCharge05(200,0);
std::vector<float>	jethelper4_jetCharge10(200,0);
std::vector<float>	jethelper4_jetProbabilityBJetTags(200,0);
std::vector<double>	jethelper4_mass(200,0);
std::vector<int>	jethelper4_nConstituents(200,0);
std::vector<size_t>	jethelper4_numberOfDaughters(200,0);
std::vector<int>	jethelper4_partonFlavour(200,0);
std::vector<double>	jethelper4_phi(200,0);
std::vector<double>	jethelper4_pt(200,0);
std::vector<double>	jethelper4_rapidity(200,0);
std::vector<float>	jethelper4_tau1(200,0);
std::vector<float>	jethelper4_tau2(200,0);
std::vector<float>	jethelper4_tau3(200,0);
std::vector<float>	jethelper4_trackCountingHighEffBJetTags(200,0);
std::vector<float>	jethelper4_trackCountingHighPurBJetTags(200,0);
std::vector<double>	jethelper4_uncor_energy(200,0);
std::vector<double>	jethelper4_uncor_et(200,0);
std::vector<double>	jethelper4_uncor_pt(200,0);
std::vector<float>	jethelper5_C2beta17(200,0);
std::vector<float>	jethelper5_chargedEmEnergyFraction(200,0);
std::vector<float>	jethelper5_chargedHadronEnergyFraction(200,0);
std::vector<int>	jethelper5_chargedMultiplicity(200,0);
std::vector<float>	jethelper5_combinedSecondaryVertexBJetTags(200,0);
std::vector<float>	jethelper5_combinedSecondaryVertexMVABJetTags(200,0);
std::vector<double>	jethelper5_energy(200,0);
std::vector<double>	jethelper5_et(200,0);
std::vector<double>	jethelper5_eta(200,0);
std::vector<int>	jethelper5_getChargedPt0(200,0);
std::vector<int>	jethelper5_getChargedPt1(200,0);
std::vector<int>	jethelper5_getChargedPt2(200,0);
std::vector<int>	jethelper5_getChargedPt3(200,0);
std::vector<int>	jethelper5_getNcharged01(200,0);
std::vector<int>	jethelper5_getNneutral01(200,0);
std::vector<int>	jethelper5_getPt0(200,0);
std::vector<int>	jethelper5_getPt1(200,0);
std::vector<int>	jethelper5_getPt2(200,0);
std::vector<int>	jethelper5_getPt3(200,0);
std::vector<float>	jethelper5_jetArea(200,0);
std::vector<float>	jethelper5_jetBProbabilityBJetTags(200,0);
std::vector<float>	jethelper5_jetCharge03(200,0);
std::vector<float>	jethelper5_jetCharge05(200,0);
std::vector<float>	jethelper5_jetCharge10(200,0);
std::vector<float>	jethelper5_jetProbabilityBJetTags(200,0);
std::vector<double>	jethelper5_mass(200,0);
std::vector<float>	jethelper5_muonEnergyFraction(200,0);
std::vector<int>	jethelper5_nConstituents(200,0);
std::vector<float>	jethelper5_neutralEmEnergyFraction(200,0);
std::vector<float>	jethelper5_neutralHadronEnergyFraction(200,0);
std::vector<int>	jethelper5_partonFlavour(200,0);
std::vector<double>	jethelper5_phi(200,0);
std::vector<float>	jethelper5_photonEnergyFraction(200,0);
std::vector<double>	jethelper5_pt(200,0);
std::vector<float>	jethelper5_qjetsvolatility(200,0);
std::vector<double>	jethelper5_rapidity(200,0);
std::vector<float>	jethelper5_tau1(200,0);
std::vector<float>	jethelper5_tau2(200,0);
std::vector<float>	jethelper5_tau3(200,0);
std::vector<float>	jethelper5_trackCountingHighEffBJetTags(200,0);
std::vector<float>	jethelper5_trackCountingHighPurBJetTags(200,0);
std::vector<double>	jethelper5_uncor_energy(200,0);
std::vector<double>	jethelper5_uncor_et(200,0);
std::vector<double>	jethelper5_uncor_pt(200,0);
std::vector<double>	jethelper6_energy(200,0);
std::vector<double>	jethelper6_et(200,0);
std::vector<double>	jethelper6_eta(200,0);
std::vector<float>	jethelper6_genC2beta17(200,0);
std::vector<float>	jethelper6_genC2beta17CHS(200,0);
std::vector<float>	jethelper6_genJetCharge03(200,0);
std::vector<float>	jethelper6_genJetCharge05(200,0);
std::vector<float>	jethelper6_genJetCharge10(200,0);
std::vector<float>	jethelper6_genNCHS(200,0);
std::vector<float>	jethelper6_genTau1(200,0);
std::vector<float>	jethelper6_genTau1CHS(200,0);
std::vector<float>	jethelper6_genTau1Pt2(200,0);
std::vector<float>	jethelper6_genTau1Pt5(200,0);
std::vector<float>	jethelper6_genTau2(200,0);
std::vector<float>	jethelper6_genTau2CHS(200,0);
std::vector<float>	jethelper6_genTau2Pt2(200,0);
std::vector<float>	jethelper6_genTau2Pt5(200,0);
std::vector<float>	jethelper6_genTau3(200,0);
std::vector<float>	jethelper6_jetArea(200,0);
std::vector<double>	jethelper6_mass(200,0);
std::vector<int>	jethelper6_nConstituents(200,0);
std::vector<int>	jethelper6_partonFlavour(200,0);
std::vector<double>	jethelper6_phi(200,0);
std::vector<double>	jethelper6_pt(200,0);
std::vector<double>	jethelper6_rapidity(200,0);
std::vector<double>	jethelper7_daughter_0_energy(200,0);
std::vector<double>	jethelper7_daughter_0_eta(200,0);
std::vector<double>	jethelper7_daughter_0_mass(200,0);
std::vector<double>	jethelper7_daughter_0_phi(200,0);
std::vector<double>	jethelper7_daughter_0_pt(200,0);
std::vector<double>	jethelper7_daughter_0_rapidity(200,0);
std::vector<double>	jethelper7_daughter_1_energy(200,0);
std::vector<double>	jethelper7_daughter_1_eta(200,0);
std::vector<double>	jethelper7_daughter_1_mass(200,0);
std::vector<double>	jethelper7_daughter_1_phi(200,0);
std::vector<double>	jethelper7_daughter_1_pt(200,0);
std::vector<double>	jethelper7_daughter_1_rapidity(200,0);
std::vector<double>	jethelper7_energy(200,0);
std::vector<double>	jethelper7_et(200,0);
std::vector<double>	jethelper7_eta(200,0);
std::vector<float>	jethelper7_genTau1(200,0);
std::vector<float>	jethelper7_genTau2(200,0);
std::vector<float>	jethelper7_genTau3(200,0);
std::vector<float>	jethelper7_jetArea(200,0);
std::vector<float>	jethelper7_jetCharge03(200,0);
std::vector<float>	jethelper7_jetCharge05(200,0);
std::vector<float>	jethelper7_jetCharge10(200,0);
std::vector<double>	jethelper7_mass(200,0);
std::vector<int>	jethelper7_nConstituents(200,0);
std::vector<size_t>	jethelper7_numberOfDaughters(200,0);
std::vector<int>	jethelper7_partonFlavour(200,0);
std::vector<double>	jethelper7_phi(200,0);
std::vector<double>	jethelper7_pt(200,0);
std::vector<double>	jethelper7_rapidity(200,0);
std::vector<float>	jethelper_chargedEmEnergyFraction(200,0);
std::vector<float>	jethelper_chargedHadronEnergyFraction(200,0);
std::vector<int>	jethelper_chargedMultiplicity(200,0);
std::vector<float>	jethelper_combinedSecondaryVertexBJetTags(200,0);
std::vector<float>	jethelper_combinedSecondaryVertexMVABJetTags(200,0);
std::vector<double>	jethelper_energy(200,0);
std::vector<double>	jethelper_et(200,0);
std::vector<double>	jethelper_eta(200,0);
std::vector<float>	jethelper_jetArea(200,0);
std::vector<float>	jethelper_jetBProbabilityBJetTags(200,0);
std::vector<float>	jethelper_jetCharge03(200,0);
std::vector<float>	jethelper_jetCharge05(200,0);
std::vector<float>	jethelper_jetCharge10(200,0);
std::vector<float>	jethelper_jetProbabilityBJetTags(200,0);
std::vector<double>	jethelper_mass(200,0);
std::vector<float>	jethelper_muonEnergyFraction(200,0);
std::vector<int>	jethelper_nConstituents(200,0);
std::vector<float>	jethelper_neutralEmEnergyFraction(200,0);
std::vector<float>	jethelper_neutralHadronEnergyFraction(200,0);
std::vector<int>	jethelper_partonFlavour(200,0);
std::vector<double>	jethelper_phi(200,0);
std::vector<float>	jethelper_photonEnergyFraction(200,0);
std::vector<double>	jethelper_pt(200,0);
std::vector<double>	jethelper_rapidity(200,0);
std::vector<float>	jethelper_trackCountingHighEffBJetTags(200,0);
std::vector<float>	jethelper_trackCountingHighPurBJetTags(200,0);
std::vector<double>	jethelper_uncor_energy(200,0);
std::vector<double>	jethelper_uncor_et(200,0);
std::vector<double>	jethelper_uncor_pt(200,0);
std::vector<int>	leafcandidate_charge(200,0);
std::vector<double>	leafcandidate_energy(200,0);
std::vector<double>	leafcandidate_et(200,0);
std::vector<double>	leafcandidate_eta(200,0);
std::vector<double>	leafcandidate_phi(200,0);
std::vector<double>	leafcandidate_pt(200,0);
double	lheeventproduct_hepeup_AQCDUP;
double	lheeventproduct_hepeup_AQEDUP;
int	lheeventproduct_hepeup_IDPRUP;
int	lheeventproduct_hepeup_NUP;
double	lheeventproduct_hepeup_SCALUP;
double	lheeventproduct_hepeup_XWGTUP;
double	lheeventproducthelper_mg;
double	lheeventproducthelper_mt1;
double	lheeventproducthelper_mz1;
std::vector<double>	met1_energy(200,0);
std::vector<double>	met1_et(200,0);
std::vector<double>	met1_mEtSig(200,0);
std::vector<double>	met1_phi(200,0);
std::vector<double>	met1_pt(200,0);
std::vector<double>	met1_significance(200,0);
std::vector<double>	met1_sumEt(200,0);
std::vector<double>	met_energy(200,0);
std::vector<double>	met_et(200,0);
std::vector<double>	met_mEtSig(200,0);
std::vector<double>	met_phi(200,0);
std::vector<double>	met_pt(200,0);
std::vector<double>	met_significance(200,0);
std::vector<double>	met_sumEt(200,0);
std::vector<float>	muonhelper_TMOneStationTight(200,0);
std::vector<int>	muonhelper_charge(200,0);
std::vector<double>	muonhelper_dB(200,0);
std::vector<double>	muonhelper_dxywrtPV(200,0);
std::vector<double>	muonhelper_dzwrtPV(200,0);
std::vector<double>	muonhelper_energy(200,0);
std::vector<double>	muonhelper_et(200,0);
std::vector<double>	muonhelper_eta(200,0);
std::vector<unsigned short>	muonhelper_globalTrack_hitPattern_numberOfValidMuonHits(200,0);
std::vector<double>	muonhelper_globalTrack_normalizedChi2(200,0);
std::vector<unsigned short>	muonhelper_innerTrack_hitPattern_numberOfValidPixelHits(200,0);
std::vector<unsigned short>	muonhelper_innerTrack_hitPattern_pixelLayersWithMeasurement(200,0);
std::vector<double>	muonhelper_innerTrack_normalizedChi2(200,0);
std::vector<int>	muonhelper_isGlobalMuon(200,0);
std::vector<int>	muonhelper_isPFMuon(200,0);
std::vector<int>	muonhelper_isTrackerMuon(200,0);
std::vector<int>	muonhelper_numberOfMatchedStations(200,0);
std::vector<float>	muonhelper_pfIsolationR04_sumChargedHadronPt(200,0);
std::vector<float>	muonhelper_pfIsolationR04_sumChargedParticlePt(200,0);
std::vector<float>	muonhelper_pfIsolationR04_sumNeutralHadronEt(200,0);
std::vector<float>	muonhelper_pfIsolationR04_sumNeutralHadronEtHighThreshold(200,0);
std::vector<float>	muonhelper_pfIsolationR04_sumPUPt(200,0);
std::vector<float>	muonhelper_pfIsolationR04_sumPhotonEt(200,0);
std::vector<float>	muonhelper_pfIsolationR04_sumPhotonEtHighThreshold(200,0);
std::vector<double>	muonhelper_phi(200,0);
std::vector<double>	muonhelper_pt(200,0);
std::vector<unsigned short>	muonhelper_track_hitPattern_trackerLayersWithMeasurement(200,0);
int	ncalomet;
int	ncalomet1;
int	ncmgbasemet;
int	ncmgbasemet1;
int	ncmgbasemet2;
int	ncmgelectron;
int	ncmgelectron1;
int	ncmgmuon;
int	ncmgmuon1;
int	ncmgpfjet;
int	ncmgtau;
int	ncmgtau1;
int	nelectronhelper;
int	ngenjet;
int	ngenparticlehelper;
int	njethelper;
int	njethelper1;
int	njethelper2;
int	njethelper3;
int	njethelper4;
int	njethelper5;
int	njethelper6;
int	njethelper7;
int	nleafcandidate;
int	nmet;
int	nmet1;
int	nmuonhelper;
int	npileupsummaryinfo;
int	ntau;
int	nvertex;
std::vector<int>	pileupsummaryinfo_getBunchCrossing(200,0);
std::vector<int>	pileupsummaryinfo_getPU_NumInteractions(200,0);
std::vector<float>	pileupsummaryinfo_getTrueNumInteractions(200,0);
double	sdouble_kt6PFJets_rho_value;
double	sdouble_vertexWeightSummer12MC53X2012ABCDData_value;
int	sint_hcallasereventfilter2012_value;
int	sint_simpleGenInfo_value;
std::vector<float>	tau_byLooseCombinedIsolationDeltaBetaCorr(200,0);
std::vector<float>	tau_byMediumCombinedIsolationDeltaBetaCorr(200,0);
std::vector<float>	tau_caloIso(200,0);
std::vector<float>	tau_ecalIso(200,0);
std::vector<double>	tau_energy(200,0);
std::vector<double>	tau_et(200,0);
std::vector<double>	tau_eta(200,0);
std::vector<float>	tau_hcalIso(200,0);
std::vector<double>	tau_phi(200,0);
std::vector<double>	tau_pt(200,0);
std::vector<float>	tau_trackIso(200,0);
int	triggerresultshelper1_CSCTightHaloFilterPath;
int	triggerresultshelper1_EcalDeadCellTriggerPrimitiveFilterPath;
int	triggerresultshelper1_HBHENoiseFilterPath;
int	triggerresultshelper1_eeBadScFilterPath;
int	triggerresultshelper1_hcalLaserEventFilterPath;
int	triggerresultshelper1_metNoiseCleaningPath;
int	triggerresultshelper1_noscrapingFilterPath;
int	triggerresultshelper1_primaryVertexFilterPath;
int	triggerresultshelper1_trackingFailureFilterPath;
int	triggerresultshelper1_trkPOGFiltersPath;
int	triggerresultshelper2_totalKinematicsFilterPath;
int	triggerresultshelper3_trackIsolationMakerFilterPath;
int	triggerresultshelper_HLT_DiPFJetAve320_v10;
int	triggerresultshelper_HLT_DiPFJetAve320_v11;
int	triggerresultshelper_HLT_DiPFJetAve320_v12;
int	triggerresultshelper_HLT_DiPFJetAve320_v2;
int	triggerresultshelper_HLT_DiPFJetAve320_v3;
int	triggerresultshelper_HLT_DiPFJetAve320_v4;
int	triggerresultshelper_HLT_DiPFJetAve320_v5;
int	triggerresultshelper_HLT_DiPFJetAve320_v6;
int	triggerresultshelper_HLT_DiPFJetAve320_v7;
int	triggerresultshelper_HLT_DiPFJetAve320_v8;
int	triggerresultshelper_HLT_DiPFJetAve320_v9;
int	triggerresultshelper_HLT_DiPFJetAve400_v10;
int	triggerresultshelper_HLT_DiPFJetAve400_v11;
int	triggerresultshelper_HLT_DiPFJetAve400_v12;
int	triggerresultshelper_HLT_DiPFJetAve400_v2;
int	triggerresultshelper_HLT_DiPFJetAve400_v3;
int	triggerresultshelper_HLT_DiPFJetAve400_v4;
int	triggerresultshelper_HLT_DiPFJetAve400_v5;
int	triggerresultshelper_HLT_DiPFJetAve400_v6;
int	triggerresultshelper_HLT_DiPFJetAve400_v7;
int	triggerresultshelper_HLT_DiPFJetAve400_v8;
int	triggerresultshelper_HLT_DiPFJetAve400_v9;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v11;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v12;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v6;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v7;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v8;
int	triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v9;
int	triggerresultshelper_HLT_HT450_v1;
int	triggerresultshelper_HLT_HT450_v2;
int	triggerresultshelper_HLT_HT450_v3;
int	triggerresultshelper_HLT_HT450_v4;
int	triggerresultshelper_HLT_HT450_v5;
int	triggerresultshelper_HLT_HT450_v6;
int	triggerresultshelper_HLT_HT450_v7;
int	triggerresultshelper_HLT_HT450_v8;
int	triggerresultshelper_HLT_HT450_v9;
int	triggerresultshelper_HLT_HT500_v1;
int	triggerresultshelper_HLT_HT500_v2;
int	triggerresultshelper_HLT_HT500_v3;
int	triggerresultshelper_HLT_HT500_v4;
int	triggerresultshelper_HLT_HT500_v5;
int	triggerresultshelper_HLT_HT500_v6;
int	triggerresultshelper_HLT_HT500_v7;
int	triggerresultshelper_HLT_HT500_v8;
int	triggerresultshelper_HLT_HT500_v9;
int	triggerresultshelper_HLT_HT550_v1;
int	triggerresultshelper_HLT_HT550_v2;
int	triggerresultshelper_HLT_HT550_v3;
int	triggerresultshelper_HLT_HT550_v4;
int	triggerresultshelper_HLT_HT550_v5;
int	triggerresultshelper_HLT_HT550_v6;
int	triggerresultshelper_HLT_HT550_v7;
int	triggerresultshelper_HLT_HT550_v8;
int	triggerresultshelper_HLT_HT550_v9;
int	triggerresultshelper_HLT_HT650_v1;
int	triggerresultshelper_HLT_HT650_v2;
int	triggerresultshelper_HLT_HT650_v3;
int	triggerresultshelper_HLT_HT650_v4;
int	triggerresultshelper_HLT_HT650_v5;
int	triggerresultshelper_HLT_HT650_v6;
int	triggerresultshelper_HLT_HT650_v7;
int	triggerresultshelper_HLT_HT650_v8;
int	triggerresultshelper_HLT_HT650_v9;
int	triggerresultshelper_HLT_HT750_v1;
int	triggerresultshelper_HLT_HT750_v2;
int	triggerresultshelper_HLT_HT750_v3;
int	triggerresultshelper_HLT_HT750_v4;
int	triggerresultshelper_HLT_HT750_v5;
int	triggerresultshelper_HLT_HT750_v6;
int	triggerresultshelper_HLT_HT750_v7;
int	triggerresultshelper_HLT_HT750_v8;
int	triggerresultshelper_HLT_HT750_v9;
int	triggerresultshelper_HLT_PFHT650_v10;
int	triggerresultshelper_HLT_PFHT650_v11;
int	triggerresultshelper_HLT_PFHT650_v5;
int	triggerresultshelper_HLT_PFHT650_v6;
int	triggerresultshelper_HLT_PFHT650_v7;
int	triggerresultshelper_HLT_PFHT650_v8;
int	triggerresultshelper_HLT_PFHT650_v9;
int	triggerresultshelper_HLT_PFHT700_v10;
int	triggerresultshelper_HLT_PFHT700_v11;
int	triggerresultshelper_HLT_PFHT700_v5;
int	triggerresultshelper_HLT_PFHT700_v6;
int	triggerresultshelper_HLT_PFHT700_v7;
int	triggerresultshelper_HLT_PFHT700_v8;
int	triggerresultshelper_HLT_PFHT700_v9;
int	triggerresultshelper_HLT_PFHT750_v10;
int	triggerresultshelper_HLT_PFHT750_v11;
int	triggerresultshelper_HLT_PFHT750_v5;
int	triggerresultshelper_HLT_PFHT750_v6;
int	triggerresultshelper_HLT_PFHT750_v7;
int	triggerresultshelper_HLT_PFHT750_v8;
int	triggerresultshelper_HLT_PFHT750_v9;
int	triggerresultshelper_HLT_PFJet320_v10;
int	triggerresultshelper_HLT_PFJet320_v11;
int	triggerresultshelper_HLT_PFJet320_v3;
int	triggerresultshelper_HLT_PFJet320_v4;
int	triggerresultshelper_HLT_PFJet320_v5;
int	triggerresultshelper_HLT_PFJet320_v6;
int	triggerresultshelper_HLT_PFJet320_v7;
int	triggerresultshelper_HLT_PFJet320_v8;
int	triggerresultshelper_HLT_PFJet320_v9;
int	triggerresultshelper_HLT_PFJet400_v10;
int	triggerresultshelper_HLT_PFJet400_v11;
int	triggerresultshelper_HLT_PFJet400_v3;
int	triggerresultshelper_HLT_PFJet400_v4;
int	triggerresultshelper_HLT_PFJet400_v5;
int	triggerresultshelper_HLT_PFJet400_v6;
int	triggerresultshelper_HLT_PFJet400_v7;
int	triggerresultshelper_HLT_PFJet400_v8;
int	triggerresultshelper_HLT_PFJet400_v9;
int	triggerresultshelper_HLT_PFNoPUHT650_v1;
int	triggerresultshelper_HLT_PFNoPUHT650_v2;
int	triggerresultshelper_HLT_PFNoPUHT650_v3;
int	triggerresultshelper_HLT_PFNoPUHT650_v4;
int	triggerresultshelper_HLT_PFNoPUHT650_v5;
int	triggerresultshelper_HLT_PFNoPUHT650_v6;
int	triggerresultshelper_HLT_PFNoPUHT700_v1;
int	triggerresultshelper_HLT_PFNoPUHT700_v2;
int	triggerresultshelper_HLT_PFNoPUHT700_v3;
int	triggerresultshelper_HLT_PFNoPUHT700_v4;
int	triggerresultshelper_HLT_PFNoPUHT700_v5;
int	triggerresultshelper_HLT_PFNoPUHT700_v6;
int	triggerresultshelper_HLT_PFNoPUHT750_v1;
int	triggerresultshelper_HLT_PFNoPUHT750_v2;
int	triggerresultshelper_HLT_PFNoPUHT750_v3;
int	triggerresultshelper_HLT_PFNoPUHT750_v4;
int	triggerresultshelper_HLT_PFNoPUHT750_v5;
int	triggerresultshelper_HLT_PFNoPUHT750_v6;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v10;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v11;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v12;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v2;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v3;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v4;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v5;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v6;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v7;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v8;
int	triggerresultshelper_prescaleHLT_DiPFJetAve320_v9;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v10;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v11;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v12;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v2;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v3;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v4;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v5;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v6;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v7;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v8;
int	triggerresultshelper_prescaleHLT_DiPFJetAve400_v9;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v11;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v12;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v6;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v7;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v8;
int	triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v9;
int	triggerresultshelper_prescaleHLT_HT450_v1;
int	triggerresultshelper_prescaleHLT_HT450_v2;
int	triggerresultshelper_prescaleHLT_HT450_v3;
int	triggerresultshelper_prescaleHLT_HT450_v4;
int	triggerresultshelper_prescaleHLT_HT450_v5;
int	triggerresultshelper_prescaleHLT_HT450_v6;
int	triggerresultshelper_prescaleHLT_HT450_v7;
int	triggerresultshelper_prescaleHLT_HT450_v8;
int	triggerresultshelper_prescaleHLT_HT450_v9;
int	triggerresultshelper_prescaleHLT_HT500_v1;
int	triggerresultshelper_prescaleHLT_HT500_v2;
int	triggerresultshelper_prescaleHLT_HT500_v3;
int	triggerresultshelper_prescaleHLT_HT500_v4;
int	triggerresultshelper_prescaleHLT_HT500_v5;
int	triggerresultshelper_prescaleHLT_HT500_v6;
int	triggerresultshelper_prescaleHLT_HT500_v7;
int	triggerresultshelper_prescaleHLT_HT500_v8;
int	triggerresultshelper_prescaleHLT_HT500_v9;
int	triggerresultshelper_prescaleHLT_HT550_v1;
int	triggerresultshelper_prescaleHLT_HT550_v2;
int	triggerresultshelper_prescaleHLT_HT550_v3;
int	triggerresultshelper_prescaleHLT_HT550_v4;
int	triggerresultshelper_prescaleHLT_HT550_v5;
int	triggerresultshelper_prescaleHLT_HT550_v6;
int	triggerresultshelper_prescaleHLT_HT550_v7;
int	triggerresultshelper_prescaleHLT_HT550_v8;
int	triggerresultshelper_prescaleHLT_HT550_v9;
int	triggerresultshelper_prescaleHLT_HT650_v1;
int	triggerresultshelper_prescaleHLT_HT650_v2;
int	triggerresultshelper_prescaleHLT_HT650_v3;
int	triggerresultshelper_prescaleHLT_HT650_v4;
int	triggerresultshelper_prescaleHLT_HT650_v5;
int	triggerresultshelper_prescaleHLT_HT650_v6;
int	triggerresultshelper_prescaleHLT_HT650_v7;
int	triggerresultshelper_prescaleHLT_HT650_v8;
int	triggerresultshelper_prescaleHLT_HT650_v9;
int	triggerresultshelper_prescaleHLT_HT750_v1;
int	triggerresultshelper_prescaleHLT_HT750_v2;
int	triggerresultshelper_prescaleHLT_HT750_v3;
int	triggerresultshelper_prescaleHLT_HT750_v4;
int	triggerresultshelper_prescaleHLT_HT750_v5;
int	triggerresultshelper_prescaleHLT_HT750_v6;
int	triggerresultshelper_prescaleHLT_HT750_v7;
int	triggerresultshelper_prescaleHLT_HT750_v8;
int	triggerresultshelper_prescaleHLT_HT750_v9;
int	triggerresultshelper_prescaleHLT_PFHT650_v10;
int	triggerresultshelper_prescaleHLT_PFHT650_v11;
int	triggerresultshelper_prescaleHLT_PFHT650_v5;
int	triggerresultshelper_prescaleHLT_PFHT650_v6;
int	triggerresultshelper_prescaleHLT_PFHT650_v7;
int	triggerresultshelper_prescaleHLT_PFHT650_v8;
int	triggerresultshelper_prescaleHLT_PFHT650_v9;
int	triggerresultshelper_prescaleHLT_PFHT700_v10;
int	triggerresultshelper_prescaleHLT_PFHT700_v11;
int	triggerresultshelper_prescaleHLT_PFHT700_v5;
int	triggerresultshelper_prescaleHLT_PFHT700_v6;
int	triggerresultshelper_prescaleHLT_PFHT700_v7;
int	triggerresultshelper_prescaleHLT_PFHT700_v8;
int	triggerresultshelper_prescaleHLT_PFHT700_v9;
int	triggerresultshelper_prescaleHLT_PFHT750_v10;
int	triggerresultshelper_prescaleHLT_PFHT750_v11;
int	triggerresultshelper_prescaleHLT_PFHT750_v5;
int	triggerresultshelper_prescaleHLT_PFHT750_v6;
int	triggerresultshelper_prescaleHLT_PFHT750_v7;
int	triggerresultshelper_prescaleHLT_PFHT750_v8;
int	triggerresultshelper_prescaleHLT_PFHT750_v9;
int	triggerresultshelper_prescaleHLT_PFJet320_v10;
int	triggerresultshelper_prescaleHLT_PFJet320_v11;
int	triggerresultshelper_prescaleHLT_PFJet320_v3;
int	triggerresultshelper_prescaleHLT_PFJet320_v4;
int	triggerresultshelper_prescaleHLT_PFJet320_v5;
int	triggerresultshelper_prescaleHLT_PFJet320_v6;
int	triggerresultshelper_prescaleHLT_PFJet320_v7;
int	triggerresultshelper_prescaleHLT_PFJet320_v8;
int	triggerresultshelper_prescaleHLT_PFJet320_v9;
int	triggerresultshelper_prescaleHLT_PFJet400_v10;
int	triggerresultshelper_prescaleHLT_PFJet400_v11;
int	triggerresultshelper_prescaleHLT_PFJet400_v3;
int	triggerresultshelper_prescaleHLT_PFJet400_v4;
int	triggerresultshelper_prescaleHLT_PFJet400_v5;
int	triggerresultshelper_prescaleHLT_PFJet400_v6;
int	triggerresultshelper_prescaleHLT_PFJet400_v7;
int	triggerresultshelper_prescaleHLT_PFJet400_v8;
int	triggerresultshelper_prescaleHLT_PFJet400_v9;
int	triggerresultshelper_prescaleHLT_PFNoPUHT650_v1;
int	triggerresultshelper_prescaleHLT_PFNoPUHT650_v2;
int	triggerresultshelper_prescaleHLT_PFNoPUHT650_v3;
int	triggerresultshelper_prescaleHLT_PFNoPUHT650_v4;
int	triggerresultshelper_prescaleHLT_PFNoPUHT650_v5;
int	triggerresultshelper_prescaleHLT_PFNoPUHT650_v6;
int	triggerresultshelper_prescaleHLT_PFNoPUHT700_v1;
int	triggerresultshelper_prescaleHLT_PFNoPUHT700_v2;
int	triggerresultshelper_prescaleHLT_PFNoPUHT700_v3;
int	triggerresultshelper_prescaleHLT_PFNoPUHT700_v4;
int	triggerresultshelper_prescaleHLT_PFNoPUHT700_v5;
int	triggerresultshelper_prescaleHLT_PFNoPUHT700_v6;
int	triggerresultshelper_prescaleHLT_PFNoPUHT750_v1;
int	triggerresultshelper_prescaleHLT_PFNoPUHT750_v2;
int	triggerresultshelper_prescaleHLT_PFNoPUHT750_v3;
int	triggerresultshelper_prescaleHLT_PFNoPUHT750_v4;
int	triggerresultshelper_prescaleHLT_PFNoPUHT750_v5;
int	triggerresultshelper_prescaleHLT_PFNoPUHT750_v6;
std::vector<double>	vertex_chi2(2,0);
std::vector<int>	vertex_isFake(2,0);
std::vector<double>	vertex_ndof(2,0);
std::vector<double>	vertex_position_Rho(2,0);
std::vector<double>	vertex_x(2,0);
std::vector<double>	vertex_xError(2,0);
std::vector<double>	vertex_y(2,0);
std::vector<double>	vertex_yError(2,0);
std::vector<double>	vertex_z(2,0);
std::vector<double>	vertex_zError(2,0);

//-----------------------------------------------------------------------------
// --- Structs can be filled by calling fillObjects()
// --- after the call to stream.read(...)
//-----------------------------------------------------------------------------
struct calomet_s
{
  bool	selected;
  double	energy;
  double	pt;
  double	phi;
  double	sumEt;
  double	mEtSig;
  double	significance;
};
std::vector<calomet_s> calomet(200);

std::ostream& operator<<(std::ostream& os, const calomet_s& o)
{
  char r[1024];
  os << "calomet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct calomet1_s
{
  bool	selected;
  double	energy;
  double	pt;
  double	phi;
  double	sumEt;
  double	mEtSig;
  double	significance;
};
std::vector<calomet1_s> calomet1(200);

std::ostream& operator<<(std::ostream& os, const calomet1_s& o)
{
  char r[1024];
  os << "calomet1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgbasemet_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	sumEt;
};
std::vector<cmgbasemet_s> cmgbasemet(200);

std::ostream& operator<<(std::ostream& os, const cmgbasemet_s& o)
{
  char r[1024];
  os << "cmgbasemet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgbasemet1_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	sumEt;
};
std::vector<cmgbasemet1_s> cmgbasemet1(200);

std::ostream& operator<<(std::ostream& os, const cmgbasemet1_s& o)
{
  char r[1024];
  os << "cmgbasemet1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgbasemet2_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	sumEt;
};
std::vector<cmgbasemet2_s> cmgbasemet2(200);

std::ostream& operator<<(std::ostream& os, const cmgbasemet2_s& o)
{
  char r[1024];
  os << "cmgbasemet2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgelectron_s
{
  bool	selected;
  int	charge;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<cmgelectron_s> cmgelectron(200);

std::ostream& operator<<(std::ostream& os, const cmgelectron_s& o)
{
  char r[1024];
  os << "cmgelectron" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgelectron1_s
{
  bool	selected;
  int	charge;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<cmgelectron1_s> cmgelectron1(200);

std::ostream& operator<<(std::ostream& os, const cmgelectron1_s& o)
{
  char r[1024];
  os << "cmgelectron1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgmuon_s
{
  bool	selected;
  int	charge;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<cmgmuon_s> cmgmuon(200);

std::ostream& operator<<(std::ostream& os, const cmgmuon_s& o)
{
  char r[1024];
  os << "cmgmuon" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgmuon1_s
{
  bool	selected;
  int	charge;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<cmgmuon1_s> cmgmuon1(200);

std::ostream& operator<<(std::ostream& os, const cmgmuon1_s& o)
{
  char r[1024];
  os << "cmgmuon1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgpfjet_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  float	jetArea;
  double	mass;
  double	rapidity;
  int	nConstituents;
  int	partonFlavour;
  double	component_0_fraction;
  double	component_0_number;
  double	component_1_fraction;
  double	component_1_number;
  double	component_2_fraction;
  double	component_2_number;
  double	component_3_fraction;
  double	component_3_number;
  double	component_4_fraction;
  double	component_4_number;
  double	component_5_fraction;
  double	component_5_number;
  double	component_6_fraction;
  double	component_6_number;
  double	component_7_fraction;
  double	component_7_number;
  double	trackCountingHighEffBJetTag;
  double	trackCountingHighPurBJetTags;
  double	jetProbabilityBJetTags;
  double	jetBProbabilityBJetTags;
  double	combinedSecondaryVertexBJetTags;
  double	combinedSecondaryVertexMVABJetTags;
};
std::vector<cmgpfjet_s> cmgpfjet(200);

std::ostream& operator<<(std::ostream& os, const cmgpfjet_s& o)
{
  char r[1024];
  os << "cmgpfjet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "component_0_fraction", (double)o.component_0_fraction); os << r;
  sprintf(r, "  %-32s: %f\n", "component_0_number", (double)o.component_0_number); os << r;
  sprintf(r, "  %-32s: %f\n", "component_1_fraction", (double)o.component_1_fraction); os << r;
  sprintf(r, "  %-32s: %f\n", "component_1_number", (double)o.component_1_number); os << r;
  sprintf(r, "  %-32s: %f\n", "component_2_fraction", (double)o.component_2_fraction); os << r;
  sprintf(r, "  %-32s: %f\n", "component_2_number", (double)o.component_2_number); os << r;
  sprintf(r, "  %-32s: %f\n", "component_3_fraction", (double)o.component_3_fraction); os << r;
  sprintf(r, "  %-32s: %f\n", "component_3_number", (double)o.component_3_number); os << r;
  sprintf(r, "  %-32s: %f\n", "component_4_fraction", (double)o.component_4_fraction); os << r;
  sprintf(r, "  %-32s: %f\n", "component_4_number", (double)o.component_4_number); os << r;
  sprintf(r, "  %-32s: %f\n", "component_5_fraction", (double)o.component_5_fraction); os << r;
  sprintf(r, "  %-32s: %f\n", "component_5_number", (double)o.component_5_number); os << r;
  sprintf(r, "  %-32s: %f\n", "component_6_fraction", (double)o.component_6_fraction); os << r;
  sprintf(r, "  %-32s: %f\n", "component_6_number", (double)o.component_6_number); os << r;
  sprintf(r, "  %-32s: %f\n", "component_7_fraction", (double)o.component_7_fraction); os << r;
  sprintf(r, "  %-32s: %f\n", "component_7_number", (double)o.component_7_number); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTag", (double)o.trackCountingHighEffBJetTag); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgtau_s
{
  bool	selected;
  int	charge;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  float	byLooseCombinedIsolationDeltaBetaCorr;
};
std::vector<cmgtau_s> cmgtau(200);

std::ostream& operator<<(std::ostream& os, const cmgtau_s& o)
{
  char r[1024];
  os << "cmgtau" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "byLooseCombinedIsolationDeltaBetaCorr", (double)o.byLooseCombinedIsolationDeltaBetaCorr); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct cmgtau1_s
{
  bool	selected;
  int	charge;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  float	byTightCombinedIsolationDeltaBetaCorr;
};
std::vector<cmgtau1_s> cmgtau1(200);

std::ostream& operator<<(std::ostream& os, const cmgtau1_s& o)
{
  char r[1024];
  os << "cmgtau1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "byTightCombinedIsolationDeltaBetaCorr", (double)o.byTightCombinedIsolationDeltaBetaCorr); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct electronhelper_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  int	charge;
  float	eSuperClusterOverP;
  float	deltaEtaSuperClusterTrackAtVtx;
  float	deltaPhiSuperClusterTrackAtVtx;
  int	isEB;
  int	isEE;
  float	scSigmaIEtaIEta;
  float	hadronicOverEm;
  double	superCluster_energy;
  unsigned short	gsfTrack_trackerExpectedHitsInner_numberOfHits;
  float	simpleEleId80relIso;
  float	simpleEleId95relIso;
  float	trackIso;
  float	ecalIso;
  float	hcalIso;
  float	caloIso;
  double	dxywrtPV;
  double	dzwrtPV;
};
std::vector<electronhelper_s> electronhelper(200);

std::ostream& operator<<(std::ostream& os, const electronhelper_s& o)
{
  char r[1024];
  os << "electronhelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "eSuperClusterOverP", (double)o.eSuperClusterOverP); os << r;
  sprintf(r, "  %-32s: %f\n", "deltaEtaSuperClusterTrackAtVtx", (double)o.deltaEtaSuperClusterTrackAtVtx); os << r;
  sprintf(r, "  %-32s: %f\n", "deltaPhiSuperClusterTrackAtVtx", (double)o.deltaPhiSuperClusterTrackAtVtx); os << r;
  sprintf(r, "  %-32s: %f\n", "isEB", (double)o.isEB); os << r;
  sprintf(r, "  %-32s: %f\n", "isEE", (double)o.isEE); os << r;
  sprintf(r, "  %-32s: %f\n", "scSigmaIEtaIEta", (double)o.scSigmaIEtaIEta); os << r;
  sprintf(r, "  %-32s: %f\n", "hadronicOverEm", (double)o.hadronicOverEm); os << r;
  sprintf(r, "  %-32s: %f\n", "superCluster_energy", (double)o.superCluster_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "gsfTrack_trackerExpectedHitsInner_numberOfHits", (double)o.gsfTrack_trackerExpectedHitsInner_numberOfHits); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleEleId80relIso", (double)o.simpleEleId80relIso); os << r;
  sprintf(r, "  %-32s: %f\n", "simpleEleId95relIso", (double)o.simpleEleId95relIso); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "caloIso", (double)o.caloIso); os << r;
  sprintf(r, "  %-32s: %f\n", "dxywrtPV", (double)o.dxywrtPV); os << r;
  sprintf(r, "  %-32s: %f\n", "dzwrtPV", (double)o.dzwrtPV); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct genjet_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  double	mass;
  double	rapidity;
  int	nConstituents;
};
std::vector<genjet_s> genjet(200);

std::ostream& operator<<(std::ostream& os, const genjet_s& o)
{
  char r[1024];
  os << "genjet" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct genparticlehelper_s
{
  bool	selected;
  int	firstMother;
  int	lastMother;
  int	firstDaughter;
  int	lastDaughter;
  int	charge;
  int	pdgId;
  int	status;
  double	pt;
  double	eta;
  double	phi;
  double	mass;
};
std::vector<genparticlehelper_s> genparticlehelper(500);

std::ostream& operator<<(std::ostream& os, const genparticlehelper_s& o)
{
  char r[1024];
  os << "genparticlehelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "firstMother", (double)o.firstMother); os << r;
  sprintf(r, "  %-32s: %f\n", "lastMother", (double)o.lastMother); os << r;
  sprintf(r, "  %-32s: %f\n", "firstDaughter", (double)o.firstDaughter); os << r;
  sprintf(r, "  %-32s: %f\n", "lastDaughter", (double)o.lastDaughter); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "pdgId", (double)o.pdgId); os << r;
  sprintf(r, "  %-32s: %f\n", "status", (double)o.status); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper_s
{
  bool	selected;
  double	energy;
  double	uncor_energy;
  double	et;
  double	uncor_et;
  double	pt;
  double	uncor_pt;
  double	phi;
  double	eta;
  double	rapidity;
  double	mass;
  float	jetArea;
  float	jetCharge03;
  float	jetCharge05;
  float	jetCharge10;
  float	chargedHadronEnergyFraction;
  float	neutralHadronEnergyFraction;
  float	chargedEmEnergyFraction;
  float	neutralEmEnergyFraction;
  float	photonEnergyFraction;
  float	muonEnergyFraction;
  int	chargedMultiplicity;
  int	nConstituents;
  int	partonFlavour;
  float	trackCountingHighEffBJetTags;
  float	trackCountingHighPurBJetTags;
  float	jetProbabilityBJetTags;
  float	jetBProbabilityBJetTags;
  float	combinedSecondaryVertexBJetTags;
  float	combinedSecondaryVertexMVABJetTags;
};
std::vector<jethelper_s> jethelper(200);

std::ostream& operator<<(std::ostream& os, const jethelper_s& o)
{
  char r[1024];
  os << "jethelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_et", (double)o.uncor_et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge03", (double)o.jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge05", (double)o.jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge10", (double)o.jetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "photonEnergyFraction", (double)o.photonEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "muonEnergyFraction", (double)o.muonEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper1_s
{
  bool	selected;
  double	energy;
  double	uncor_energy;
  double	et;
  double	uncor_et;
  double	pt;
  double	uncor_pt;
  double	phi;
  double	eta;
  double	rapidity;
  double	mass;
  float	jetArea;
  float	jetCharge03;
  float	jetCharge05;
  float	jetCharge10;
  float	chargedHadronEnergyFraction;
  float	neutralHadronEnergyFraction;
  float	chargedEmEnergyFraction;
  float	neutralEmEnergyFraction;
  float	photonEnergyFraction;
  float	muonEnergyFraction;
  int	chargedMultiplicity;
  int	nConstituents;
  int	partonFlavour;
  float	trackCountingHighEffBJetTags;
  float	trackCountingHighPurBJetTags;
  float	jetProbabilityBJetTags;
  float	jetBProbabilityBJetTags;
  float	combinedSecondaryVertexBJetTags;
  float	combinedSecondaryVertexMVABJetTags;
};
std::vector<jethelper1_s> jethelper1(200);

std::ostream& operator<<(std::ostream& os, const jethelper1_s& o)
{
  char r[1024];
  os << "jethelper1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_et", (double)o.uncor_et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge03", (double)o.jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge05", (double)o.jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge10", (double)o.jetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "photonEnergyFraction", (double)o.photonEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "muonEnergyFraction", (double)o.muonEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper2_s
{
  bool	selected;
  double	energy;
  double	uncor_energy;
  double	et;
  double	uncor_et;
  double	pt;
  double	uncor_pt;
  double	phi;
  double	eta;
  double	rapidity;
  double	mass;
  float	jetArea;
  float	jetCharge03;
  float	jetCharge05;
  float	jetCharge10;
  int	nConstituents;
  int	partonFlavour;
  float	trackCountingHighEffBJetTags;
  float	trackCountingHighPurBJetTags;
  float	jetProbabilityBJetTags;
  float	jetBProbabilityBJetTags;
  float	combinedSecondaryVertexBJetTags;
  float	combinedSecondaryVertexMVABJetTags;
  size_t	numberOfDaughters;
  double	daughter_0_energy;
  double	daughter_0_pt;
  double	daughter_0_eta;
  double	daughter_0_rapidity;
  double	daughter_0_phi;
  double	daughter_0_mass;
  double	daughter_1_energy;
  double	daughter_1_pt;
  double	daughter_1_eta;
  double	daughter_1_rapidity;
  double	daughter_1_phi;
  double	daughter_1_mass;
  float	daughter_0_jetCharge03;
  float	daughter_0_jetCharge05;
  float	daughter_0_jetCharge10;
  float	daughter_1_jetCharge03;
  float	daughter_1_jetCharge05;
  float	daughter_1_jetCharge10;
};
std::vector<jethelper2_s> jethelper2(200);

std::ostream& operator<<(std::ostream& os, const jethelper2_s& o)
{
  char r[1024];
  os << "jethelper2" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_et", (double)o.uncor_et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge03", (double)o.jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge05", (double)o.jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge10", (double)o.jetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfDaughters", (double)o.numberOfDaughters); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_energy", (double)o.daughter_0_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_pt", (double)o.daughter_0_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_eta", (double)o.daughter_0_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_rapidity", (double)o.daughter_0_rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_phi", (double)o.daughter_0_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_mass", (double)o.daughter_0_mass); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_energy", (double)o.daughter_1_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_pt", (double)o.daughter_1_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_eta", (double)o.daughter_1_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_rapidity", (double)o.daughter_1_rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_phi", (double)o.daughter_1_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_mass", (double)o.daughter_1_mass); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_jetCharge03", (double)o.daughter_0_jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_jetCharge05", (double)o.daughter_0_jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_jetCharge10", (double)o.daughter_0_jetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_jetCharge03", (double)o.daughter_1_jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_jetCharge05", (double)o.daughter_1_jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_jetCharge10", (double)o.daughter_1_jetCharge10); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper3_s
{
  bool	selected;
  double	energy;
  double	uncor_energy;
  double	et;
  double	uncor_et;
  double	pt;
  double	uncor_pt;
  double	phi;
  double	eta;
  double	rapidity;
  double	mass;
  float	jetArea;
  float	jetCharge03;
  float	jetCharge05;
  float	jetCharge10;
  float	chargedHadronEnergyFraction;
  float	neutralHadronEnergyFraction;
  float	chargedEmEnergyFraction;
  float	neutralEmEnergyFraction;
  float	photonEnergyFraction;
  float	muonEnergyFraction;
  int	chargedMultiplicity;
  int	nConstituents;
  int	partonFlavour;
  float	trackCountingHighEffBJetTags;
  float	trackCountingHighPurBJetTags;
  float	jetProbabilityBJetTags;
  float	jetBProbabilityBJetTags;
  float	combinedSecondaryVertexBJetTags;
  float	combinedSecondaryVertexMVABJetTags;
  float	qjetsvolatility;
  float	tau1;
  float	tau2;
  float	tau3;
  float	C2beta17;
};
std::vector<jethelper3_s> jethelper3(200);

std::ostream& operator<<(std::ostream& os, const jethelper3_s& o)
{
  char r[1024];
  os << "jethelper3" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_et", (double)o.uncor_et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge03", (double)o.jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge05", (double)o.jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge10", (double)o.jetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "photonEnergyFraction", (double)o.photonEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "muonEnergyFraction", (double)o.muonEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "qjetsvolatility", (double)o.qjetsvolatility); os << r;
  sprintf(r, "  %-32s: %f\n", "tau1", (double)o.tau1); os << r;
  sprintf(r, "  %-32s: %f\n", "tau2", (double)o.tau2); os << r;
  sprintf(r, "  %-32s: %f\n", "tau3", (double)o.tau3); os << r;
  sprintf(r, "  %-32s: %f\n", "C2beta17", (double)o.C2beta17); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper4_s
{
  bool	selected;
  double	energy;
  double	uncor_energy;
  double	et;
  double	uncor_et;
  double	pt;
  double	uncor_pt;
  double	phi;
  double	eta;
  double	rapidity;
  double	mass;
  float	jetArea;
  float	jetCharge03;
  float	jetCharge05;
  float	jetCharge10;
  int	nConstituents;
  int	partonFlavour;
  float	trackCountingHighEffBJetTags;
  float	trackCountingHighPurBJetTags;
  float	jetProbabilityBJetTags;
  float	jetBProbabilityBJetTags;
  float	combinedSecondaryVertexBJetTags;
  float	combinedSecondaryVertexMVABJetTags;
  size_t	numberOfDaughters;
  double	daughter_0_energy;
  double	daughter_0_pt;
  double	daughter_0_eta;
  double	daughter_0_rapidity;
  double	daughter_0_phi;
  double	daughter_0_mass;
  double	daughter_1_energy;
  double	daughter_1_pt;
  double	daughter_1_eta;
  double	daughter_1_rapidity;
  double	daughter_1_phi;
  double	daughter_1_mass;
  float	tau1;
  float	tau2;
  float	tau3;
  float	daughter_0_jetCharge03;
  float	daughter_0_jetCharge05;
  float	daughter_0_jetCharge10;
  float	daughter_1_jetCharge03;
  float	daughter_1_jetCharge05;
  float	daughter_1_jetCharge10;
};
std::vector<jethelper4_s> jethelper4(200);

std::ostream& operator<<(std::ostream& os, const jethelper4_s& o)
{
  char r[1024];
  os << "jethelper4" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_et", (double)o.uncor_et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge03", (double)o.jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge05", (double)o.jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge10", (double)o.jetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfDaughters", (double)o.numberOfDaughters); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_energy", (double)o.daughter_0_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_pt", (double)o.daughter_0_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_eta", (double)o.daughter_0_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_rapidity", (double)o.daughter_0_rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_phi", (double)o.daughter_0_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_mass", (double)o.daughter_0_mass); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_energy", (double)o.daughter_1_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_pt", (double)o.daughter_1_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_eta", (double)o.daughter_1_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_rapidity", (double)o.daughter_1_rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_phi", (double)o.daughter_1_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_mass", (double)o.daughter_1_mass); os << r;
  sprintf(r, "  %-32s: %f\n", "tau1", (double)o.tau1); os << r;
  sprintf(r, "  %-32s: %f\n", "tau2", (double)o.tau2); os << r;
  sprintf(r, "  %-32s: %f\n", "tau3", (double)o.tau3); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_jetCharge03", (double)o.daughter_0_jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_jetCharge05", (double)o.daughter_0_jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_jetCharge10", (double)o.daughter_0_jetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_jetCharge03", (double)o.daughter_1_jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_jetCharge05", (double)o.daughter_1_jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_jetCharge10", (double)o.daughter_1_jetCharge10); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper5_s
{
  bool	selected;
  double	energy;
  double	uncor_energy;
  double	et;
  double	uncor_et;
  double	pt;
  double	uncor_pt;
  double	phi;
  double	eta;
  double	rapidity;
  double	mass;
  float	jetArea;
  float	jetCharge03;
  float	jetCharge05;
  float	jetCharge10;
  float	chargedHadronEnergyFraction;
  float	neutralHadronEnergyFraction;
  float	chargedEmEnergyFraction;
  float	neutralEmEnergyFraction;
  float	photonEnergyFraction;
  float	muonEnergyFraction;
  int	chargedMultiplicity;
  int	nConstituents;
  int	partonFlavour;
  float	trackCountingHighEffBJetTags;
  float	trackCountingHighPurBJetTags;
  float	jetProbabilityBJetTags;
  float	jetBProbabilityBJetTags;
  float	combinedSecondaryVertexBJetTags;
  float	combinedSecondaryVertexMVABJetTags;
  float	qjetsvolatility;
  float	tau1;
  float	tau2;
  float	tau3;
  float	C2beta17;
  int	getNcharged01;
  int	getNneutral01;
  int	getChargedPt0;
  int	getChargedPt1;
  int	getChargedPt2;
  int	getChargedPt3;
  int	getPt0;
  int	getPt1;
  int	getPt2;
  int	getPt3;
};
std::vector<jethelper5_s> jethelper5(200);

std::ostream& operator<<(std::ostream& os, const jethelper5_s& o)
{
  char r[1024];
  os << "jethelper5" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_energy", (double)o.uncor_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_et", (double)o.uncor_et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "uncor_pt", (double)o.uncor_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge03", (double)o.jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge05", (double)o.jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge10", (double)o.jetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedHadronEnergyFraction", (double)o.chargedHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralHadronEnergyFraction", (double)o.neutralHadronEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedEmEnergyFraction", (double)o.chargedEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "neutralEmEnergyFraction", (double)o.neutralEmEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "photonEnergyFraction", (double)o.photonEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "muonEnergyFraction", (double)o.muonEnergyFraction); os << r;
  sprintf(r, "  %-32s: %f\n", "chargedMultiplicity", (double)o.chargedMultiplicity); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighEffBJetTags", (double)o.trackCountingHighEffBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "trackCountingHighPurBJetTags", (double)o.trackCountingHighPurBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetProbabilityBJetTags", (double)o.jetProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "jetBProbabilityBJetTags", (double)o.jetBProbabilityBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexBJetTags", (double)o.combinedSecondaryVertexBJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "combinedSecondaryVertexMVABJetTags", (double)o.combinedSecondaryVertexMVABJetTags); os << r;
  sprintf(r, "  %-32s: %f\n", "qjetsvolatility", (double)o.qjetsvolatility); os << r;
  sprintf(r, "  %-32s: %f\n", "tau1", (double)o.tau1); os << r;
  sprintf(r, "  %-32s: %f\n", "tau2", (double)o.tau2); os << r;
  sprintf(r, "  %-32s: %f\n", "tau3", (double)o.tau3); os << r;
  sprintf(r, "  %-32s: %f\n", "C2beta17", (double)o.C2beta17); os << r;
  sprintf(r, "  %-32s: %f\n", "getNcharged01", (double)o.getNcharged01); os << r;
  sprintf(r, "  %-32s: %f\n", "getNneutral01", (double)o.getNneutral01); os << r;
  sprintf(r, "  %-32s: %f\n", "getChargedPt0", (double)o.getChargedPt0); os << r;
  sprintf(r, "  %-32s: %f\n", "getChargedPt1", (double)o.getChargedPt1); os << r;
  sprintf(r, "  %-32s: %f\n", "getChargedPt2", (double)o.getChargedPt2); os << r;
  sprintf(r, "  %-32s: %f\n", "getChargedPt3", (double)o.getChargedPt3); os << r;
  sprintf(r, "  %-32s: %f\n", "getPt0", (double)o.getPt0); os << r;
  sprintf(r, "  %-32s: %f\n", "getPt1", (double)o.getPt1); os << r;
  sprintf(r, "  %-32s: %f\n", "getPt2", (double)o.getPt2); os << r;
  sprintf(r, "  %-32s: %f\n", "getPt3", (double)o.getPt3); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper6_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  double	rapidity;
  double	mass;
  float	jetArea;
  float	genJetCharge03;
  float	genJetCharge05;
  float	genJetCharge10;
  int	nConstituents;
  int	partonFlavour;
  float	genC2beta17;
  float	genC2beta17CHS;
  float	genTau1;
  float	genTau2;
  float	genTau3;
  float	genTau1Pt2;
  float	genTau2Pt2;
  float	genTau1Pt5;
  float	genTau2Pt5;
  float	genTau1CHS;
  float	genTau2CHS;
  float	genNCHS;
};
std::vector<jethelper6_s> jethelper6(200);

std::ostream& operator<<(std::ostream& os, const jethelper6_s& o)
{
  char r[1024];
  os << "jethelper6" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "genJetCharge03", (double)o.genJetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "genJetCharge05", (double)o.genJetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "genJetCharge10", (double)o.genJetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "genC2beta17", (double)o.genC2beta17); os << r;
  sprintf(r, "  %-32s: %f\n", "genC2beta17CHS", (double)o.genC2beta17CHS); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau1", (double)o.genTau1); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau2", (double)o.genTau2); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau3", (double)o.genTau3); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau1Pt2", (double)o.genTau1Pt2); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau2Pt2", (double)o.genTau2Pt2); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau1Pt5", (double)o.genTau1Pt5); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau2Pt5", (double)o.genTau2Pt5); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau1CHS", (double)o.genTau1CHS); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau2CHS", (double)o.genTau2CHS); os << r;
  sprintf(r, "  %-32s: %f\n", "genNCHS", (double)o.genNCHS); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct jethelper7_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  double	rapidity;
  double	mass;
  float	jetArea;
  float	jetCharge03;
  float	jetCharge05;
  float	jetCharge10;
  int	nConstituents;
  int	partonFlavour;
  size_t	numberOfDaughters;
  double	daughter_0_energy;
  double	daughter_0_pt;
  double	daughter_0_eta;
  double	daughter_0_rapidity;
  double	daughter_0_phi;
  double	daughter_0_mass;
  double	daughter_1_energy;
  double	daughter_1_pt;
  double	daughter_1_eta;
  double	daughter_1_rapidity;
  double	daughter_1_phi;
  double	daughter_1_mass;
  float	genTau1;
  float	genTau2;
  float	genTau3;
};
std::vector<jethelper7_s> jethelper7(200);

std::ostream& operator<<(std::ostream& os, const jethelper7_s& o)
{
  char r[1024];
  os << "jethelper7" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "rapidity", (double)o.rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "mass", (double)o.mass); os << r;
  sprintf(r, "  %-32s: %f\n", "jetArea", (double)o.jetArea); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge03", (double)o.jetCharge03); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge05", (double)o.jetCharge05); os << r;
  sprintf(r, "  %-32s: %f\n", "jetCharge10", (double)o.jetCharge10); os << r;
  sprintf(r, "  %-32s: %f\n", "nConstituents", (double)o.nConstituents); os << r;
  sprintf(r, "  %-32s: %f\n", "partonFlavour", (double)o.partonFlavour); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfDaughters", (double)o.numberOfDaughters); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_energy", (double)o.daughter_0_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_pt", (double)o.daughter_0_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_eta", (double)o.daughter_0_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_rapidity", (double)o.daughter_0_rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_phi", (double)o.daughter_0_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_0_mass", (double)o.daughter_0_mass); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_energy", (double)o.daughter_1_energy); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_pt", (double)o.daughter_1_pt); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_eta", (double)o.daughter_1_eta); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_rapidity", (double)o.daughter_1_rapidity); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_phi", (double)o.daughter_1_phi); os << r;
  sprintf(r, "  %-32s: %f\n", "daughter_1_mass", (double)o.daughter_1_mass); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau1", (double)o.genTau1); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau2", (double)o.genTau2); os << r;
  sprintf(r, "  %-32s: %f\n", "genTau3", (double)o.genTau3); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct leafcandidate_s
{
  bool	selected;
  int	charge;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
};
std::vector<leafcandidate_s> leafcandidate(200);

std::ostream& operator<<(std::ostream& os, const leafcandidate_s& o)
{
  char r[1024];
  os << "leafcandidate" << std::endl;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	sumEt;
  double	mEtSig;
  double	significance;
};
std::vector<met_s> met(200);

std::ostream& operator<<(std::ostream& os, const met_s& o)
{
  char r[1024];
  os << "met" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct met1_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	sumEt;
  double	mEtSig;
  double	significance;
};
std::vector<met1_s> met1(200);

std::ostream& operator<<(std::ostream& os, const met1_s& o)
{
  char r[1024];
  os << "met1" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "sumEt", (double)o.sumEt); os << r;
  sprintf(r, "  %-32s: %f\n", "mEtSig", (double)o.mEtSig); os << r;
  sprintf(r, "  %-32s: %f\n", "significance", (double)o.significance); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct muonhelper_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  int	charge;
  float	TMOneStationTight;
  int	isGlobalMuon;
  int	isTrackerMuon;
  int	isPFMuon;
  unsigned short	track_hitPattern_trackerLayersWithMeasurement;
  unsigned short	innerTrack_hitPattern_pixelLayersWithMeasurement;
  double	innerTrack_normalizedChi2;
  double	globalTrack_normalizedChi2;
  unsigned short	globalTrack_hitPattern_numberOfValidMuonHits;
  int	numberOfMatchedStations;
  double	dB;
  unsigned short	innerTrack_hitPattern_numberOfValidPixelHits;
  float	pfIsolationR04_sumChargedHadronPt;
  float	pfIsolationR04_sumChargedParticlePt;
  float	pfIsolationR04_sumNeutralHadronEt;
  float	pfIsolationR04_sumPhotonEt;
  float	pfIsolationR04_sumNeutralHadronEtHighThreshold;
  float	pfIsolationR04_sumPhotonEtHighThreshold;
  float	pfIsolationR04_sumPUPt;
  double	dxywrtPV;
  double	dzwrtPV;
};
std::vector<muonhelper_s> muonhelper(200);

std::ostream& operator<<(std::ostream& os, const muonhelper_s& o)
{
  char r[1024];
  os << "muonhelper" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "charge", (double)o.charge); os << r;
  sprintf(r, "  %-32s: %f\n", "TMOneStationTight", (double)o.TMOneStationTight); os << r;
  sprintf(r, "  %-32s: %f\n", "isGlobalMuon", (double)o.isGlobalMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isTrackerMuon", (double)o.isTrackerMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "isPFMuon", (double)o.isPFMuon); os << r;
  sprintf(r, "  %-32s: %f\n", "track_hitPattern_trackerLayersWithMeasurement", (double)o.track_hitPattern_trackerLayersWithMeasurement); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_hitPattern_pixelLayersWithMeasurement", (double)o.innerTrack_hitPattern_pixelLayersWithMeasurement); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_normalizedChi2", (double)o.innerTrack_normalizedChi2); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_normalizedChi2", (double)o.globalTrack_normalizedChi2); os << r;
  sprintf(r, "  %-32s: %f\n", "globalTrack_hitPattern_numberOfValidMuonHits", (double)o.globalTrack_hitPattern_numberOfValidMuonHits); os << r;
  sprintf(r, "  %-32s: %f\n", "numberOfMatchedStations", (double)o.numberOfMatchedStations); os << r;
  sprintf(r, "  %-32s: %f\n", "dB", (double)o.dB); os << r;
  sprintf(r, "  %-32s: %f\n", "innerTrack_hitPattern_numberOfValidPixelHits", (double)o.innerTrack_hitPattern_numberOfValidPixelHits); os << r;
  sprintf(r, "  %-32s: %f\n", "pfIsolationR04_sumChargedHadronPt", (double)o.pfIsolationR04_sumChargedHadronPt); os << r;
  sprintf(r, "  %-32s: %f\n", "pfIsolationR04_sumChargedParticlePt", (double)o.pfIsolationR04_sumChargedParticlePt); os << r;
  sprintf(r, "  %-32s: %f\n", "pfIsolationR04_sumNeutralHadronEt", (double)o.pfIsolationR04_sumNeutralHadronEt); os << r;
  sprintf(r, "  %-32s: %f\n", "pfIsolationR04_sumPhotonEt", (double)o.pfIsolationR04_sumPhotonEt); os << r;
  sprintf(r, "  %-32s: %f\n", "pfIsolationR04_sumNeutralHadronEtHighThreshold", (double)o.pfIsolationR04_sumNeutralHadronEtHighThreshold); os << r;
  sprintf(r, "  %-32s: %f\n", "pfIsolationR04_sumPhotonEtHighThreshold", (double)o.pfIsolationR04_sumPhotonEtHighThreshold); os << r;
  sprintf(r, "  %-32s: %f\n", "pfIsolationR04_sumPUPt", (double)o.pfIsolationR04_sumPUPt); os << r;
  sprintf(r, "  %-32s: %f\n", "dxywrtPV", (double)o.dxywrtPV); os << r;
  sprintf(r, "  %-32s: %f\n", "dzwrtPV", (double)o.dzwrtPV); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct pileupsummaryinfo_s
{
  bool	selected;
  int	getPU_NumInteractions;
  int	getBunchCrossing;
  float	getTrueNumInteractions;
};
std::vector<pileupsummaryinfo_s> pileupsummaryinfo(200);

std::ostream& operator<<(std::ostream& os, const pileupsummaryinfo_s& o)
{
  char r[1024];
  os << "pileupsummaryinfo" << std::endl;
  sprintf(r, "  %-32s: %f\n", "getPU_NumInteractions", (double)o.getPU_NumInteractions); os << r;
  sprintf(r, "  %-32s: %f\n", "getBunchCrossing", (double)o.getBunchCrossing); os << r;
  sprintf(r, "  %-32s: %f\n", "getTrueNumInteractions", (double)o.getTrueNumInteractions); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct tau_s
{
  bool	selected;
  double	energy;
  double	et;
  double	pt;
  double	phi;
  double	eta;
  float	byLooseCombinedIsolationDeltaBetaCorr;
  float	byMediumCombinedIsolationDeltaBetaCorr;
  float	trackIso;
  float	ecalIso;
  float	hcalIso;
  float	caloIso;
};
std::vector<tau_s> tau(200);

std::ostream& operator<<(std::ostream& os, const tau_s& o)
{
  char r[1024];
  os << "tau" << std::endl;
  sprintf(r, "  %-32s: %f\n", "energy", (double)o.energy); os << r;
  sprintf(r, "  %-32s: %f\n", "et", (double)o.et); os << r;
  sprintf(r, "  %-32s: %f\n", "pt", (double)o.pt); os << r;
  sprintf(r, "  %-32s: %f\n", "phi", (double)o.phi); os << r;
  sprintf(r, "  %-32s: %f\n", "eta", (double)o.eta); os << r;
  sprintf(r, "  %-32s: %f\n", "byLooseCombinedIsolationDeltaBetaCorr", (double)o.byLooseCombinedIsolationDeltaBetaCorr); os << r;
  sprintf(r, "  %-32s: %f\n", "byMediumCombinedIsolationDeltaBetaCorr", (double)o.byMediumCombinedIsolationDeltaBetaCorr); os << r;
  sprintf(r, "  %-32s: %f\n", "trackIso", (double)o.trackIso); os << r;
  sprintf(r, "  %-32s: %f\n", "ecalIso", (double)o.ecalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "hcalIso", (double)o.hcalIso); os << r;
  sprintf(r, "  %-32s: %f\n", "caloIso", (double)o.caloIso); os << r;
  return os;
}
//-----------------------------------------------------------------------------
struct vertex_s
{
  bool	selected;
  int	isFake;
  double	chi2;
  double	ndof;
  double	x;
  double	y;
  double	z;
  double	xError;
  double	yError;
  double	zError;
  double	position_Rho;
};
std::vector<vertex_s> vertex(2);

std::ostream& operator<<(std::ostream& os, const vertex_s& o)
{
  char r[1024];
  os << "vertex" << std::endl;
  sprintf(r, "  %-32s: %f\n", "isFake", (double)o.isFake); os << r;
  sprintf(r, "  %-32s: %f\n", "chi2", (double)o.chi2); os << r;
  sprintf(r, "  %-32s: %f\n", "ndof", (double)o.ndof); os << r;
  sprintf(r, "  %-32s: %f\n", "x", (double)o.x); os << r;
  sprintf(r, "  %-32s: %f\n", "y", (double)o.y); os << r;
  sprintf(r, "  %-32s: %f\n", "z", (double)o.z); os << r;
  sprintf(r, "  %-32s: %f\n", "xError", (double)o.xError); os << r;
  sprintf(r, "  %-32s: %f\n", "yError", (double)o.yError); os << r;
  sprintf(r, "  %-32s: %f\n", "zError", (double)o.zError); os << r;
  sprintf(r, "  %-32s: %f\n", "position_Rho", (double)o.position_Rho); os << r;
  return os;
}
//-----------------------------------------------------------------------------

static bool fillObjectsCalled = false;
void fillObjects()
{
  fillObjectsCalled = true;

  calomet.resize(calomet_energy.size());
  for(unsigned int i=0; i < calomet.size(); ++i)
    {
      calomet[i].selected	= false;
      calomet[i].energy	= calomet_energy[i];
      calomet[i].pt	= calomet_pt[i];
      calomet[i].phi	= calomet_phi[i];
      calomet[i].sumEt	= calomet_sumEt[i];
      calomet[i].mEtSig	= calomet_mEtSig[i];
      calomet[i].significance	= calomet_significance[i];
    }

  calomet1.resize(calomet1_energy.size());
  for(unsigned int i=0; i < calomet1.size(); ++i)
    {
      calomet1[i].selected	= false;
      calomet1[i].energy	= calomet1_energy[i];
      calomet1[i].pt	= calomet1_pt[i];
      calomet1[i].phi	= calomet1_phi[i];
      calomet1[i].sumEt	= calomet1_sumEt[i];
      calomet1[i].mEtSig	= calomet1_mEtSig[i];
      calomet1[i].significance	= calomet1_significance[i];
    }

  cmgbasemet.resize(cmgbasemet_energy.size());
  for(unsigned int i=0; i < cmgbasemet.size(); ++i)
    {
      cmgbasemet[i].selected	= false;
      cmgbasemet[i].energy	= cmgbasemet_energy[i];
      cmgbasemet[i].et	= cmgbasemet_et[i];
      cmgbasemet[i].pt	= cmgbasemet_pt[i];
      cmgbasemet[i].phi	= cmgbasemet_phi[i];
      cmgbasemet[i].sumEt	= cmgbasemet_sumEt[i];
    }

  cmgbasemet1.resize(cmgbasemet1_energy.size());
  for(unsigned int i=0; i < cmgbasemet1.size(); ++i)
    {
      cmgbasemet1[i].selected	= false;
      cmgbasemet1[i].energy	= cmgbasemet1_energy[i];
      cmgbasemet1[i].et	= cmgbasemet1_et[i];
      cmgbasemet1[i].pt	= cmgbasemet1_pt[i];
      cmgbasemet1[i].phi	= cmgbasemet1_phi[i];
      cmgbasemet1[i].sumEt	= cmgbasemet1_sumEt[i];
    }

  cmgbasemet2.resize(cmgbasemet2_energy.size());
  for(unsigned int i=0; i < cmgbasemet2.size(); ++i)
    {
      cmgbasemet2[i].selected	= false;
      cmgbasemet2[i].energy	= cmgbasemet2_energy[i];
      cmgbasemet2[i].et	= cmgbasemet2_et[i];
      cmgbasemet2[i].pt	= cmgbasemet2_pt[i];
      cmgbasemet2[i].phi	= cmgbasemet2_phi[i];
      cmgbasemet2[i].sumEt	= cmgbasemet2_sumEt[i];
    }

  cmgelectron.resize(cmgelectron_charge.size());
  for(unsigned int i=0; i < cmgelectron.size(); ++i)
    {
      cmgelectron[i].selected	= false;
      cmgelectron[i].charge	= cmgelectron_charge[i];
      cmgelectron[i].energy	= cmgelectron_energy[i];
      cmgelectron[i].et	= cmgelectron_et[i];
      cmgelectron[i].pt	= cmgelectron_pt[i];
      cmgelectron[i].phi	= cmgelectron_phi[i];
      cmgelectron[i].eta	= cmgelectron_eta[i];
    }

  cmgelectron1.resize(cmgelectron1_charge.size());
  for(unsigned int i=0; i < cmgelectron1.size(); ++i)
    {
      cmgelectron1[i].selected	= false;
      cmgelectron1[i].charge	= cmgelectron1_charge[i];
      cmgelectron1[i].energy	= cmgelectron1_energy[i];
      cmgelectron1[i].et	= cmgelectron1_et[i];
      cmgelectron1[i].pt	= cmgelectron1_pt[i];
      cmgelectron1[i].phi	= cmgelectron1_phi[i];
      cmgelectron1[i].eta	= cmgelectron1_eta[i];
    }

  cmgmuon.resize(cmgmuon_charge.size());
  for(unsigned int i=0; i < cmgmuon.size(); ++i)
    {
      cmgmuon[i].selected	= false;
      cmgmuon[i].charge	= cmgmuon_charge[i];
      cmgmuon[i].energy	= cmgmuon_energy[i];
      cmgmuon[i].et	= cmgmuon_et[i];
      cmgmuon[i].pt	= cmgmuon_pt[i];
      cmgmuon[i].phi	= cmgmuon_phi[i];
      cmgmuon[i].eta	= cmgmuon_eta[i];
    }

  cmgmuon1.resize(cmgmuon1_charge.size());
  for(unsigned int i=0; i < cmgmuon1.size(); ++i)
    {
      cmgmuon1[i].selected	= false;
      cmgmuon1[i].charge	= cmgmuon1_charge[i];
      cmgmuon1[i].energy	= cmgmuon1_energy[i];
      cmgmuon1[i].et	= cmgmuon1_et[i];
      cmgmuon1[i].pt	= cmgmuon1_pt[i];
      cmgmuon1[i].phi	= cmgmuon1_phi[i];
      cmgmuon1[i].eta	= cmgmuon1_eta[i];
    }

  cmgpfjet.resize(cmgpfjet_energy.size());
  for(unsigned int i=0; i < cmgpfjet.size(); ++i)
    {
      cmgpfjet[i].selected	= false;
      cmgpfjet[i].energy	= cmgpfjet_energy[i];
      cmgpfjet[i].et	= cmgpfjet_et[i];
      cmgpfjet[i].pt	= cmgpfjet_pt[i];
      cmgpfjet[i].phi	= cmgpfjet_phi[i];
      cmgpfjet[i].eta	= cmgpfjet_eta[i];
      cmgpfjet[i].jetArea	= cmgpfjet_jetArea[i];
      cmgpfjet[i].mass	= cmgpfjet_mass[i];
      cmgpfjet[i].rapidity	= cmgpfjet_rapidity[i];
      cmgpfjet[i].nConstituents	= cmgpfjet_nConstituents[i];
      cmgpfjet[i].partonFlavour	= cmgpfjet_partonFlavour[i];
      cmgpfjet[i].component_0_fraction	= cmgpfjet_component_0_fraction[i];
      cmgpfjet[i].component_0_number	= cmgpfjet_component_0_number[i];
      cmgpfjet[i].component_1_fraction	= cmgpfjet_component_1_fraction[i];
      cmgpfjet[i].component_1_number	= cmgpfjet_component_1_number[i];
      cmgpfjet[i].component_2_fraction	= cmgpfjet_component_2_fraction[i];
      cmgpfjet[i].component_2_number	= cmgpfjet_component_2_number[i];
      cmgpfjet[i].component_3_fraction	= cmgpfjet_component_3_fraction[i];
      cmgpfjet[i].component_3_number	= cmgpfjet_component_3_number[i];
      cmgpfjet[i].component_4_fraction	= cmgpfjet_component_4_fraction[i];
      cmgpfjet[i].component_4_number	= cmgpfjet_component_4_number[i];
      cmgpfjet[i].component_5_fraction	= cmgpfjet_component_5_fraction[i];
      cmgpfjet[i].component_5_number	= cmgpfjet_component_5_number[i];
      cmgpfjet[i].component_6_fraction	= cmgpfjet_component_6_fraction[i];
      cmgpfjet[i].component_6_number	= cmgpfjet_component_6_number[i];
      cmgpfjet[i].component_7_fraction	= cmgpfjet_component_7_fraction[i];
      cmgpfjet[i].component_7_number	= cmgpfjet_component_7_number[i];
      cmgpfjet[i].trackCountingHighEffBJetTag	= cmgpfjet_trackCountingHighEffBJetTag[i];
      cmgpfjet[i].trackCountingHighPurBJetTags	= cmgpfjet_trackCountingHighPurBJetTags[i];
      cmgpfjet[i].jetProbabilityBJetTags	= cmgpfjet_jetProbabilityBJetTags[i];
      cmgpfjet[i].jetBProbabilityBJetTags	= cmgpfjet_jetBProbabilityBJetTags[i];
      cmgpfjet[i].combinedSecondaryVertexBJetTags	= cmgpfjet_combinedSecondaryVertexBJetTags[i];
      cmgpfjet[i].combinedSecondaryVertexMVABJetTags	= cmgpfjet_combinedSecondaryVertexMVABJetTags[i];
    }

  cmgtau.resize(cmgtau_charge.size());
  for(unsigned int i=0; i < cmgtau.size(); ++i)
    {
      cmgtau[i].selected	= false;
      cmgtau[i].charge	= cmgtau_charge[i];
      cmgtau[i].energy	= cmgtau_energy[i];
      cmgtau[i].et	= cmgtau_et[i];
      cmgtau[i].pt	= cmgtau_pt[i];
      cmgtau[i].phi	= cmgtau_phi[i];
      cmgtau[i].eta	= cmgtau_eta[i];
      cmgtau[i].byLooseCombinedIsolationDeltaBetaCorr	= cmgtau_byLooseCombinedIsolationDeltaBetaCorr[i];
    }

  cmgtau1.resize(cmgtau1_charge.size());
  for(unsigned int i=0; i < cmgtau1.size(); ++i)
    {
      cmgtau1[i].selected	= false;
      cmgtau1[i].charge	= cmgtau1_charge[i];
      cmgtau1[i].energy	= cmgtau1_energy[i];
      cmgtau1[i].et	= cmgtau1_et[i];
      cmgtau1[i].pt	= cmgtau1_pt[i];
      cmgtau1[i].phi	= cmgtau1_phi[i];
      cmgtau1[i].eta	= cmgtau1_eta[i];
      cmgtau1[i].byTightCombinedIsolationDeltaBetaCorr	= cmgtau1_byTightCombinedIsolationDeltaBetaCorr[i];
    }

  electronhelper.resize(electronhelper_energy.size());
  for(unsigned int i=0; i < electronhelper.size(); ++i)
    {
      electronhelper[i].selected	= false;
      electronhelper[i].energy	= electronhelper_energy[i];
      electronhelper[i].et	= electronhelper_et[i];
      electronhelper[i].pt	= electronhelper_pt[i];
      electronhelper[i].phi	= electronhelper_phi[i];
      electronhelper[i].eta	= electronhelper_eta[i];
      electronhelper[i].charge	= electronhelper_charge[i];
      electronhelper[i].eSuperClusterOverP	= electronhelper_eSuperClusterOverP[i];
      electronhelper[i].deltaEtaSuperClusterTrackAtVtx	= electronhelper_deltaEtaSuperClusterTrackAtVtx[i];
      electronhelper[i].deltaPhiSuperClusterTrackAtVtx	= electronhelper_deltaPhiSuperClusterTrackAtVtx[i];
      electronhelper[i].isEB	= electronhelper_isEB[i];
      electronhelper[i].isEE	= electronhelper_isEE[i];
      electronhelper[i].scSigmaIEtaIEta	= electronhelper_scSigmaIEtaIEta[i];
      electronhelper[i].hadronicOverEm	= electronhelper_hadronicOverEm[i];
      electronhelper[i].superCluster_energy	= electronhelper_superCluster_energy[i];
      electronhelper[i].gsfTrack_trackerExpectedHitsInner_numberOfHits	= electronhelper_gsfTrack_trackerExpectedHitsInner_numberOfHits[i];
      electronhelper[i].simpleEleId80relIso	= electronhelper_simpleEleId80relIso[i];
      electronhelper[i].simpleEleId95relIso	= electronhelper_simpleEleId95relIso[i];
      electronhelper[i].trackIso	= electronhelper_trackIso[i];
      electronhelper[i].ecalIso	= electronhelper_ecalIso[i];
      electronhelper[i].hcalIso	= electronhelper_hcalIso[i];
      electronhelper[i].caloIso	= electronhelper_caloIso[i];
      electronhelper[i].dxywrtPV	= electronhelper_dxywrtPV[i];
      electronhelper[i].dzwrtPV	= electronhelper_dzwrtPV[i];
    }

  genjet.resize(genjet_energy.size());
  for(unsigned int i=0; i < genjet.size(); ++i)
    {
      genjet[i].selected	= false;
      genjet[i].energy	= genjet_energy[i];
      genjet[i].et	= genjet_et[i];
      genjet[i].pt	= genjet_pt[i];
      genjet[i].phi	= genjet_phi[i];
      genjet[i].eta	= genjet_eta[i];
      genjet[i].mass	= genjet_mass[i];
      genjet[i].rapidity	= genjet_rapidity[i];
      genjet[i].nConstituents	= genjet_nConstituents[i];
    }

  genparticlehelper.resize(genparticlehelper_firstMother.size());
  for(unsigned int i=0; i < genparticlehelper.size(); ++i)
    {
      genparticlehelper[i].selected	= false;
      genparticlehelper[i].firstMother	= genparticlehelper_firstMother[i];
      genparticlehelper[i].lastMother	= genparticlehelper_lastMother[i];
      genparticlehelper[i].firstDaughter	= genparticlehelper_firstDaughter[i];
      genparticlehelper[i].lastDaughter	= genparticlehelper_lastDaughter[i];
      genparticlehelper[i].charge	= genparticlehelper_charge[i];
      genparticlehelper[i].pdgId	= genparticlehelper_pdgId[i];
      genparticlehelper[i].status	= genparticlehelper_status[i];
      genparticlehelper[i].pt	= genparticlehelper_pt[i];
      genparticlehelper[i].eta	= genparticlehelper_eta[i];
      genparticlehelper[i].phi	= genparticlehelper_phi[i];
      genparticlehelper[i].mass	= genparticlehelper_mass[i];
    }

  jethelper.resize(jethelper_energy.size());
  for(unsigned int i=0; i < jethelper.size(); ++i)
    {
      jethelper[i].selected	= false;
      jethelper[i].energy	= jethelper_energy[i];
      jethelper[i].uncor_energy	= jethelper_uncor_energy[i];
      jethelper[i].et	= jethelper_et[i];
      jethelper[i].uncor_et	= jethelper_uncor_et[i];
      jethelper[i].pt	= jethelper_pt[i];
      jethelper[i].uncor_pt	= jethelper_uncor_pt[i];
      jethelper[i].phi	= jethelper_phi[i];
      jethelper[i].eta	= jethelper_eta[i];
      jethelper[i].rapidity	= jethelper_rapidity[i];
      jethelper[i].mass	= jethelper_mass[i];
      jethelper[i].jetArea	= jethelper_jetArea[i];
      jethelper[i].jetCharge03	= jethelper_jetCharge03[i];
      jethelper[i].jetCharge05	= jethelper_jetCharge05[i];
      jethelper[i].jetCharge10	= jethelper_jetCharge10[i];
      jethelper[i].chargedHadronEnergyFraction	= jethelper_chargedHadronEnergyFraction[i];
      jethelper[i].neutralHadronEnergyFraction	= jethelper_neutralHadronEnergyFraction[i];
      jethelper[i].chargedEmEnergyFraction	= jethelper_chargedEmEnergyFraction[i];
      jethelper[i].neutralEmEnergyFraction	= jethelper_neutralEmEnergyFraction[i];
      jethelper[i].photonEnergyFraction	= jethelper_photonEnergyFraction[i];
      jethelper[i].muonEnergyFraction	= jethelper_muonEnergyFraction[i];
      jethelper[i].chargedMultiplicity	= jethelper_chargedMultiplicity[i];
      jethelper[i].nConstituents	= jethelper_nConstituents[i];
      jethelper[i].partonFlavour	= jethelper_partonFlavour[i];
      jethelper[i].trackCountingHighEffBJetTags	= jethelper_trackCountingHighEffBJetTags[i];
      jethelper[i].trackCountingHighPurBJetTags	= jethelper_trackCountingHighPurBJetTags[i];
      jethelper[i].jetProbabilityBJetTags	= jethelper_jetProbabilityBJetTags[i];
      jethelper[i].jetBProbabilityBJetTags	= jethelper_jetBProbabilityBJetTags[i];
      jethelper[i].combinedSecondaryVertexBJetTags	= jethelper_combinedSecondaryVertexBJetTags[i];
      jethelper[i].combinedSecondaryVertexMVABJetTags	= jethelper_combinedSecondaryVertexMVABJetTags[i];
    }

  jethelper1.resize(jethelper1_energy.size());
  for(unsigned int i=0; i < jethelper1.size(); ++i)
    {
      jethelper1[i].selected	= false;
      jethelper1[i].energy	= jethelper1_energy[i];
      jethelper1[i].uncor_energy	= jethelper1_uncor_energy[i];
      jethelper1[i].et	= jethelper1_et[i];
      jethelper1[i].uncor_et	= jethelper1_uncor_et[i];
      jethelper1[i].pt	= jethelper1_pt[i];
      jethelper1[i].uncor_pt	= jethelper1_uncor_pt[i];
      jethelper1[i].phi	= jethelper1_phi[i];
      jethelper1[i].eta	= jethelper1_eta[i];
      jethelper1[i].rapidity	= jethelper1_rapidity[i];
      jethelper1[i].mass	= jethelper1_mass[i];
      jethelper1[i].jetArea	= jethelper1_jetArea[i];
      jethelper1[i].jetCharge03	= jethelper1_jetCharge03[i];
      jethelper1[i].jetCharge05	= jethelper1_jetCharge05[i];
      jethelper1[i].jetCharge10	= jethelper1_jetCharge10[i];
      jethelper1[i].chargedHadronEnergyFraction	= jethelper1_chargedHadronEnergyFraction[i];
      jethelper1[i].neutralHadronEnergyFraction	= jethelper1_neutralHadronEnergyFraction[i];
      jethelper1[i].chargedEmEnergyFraction	= jethelper1_chargedEmEnergyFraction[i];
      jethelper1[i].neutralEmEnergyFraction	= jethelper1_neutralEmEnergyFraction[i];
      jethelper1[i].photonEnergyFraction	= jethelper1_photonEnergyFraction[i];
      jethelper1[i].muonEnergyFraction	= jethelper1_muonEnergyFraction[i];
      jethelper1[i].chargedMultiplicity	= jethelper1_chargedMultiplicity[i];
      jethelper1[i].nConstituents	= jethelper1_nConstituents[i];
      jethelper1[i].partonFlavour	= jethelper1_partonFlavour[i];
      jethelper1[i].trackCountingHighEffBJetTags	= jethelper1_trackCountingHighEffBJetTags[i];
      jethelper1[i].trackCountingHighPurBJetTags	= jethelper1_trackCountingHighPurBJetTags[i];
      jethelper1[i].jetProbabilityBJetTags	= jethelper1_jetProbabilityBJetTags[i];
      jethelper1[i].jetBProbabilityBJetTags	= jethelper1_jetBProbabilityBJetTags[i];
      jethelper1[i].combinedSecondaryVertexBJetTags	= jethelper1_combinedSecondaryVertexBJetTags[i];
      jethelper1[i].combinedSecondaryVertexMVABJetTags	= jethelper1_combinedSecondaryVertexMVABJetTags[i];
    }

  jethelper2.resize(jethelper2_energy.size());
  for(unsigned int i=0; i < jethelper2.size(); ++i)
    {
      jethelper2[i].selected	= false;
      jethelper2[i].energy	= jethelper2_energy[i];
      jethelper2[i].uncor_energy	= jethelper2_uncor_energy[i];
      jethelper2[i].et	= jethelper2_et[i];
      jethelper2[i].uncor_et	= jethelper2_uncor_et[i];
      jethelper2[i].pt	= jethelper2_pt[i];
      jethelper2[i].uncor_pt	= jethelper2_uncor_pt[i];
      jethelper2[i].phi	= jethelper2_phi[i];
      jethelper2[i].eta	= jethelper2_eta[i];
      jethelper2[i].rapidity	= jethelper2_rapidity[i];
      jethelper2[i].mass	= jethelper2_mass[i];
      jethelper2[i].jetArea	= jethelper2_jetArea[i];
      jethelper2[i].jetCharge03	= jethelper2_jetCharge03[i];
      jethelper2[i].jetCharge05	= jethelper2_jetCharge05[i];
      jethelper2[i].jetCharge10	= jethelper2_jetCharge10[i];
      jethelper2[i].nConstituents	= jethelper2_nConstituents[i];
      jethelper2[i].partonFlavour	= jethelper2_partonFlavour[i];
      jethelper2[i].trackCountingHighEffBJetTags	= jethelper2_trackCountingHighEffBJetTags[i];
      jethelper2[i].trackCountingHighPurBJetTags	= jethelper2_trackCountingHighPurBJetTags[i];
      jethelper2[i].jetProbabilityBJetTags	= jethelper2_jetProbabilityBJetTags[i];
      jethelper2[i].jetBProbabilityBJetTags	= jethelper2_jetBProbabilityBJetTags[i];
      jethelper2[i].combinedSecondaryVertexBJetTags	= jethelper2_combinedSecondaryVertexBJetTags[i];
      jethelper2[i].combinedSecondaryVertexMVABJetTags	= jethelper2_combinedSecondaryVertexMVABJetTags[i];
      jethelper2[i].numberOfDaughters	= jethelper2_numberOfDaughters[i];
      jethelper2[i].daughter_0_energy	= jethelper2_daughter_0_energy[i];
      jethelper2[i].daughter_0_pt	= jethelper2_daughter_0_pt[i];
      jethelper2[i].daughter_0_eta	= jethelper2_daughter_0_eta[i];
      jethelper2[i].daughter_0_rapidity	= jethelper2_daughter_0_rapidity[i];
      jethelper2[i].daughter_0_phi	= jethelper2_daughter_0_phi[i];
      jethelper2[i].daughter_0_mass	= jethelper2_daughter_0_mass[i];
      jethelper2[i].daughter_1_energy	= jethelper2_daughter_1_energy[i];
      jethelper2[i].daughter_1_pt	= jethelper2_daughter_1_pt[i];
      jethelper2[i].daughter_1_eta	= jethelper2_daughter_1_eta[i];
      jethelper2[i].daughter_1_rapidity	= jethelper2_daughter_1_rapidity[i];
      jethelper2[i].daughter_1_phi	= jethelper2_daughter_1_phi[i];
      jethelper2[i].daughter_1_mass	= jethelper2_daughter_1_mass[i];
      jethelper2[i].daughter_0_jetCharge03	= jethelper2_daughter_0_jetCharge03[i];
      jethelper2[i].daughter_0_jetCharge05	= jethelper2_daughter_0_jetCharge05[i];
      jethelper2[i].daughter_0_jetCharge10	= jethelper2_daughter_0_jetCharge10[i];
      jethelper2[i].daughter_1_jetCharge03	= jethelper2_daughter_1_jetCharge03[i];
      jethelper2[i].daughter_1_jetCharge05	= jethelper2_daughter_1_jetCharge05[i];
      jethelper2[i].daughter_1_jetCharge10	= jethelper2_daughter_1_jetCharge10[i];
    }

  jethelper3.resize(jethelper3_energy.size());
  for(unsigned int i=0; i < jethelper3.size(); ++i)
    {
      jethelper3[i].selected	= false;
      jethelper3[i].energy	= jethelper3_energy[i];
      jethelper3[i].uncor_energy	= jethelper3_uncor_energy[i];
      jethelper3[i].et	= jethelper3_et[i];
      jethelper3[i].uncor_et	= jethelper3_uncor_et[i];
      jethelper3[i].pt	= jethelper3_pt[i];
      jethelper3[i].uncor_pt	= jethelper3_uncor_pt[i];
      jethelper3[i].phi	= jethelper3_phi[i];
      jethelper3[i].eta	= jethelper3_eta[i];
      jethelper3[i].rapidity	= jethelper3_rapidity[i];
      jethelper3[i].mass	= jethelper3_mass[i];
      jethelper3[i].jetArea	= jethelper3_jetArea[i];
      jethelper3[i].jetCharge03	= jethelper3_jetCharge03[i];
      jethelper3[i].jetCharge05	= jethelper3_jetCharge05[i];
      jethelper3[i].jetCharge10	= jethelper3_jetCharge10[i];
      jethelper3[i].chargedHadronEnergyFraction	= jethelper3_chargedHadronEnergyFraction[i];
      jethelper3[i].neutralHadronEnergyFraction	= jethelper3_neutralHadronEnergyFraction[i];
      jethelper3[i].chargedEmEnergyFraction	= jethelper3_chargedEmEnergyFraction[i];
      jethelper3[i].neutralEmEnergyFraction	= jethelper3_neutralEmEnergyFraction[i];
      jethelper3[i].photonEnergyFraction	= jethelper3_photonEnergyFraction[i];
      jethelper3[i].muonEnergyFraction	= jethelper3_muonEnergyFraction[i];
      jethelper3[i].chargedMultiplicity	= jethelper3_chargedMultiplicity[i];
      jethelper3[i].nConstituents	= jethelper3_nConstituents[i];
      jethelper3[i].partonFlavour	= jethelper3_partonFlavour[i];
      jethelper3[i].trackCountingHighEffBJetTags	= jethelper3_trackCountingHighEffBJetTags[i];
      jethelper3[i].trackCountingHighPurBJetTags	= jethelper3_trackCountingHighPurBJetTags[i];
      jethelper3[i].jetProbabilityBJetTags	= jethelper3_jetProbabilityBJetTags[i];
      jethelper3[i].jetBProbabilityBJetTags	= jethelper3_jetBProbabilityBJetTags[i];
      jethelper3[i].combinedSecondaryVertexBJetTags	= jethelper3_combinedSecondaryVertexBJetTags[i];
      jethelper3[i].combinedSecondaryVertexMVABJetTags	= jethelper3_combinedSecondaryVertexMVABJetTags[i];
      jethelper3[i].qjetsvolatility	= jethelper3_qjetsvolatility[i];
      jethelper3[i].tau1	= jethelper3_tau1[i];
      jethelper3[i].tau2	= jethelper3_tau2[i];
      jethelper3[i].tau3	= jethelper3_tau3[i];
      jethelper3[i].C2beta17	= jethelper3_C2beta17[i];
    }

  jethelper4.resize(jethelper4_energy.size());
  for(unsigned int i=0; i < jethelper4.size(); ++i)
    {
      jethelper4[i].selected	= false;
      jethelper4[i].energy	= jethelper4_energy[i];
      jethelper4[i].uncor_energy	= jethelper4_uncor_energy[i];
      jethelper4[i].et	= jethelper4_et[i];
      jethelper4[i].uncor_et	= jethelper4_uncor_et[i];
      jethelper4[i].pt	= jethelper4_pt[i];
      jethelper4[i].uncor_pt	= jethelper4_uncor_pt[i];
      jethelper4[i].phi	= jethelper4_phi[i];
      jethelper4[i].eta	= jethelper4_eta[i];
      jethelper4[i].rapidity	= jethelper4_rapidity[i];
      jethelper4[i].mass	= jethelper4_mass[i];
      jethelper4[i].jetArea	= jethelper4_jetArea[i];
      jethelper4[i].jetCharge03	= jethelper4_jetCharge03[i];
      jethelper4[i].jetCharge05	= jethelper4_jetCharge05[i];
      jethelper4[i].jetCharge10	= jethelper4_jetCharge10[i];
      jethelper4[i].nConstituents	= jethelper4_nConstituents[i];
      jethelper4[i].partonFlavour	= jethelper4_partonFlavour[i];
      jethelper4[i].trackCountingHighEffBJetTags	= jethelper4_trackCountingHighEffBJetTags[i];
      jethelper4[i].trackCountingHighPurBJetTags	= jethelper4_trackCountingHighPurBJetTags[i];
      jethelper4[i].jetProbabilityBJetTags	= jethelper4_jetProbabilityBJetTags[i];
      jethelper4[i].jetBProbabilityBJetTags	= jethelper4_jetBProbabilityBJetTags[i];
      jethelper4[i].combinedSecondaryVertexBJetTags	= jethelper4_combinedSecondaryVertexBJetTags[i];
      jethelper4[i].combinedSecondaryVertexMVABJetTags	= jethelper4_combinedSecondaryVertexMVABJetTags[i];
      jethelper4[i].numberOfDaughters	= jethelper4_numberOfDaughters[i];
      jethelper4[i].daughter_0_energy	= jethelper4_daughter_0_energy[i];
      jethelper4[i].daughter_0_pt	= jethelper4_daughter_0_pt[i];
      jethelper4[i].daughter_0_eta	= jethelper4_daughter_0_eta[i];
      jethelper4[i].daughter_0_rapidity	= jethelper4_daughter_0_rapidity[i];
      jethelper4[i].daughter_0_phi	= jethelper4_daughter_0_phi[i];
      jethelper4[i].daughter_0_mass	= jethelper4_daughter_0_mass[i];
      jethelper4[i].daughter_1_energy	= jethelper4_daughter_1_energy[i];
      jethelper4[i].daughter_1_pt	= jethelper4_daughter_1_pt[i];
      jethelper4[i].daughter_1_eta	= jethelper4_daughter_1_eta[i];
      jethelper4[i].daughter_1_rapidity	= jethelper4_daughter_1_rapidity[i];
      jethelper4[i].daughter_1_phi	= jethelper4_daughter_1_phi[i];
      jethelper4[i].daughter_1_mass	= jethelper4_daughter_1_mass[i];
      jethelper4[i].tau1	= jethelper4_tau1[i];
      jethelper4[i].tau2	= jethelper4_tau2[i];
      jethelper4[i].tau3	= jethelper4_tau3[i];
      jethelper4[i].daughter_0_jetCharge03	= jethelper4_daughter_0_jetCharge03[i];
      jethelper4[i].daughter_0_jetCharge05	= jethelper4_daughter_0_jetCharge05[i];
      jethelper4[i].daughter_0_jetCharge10	= jethelper4_daughter_0_jetCharge10[i];
      jethelper4[i].daughter_1_jetCharge03	= jethelper4_daughter_1_jetCharge03[i];
      jethelper4[i].daughter_1_jetCharge05	= jethelper4_daughter_1_jetCharge05[i];
      jethelper4[i].daughter_1_jetCharge10	= jethelper4_daughter_1_jetCharge10[i];
    }

  jethelper5.resize(jethelper5_energy.size());
  for(unsigned int i=0; i < jethelper5.size(); ++i)
    {
      jethelper5[i].selected	= false;
      jethelper5[i].energy	= jethelper5_energy[i];
      jethelper5[i].uncor_energy	= jethelper5_uncor_energy[i];
      jethelper5[i].et	= jethelper5_et[i];
      jethelper5[i].uncor_et	= jethelper5_uncor_et[i];
      jethelper5[i].pt	= jethelper5_pt[i];
      jethelper5[i].uncor_pt	= jethelper5_uncor_pt[i];
      jethelper5[i].phi	= jethelper5_phi[i];
      jethelper5[i].eta	= jethelper5_eta[i];
      jethelper5[i].rapidity	= jethelper5_rapidity[i];
      jethelper5[i].mass	= jethelper5_mass[i];
      jethelper5[i].jetArea	= jethelper5_jetArea[i];
      jethelper5[i].jetCharge03	= jethelper5_jetCharge03[i];
      jethelper5[i].jetCharge05	= jethelper5_jetCharge05[i];
      jethelper5[i].jetCharge10	= jethelper5_jetCharge10[i];
      jethelper5[i].chargedHadronEnergyFraction	= jethelper5_chargedHadronEnergyFraction[i];
      jethelper5[i].neutralHadronEnergyFraction	= jethelper5_neutralHadronEnergyFraction[i];
      jethelper5[i].chargedEmEnergyFraction	= jethelper5_chargedEmEnergyFraction[i];
      jethelper5[i].neutralEmEnergyFraction	= jethelper5_neutralEmEnergyFraction[i];
      jethelper5[i].photonEnergyFraction	= jethelper5_photonEnergyFraction[i];
      jethelper5[i].muonEnergyFraction	= jethelper5_muonEnergyFraction[i];
      jethelper5[i].chargedMultiplicity	= jethelper5_chargedMultiplicity[i];
      jethelper5[i].nConstituents	= jethelper5_nConstituents[i];
      jethelper5[i].partonFlavour	= jethelper5_partonFlavour[i];
      jethelper5[i].trackCountingHighEffBJetTags	= jethelper5_trackCountingHighEffBJetTags[i];
      jethelper5[i].trackCountingHighPurBJetTags	= jethelper5_trackCountingHighPurBJetTags[i];
      jethelper5[i].jetProbabilityBJetTags	= jethelper5_jetProbabilityBJetTags[i];
      jethelper5[i].jetBProbabilityBJetTags	= jethelper5_jetBProbabilityBJetTags[i];
      jethelper5[i].combinedSecondaryVertexBJetTags	= jethelper5_combinedSecondaryVertexBJetTags[i];
      jethelper5[i].combinedSecondaryVertexMVABJetTags	= jethelper5_combinedSecondaryVertexMVABJetTags[i];
      jethelper5[i].qjetsvolatility	= jethelper5_qjetsvolatility[i];
      jethelper5[i].tau1	= jethelper5_tau1[i];
      jethelper5[i].tau2	= jethelper5_tau2[i];
      jethelper5[i].tau3	= jethelper5_tau3[i];
      jethelper5[i].C2beta17	= jethelper5_C2beta17[i];
      jethelper5[i].getNcharged01	= jethelper5_getNcharged01[i];
      jethelper5[i].getNneutral01	= jethelper5_getNneutral01[i];
      jethelper5[i].getChargedPt0	= jethelper5_getChargedPt0[i];
      jethelper5[i].getChargedPt1	= jethelper5_getChargedPt1[i];
      jethelper5[i].getChargedPt2	= jethelper5_getChargedPt2[i];
      jethelper5[i].getChargedPt3	= jethelper5_getChargedPt3[i];
      jethelper5[i].getPt0	= jethelper5_getPt0[i];
      jethelper5[i].getPt1	= jethelper5_getPt1[i];
      jethelper5[i].getPt2	= jethelper5_getPt2[i];
      jethelper5[i].getPt3	= jethelper5_getPt3[i];
    }

  jethelper6.resize(jethelper6_energy.size());
  for(unsigned int i=0; i < jethelper6.size(); ++i)
    {
      jethelper6[i].selected	= false;
      jethelper6[i].energy	= jethelper6_energy[i];
      jethelper6[i].et	= jethelper6_et[i];
      jethelper6[i].pt	= jethelper6_pt[i];
      jethelper6[i].phi	= jethelper6_phi[i];
      jethelper6[i].eta	= jethelper6_eta[i];
      jethelper6[i].rapidity	= jethelper6_rapidity[i];
      jethelper6[i].mass	= jethelper6_mass[i];
      jethelper6[i].jetArea	= jethelper6_jetArea[i];
      jethelper6[i].genJetCharge03	= jethelper6_genJetCharge03[i];
      jethelper6[i].genJetCharge05	= jethelper6_genJetCharge05[i];
      jethelper6[i].genJetCharge10	= jethelper6_genJetCharge10[i];
      jethelper6[i].nConstituents	= jethelper6_nConstituents[i];
      jethelper6[i].partonFlavour	= jethelper6_partonFlavour[i];
      jethelper6[i].genC2beta17	= jethelper6_genC2beta17[i];
      jethelper6[i].genC2beta17CHS	= jethelper6_genC2beta17CHS[i];
      jethelper6[i].genTau1	= jethelper6_genTau1[i];
      jethelper6[i].genTau2	= jethelper6_genTau2[i];
      jethelper6[i].genTau3	= jethelper6_genTau3[i];
      jethelper6[i].genTau1Pt2	= jethelper6_genTau1Pt2[i];
      jethelper6[i].genTau2Pt2	= jethelper6_genTau2Pt2[i];
      jethelper6[i].genTau1Pt5	= jethelper6_genTau1Pt5[i];
      jethelper6[i].genTau2Pt5	= jethelper6_genTau2Pt5[i];
      jethelper6[i].genTau1CHS	= jethelper6_genTau1CHS[i];
      jethelper6[i].genTau2CHS	= jethelper6_genTau2CHS[i];
      jethelper6[i].genNCHS	= jethelper6_genNCHS[i];
    }

  jethelper7.resize(jethelper7_energy.size());
  for(unsigned int i=0; i < jethelper7.size(); ++i)
    {
      jethelper7[i].selected	= false;
      jethelper7[i].energy	= jethelper7_energy[i];
      jethelper7[i].et	= jethelper7_et[i];
      jethelper7[i].pt	= jethelper7_pt[i];
      jethelper7[i].phi	= jethelper7_phi[i];
      jethelper7[i].eta	= jethelper7_eta[i];
      jethelper7[i].rapidity	= jethelper7_rapidity[i];
      jethelper7[i].mass	= jethelper7_mass[i];
      jethelper7[i].jetArea	= jethelper7_jetArea[i];
      jethelper7[i].jetCharge03	= jethelper7_jetCharge03[i];
      jethelper7[i].jetCharge05	= jethelper7_jetCharge05[i];
      jethelper7[i].jetCharge10	= jethelper7_jetCharge10[i];
      jethelper7[i].nConstituents	= jethelper7_nConstituents[i];
      jethelper7[i].partonFlavour	= jethelper7_partonFlavour[i];
      jethelper7[i].numberOfDaughters	= jethelper7_numberOfDaughters[i];
      jethelper7[i].daughter_0_energy	= jethelper7_daughter_0_energy[i];
      jethelper7[i].daughter_0_pt	= jethelper7_daughter_0_pt[i];
      jethelper7[i].daughter_0_eta	= jethelper7_daughter_0_eta[i];
      jethelper7[i].daughter_0_rapidity	= jethelper7_daughter_0_rapidity[i];
      jethelper7[i].daughter_0_phi	= jethelper7_daughter_0_phi[i];
      jethelper7[i].daughter_0_mass	= jethelper7_daughter_0_mass[i];
      jethelper7[i].daughter_1_energy	= jethelper7_daughter_1_energy[i];
      jethelper7[i].daughter_1_pt	= jethelper7_daughter_1_pt[i];
      jethelper7[i].daughter_1_eta	= jethelper7_daughter_1_eta[i];
      jethelper7[i].daughter_1_rapidity	= jethelper7_daughter_1_rapidity[i];
      jethelper7[i].daughter_1_phi	= jethelper7_daughter_1_phi[i];
      jethelper7[i].daughter_1_mass	= jethelper7_daughter_1_mass[i];
      jethelper7[i].genTau1	= jethelper7_genTau1[i];
      jethelper7[i].genTau2	= jethelper7_genTau2[i];
      jethelper7[i].genTau3	= jethelper7_genTau3[i];
    }

  leafcandidate.resize(leafcandidate_charge.size());
  for(unsigned int i=0; i < leafcandidate.size(); ++i)
    {
      leafcandidate[i].selected	= false;
      leafcandidate[i].charge	= leafcandidate_charge[i];
      leafcandidate[i].energy	= leafcandidate_energy[i];
      leafcandidate[i].et	= leafcandidate_et[i];
      leafcandidate[i].pt	= leafcandidate_pt[i];
      leafcandidate[i].phi	= leafcandidate_phi[i];
      leafcandidate[i].eta	= leafcandidate_eta[i];
    }

  met.resize(met_energy.size());
  for(unsigned int i=0; i < met.size(); ++i)
    {
      met[i].selected	= false;
      met[i].energy	= met_energy[i];
      met[i].et	= met_et[i];
      met[i].pt	= met_pt[i];
      met[i].phi	= met_phi[i];
      met[i].sumEt	= met_sumEt[i];
      met[i].mEtSig	= met_mEtSig[i];
      met[i].significance	= met_significance[i];
    }

  met1.resize(met1_energy.size());
  for(unsigned int i=0; i < met1.size(); ++i)
    {
      met1[i].selected	= false;
      met1[i].energy	= met1_energy[i];
      met1[i].et	= met1_et[i];
      met1[i].pt	= met1_pt[i];
      met1[i].phi	= met1_phi[i];
      met1[i].sumEt	= met1_sumEt[i];
      met1[i].mEtSig	= met1_mEtSig[i];
      met1[i].significance	= met1_significance[i];
    }

  muonhelper.resize(muonhelper_energy.size());
  for(unsigned int i=0; i < muonhelper.size(); ++i)
    {
      muonhelper[i].selected	= false;
      muonhelper[i].energy	= muonhelper_energy[i];
      muonhelper[i].et	= muonhelper_et[i];
      muonhelper[i].pt	= muonhelper_pt[i];
      muonhelper[i].phi	= muonhelper_phi[i];
      muonhelper[i].eta	= muonhelper_eta[i];
      muonhelper[i].charge	= muonhelper_charge[i];
      muonhelper[i].TMOneStationTight	= muonhelper_TMOneStationTight[i];
      muonhelper[i].isGlobalMuon	= muonhelper_isGlobalMuon[i];
      muonhelper[i].isTrackerMuon	= muonhelper_isTrackerMuon[i];
      muonhelper[i].isPFMuon	= muonhelper_isPFMuon[i];
      muonhelper[i].track_hitPattern_trackerLayersWithMeasurement	= muonhelper_track_hitPattern_trackerLayersWithMeasurement[i];
      muonhelper[i].innerTrack_hitPattern_pixelLayersWithMeasurement	= muonhelper_innerTrack_hitPattern_pixelLayersWithMeasurement[i];
      muonhelper[i].innerTrack_normalizedChi2	= muonhelper_innerTrack_normalizedChi2[i];
      muonhelper[i].globalTrack_normalizedChi2	= muonhelper_globalTrack_normalizedChi2[i];
      muonhelper[i].globalTrack_hitPattern_numberOfValidMuonHits	= muonhelper_globalTrack_hitPattern_numberOfValidMuonHits[i];
      muonhelper[i].numberOfMatchedStations	= muonhelper_numberOfMatchedStations[i];
      muonhelper[i].dB	= muonhelper_dB[i];
      muonhelper[i].innerTrack_hitPattern_numberOfValidPixelHits	= muonhelper_innerTrack_hitPattern_numberOfValidPixelHits[i];
      muonhelper[i].pfIsolationR04_sumChargedHadronPt	= muonhelper_pfIsolationR04_sumChargedHadronPt[i];
      muonhelper[i].pfIsolationR04_sumChargedParticlePt	= muonhelper_pfIsolationR04_sumChargedParticlePt[i];
      muonhelper[i].pfIsolationR04_sumNeutralHadronEt	= muonhelper_pfIsolationR04_sumNeutralHadronEt[i];
      muonhelper[i].pfIsolationR04_sumPhotonEt	= muonhelper_pfIsolationR04_sumPhotonEt[i];
      muonhelper[i].pfIsolationR04_sumNeutralHadronEtHighThreshold	= muonhelper_pfIsolationR04_sumNeutralHadronEtHighThreshold[i];
      muonhelper[i].pfIsolationR04_sumPhotonEtHighThreshold	= muonhelper_pfIsolationR04_sumPhotonEtHighThreshold[i];
      muonhelper[i].pfIsolationR04_sumPUPt	= muonhelper_pfIsolationR04_sumPUPt[i];
      muonhelper[i].dxywrtPV	= muonhelper_dxywrtPV[i];
      muonhelper[i].dzwrtPV	= muonhelper_dzwrtPV[i];
    }

  pileupsummaryinfo.resize(pileupsummaryinfo_getPU_NumInteractions.size());
  for(unsigned int i=0; i < pileupsummaryinfo.size(); ++i)
    {
      pileupsummaryinfo[i].selected	= false;
      pileupsummaryinfo[i].getPU_NumInteractions	= pileupsummaryinfo_getPU_NumInteractions[i];
      pileupsummaryinfo[i].getBunchCrossing	= pileupsummaryinfo_getBunchCrossing[i];
      pileupsummaryinfo[i].getTrueNumInteractions	= pileupsummaryinfo_getTrueNumInteractions[i];
    }

  tau.resize(tau_energy.size());
  for(unsigned int i=0; i < tau.size(); ++i)
    {
      tau[i].selected	= false;
      tau[i].energy	= tau_energy[i];
      tau[i].et	= tau_et[i];
      tau[i].pt	= tau_pt[i];
      tau[i].phi	= tau_phi[i];
      tau[i].eta	= tau_eta[i];
      tau[i].byLooseCombinedIsolationDeltaBetaCorr	= tau_byLooseCombinedIsolationDeltaBetaCorr[i];
      tau[i].byMediumCombinedIsolationDeltaBetaCorr	= tau_byMediumCombinedIsolationDeltaBetaCorr[i];
      tau[i].trackIso	= tau_trackIso[i];
      tau[i].ecalIso	= tau_ecalIso[i];
      tau[i].hcalIso	= tau_hcalIso[i];
      tau[i].caloIso	= tau_caloIso[i];
    }

  vertex.resize(vertex_isFake.size());
  for(unsigned int i=0; i < vertex.size(); ++i)
    {
      vertex[i].selected	= false;
      vertex[i].isFake	= vertex_isFake[i];
      vertex[i].chi2	= vertex_chi2[i];
      vertex[i].ndof	= vertex_ndof[i];
      vertex[i].x	= vertex_x[i];
      vertex[i].y	= vertex_y[i];
      vertex[i].z	= vertex_z[i];
      vertex[i].xError	= vertex_xError[i];
      vertex[i].yError	= vertex_yError[i];
      vertex[i].zError	= vertex_zError[i];
      vertex[i].position_Rho	= vertex_position_Rho[i];
    }
}

//-----------------------------------------------------------------------------
// --- Call saveSelectedObjects() just before call to addEvent if
// --- you wish to save only the selected objects
//-----------------------------------------------------------------------------
void saveSelectedObjects()
{
  if ( ! fillObjectsCalled ) return;
  int n = 0;

  n = 0;
  for(unsigned int i=0; i < calomet.size(); ++i)
    {
      if ( ! calomet[i].selected ) continue;
      calomet_energy[n]	= calomet[i].energy;
      calomet_pt[n]	= calomet[i].pt;
      calomet_phi[n]	= calomet[i].phi;
      calomet_sumEt[n]	= calomet[i].sumEt;
      calomet_mEtSig[n]	= calomet[i].mEtSig;
      calomet_significance[n]	= calomet[i].significance;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < calomet1.size(); ++i)
    {
      if ( ! calomet1[i].selected ) continue;
      calomet1_energy[n]	= calomet1[i].energy;
      calomet1_pt[n]	= calomet1[i].pt;
      calomet1_phi[n]	= calomet1[i].phi;
      calomet1_sumEt[n]	= calomet1[i].sumEt;
      calomet1_mEtSig[n]	= calomet1[i].mEtSig;
      calomet1_significance[n]	= calomet1[i].significance;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgbasemet.size(); ++i)
    {
      if ( ! cmgbasemet[i].selected ) continue;
      cmgbasemet_energy[n]	= cmgbasemet[i].energy;
      cmgbasemet_et[n]	= cmgbasemet[i].et;
      cmgbasemet_pt[n]	= cmgbasemet[i].pt;
      cmgbasemet_phi[n]	= cmgbasemet[i].phi;
      cmgbasemet_sumEt[n]	= cmgbasemet[i].sumEt;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgbasemet1.size(); ++i)
    {
      if ( ! cmgbasemet1[i].selected ) continue;
      cmgbasemet1_energy[n]	= cmgbasemet1[i].energy;
      cmgbasemet1_et[n]	= cmgbasemet1[i].et;
      cmgbasemet1_pt[n]	= cmgbasemet1[i].pt;
      cmgbasemet1_phi[n]	= cmgbasemet1[i].phi;
      cmgbasemet1_sumEt[n]	= cmgbasemet1[i].sumEt;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgbasemet2.size(); ++i)
    {
      if ( ! cmgbasemet2[i].selected ) continue;
      cmgbasemet2_energy[n]	= cmgbasemet2[i].energy;
      cmgbasemet2_et[n]	= cmgbasemet2[i].et;
      cmgbasemet2_pt[n]	= cmgbasemet2[i].pt;
      cmgbasemet2_phi[n]	= cmgbasemet2[i].phi;
      cmgbasemet2_sumEt[n]	= cmgbasemet2[i].sumEt;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgelectron.size(); ++i)
    {
      if ( ! cmgelectron[i].selected ) continue;
      cmgelectron_charge[n]	= cmgelectron[i].charge;
      cmgelectron_energy[n]	= cmgelectron[i].energy;
      cmgelectron_et[n]	= cmgelectron[i].et;
      cmgelectron_pt[n]	= cmgelectron[i].pt;
      cmgelectron_phi[n]	= cmgelectron[i].phi;
      cmgelectron_eta[n]	= cmgelectron[i].eta;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgelectron1.size(); ++i)
    {
      if ( ! cmgelectron1[i].selected ) continue;
      cmgelectron1_charge[n]	= cmgelectron1[i].charge;
      cmgelectron1_energy[n]	= cmgelectron1[i].energy;
      cmgelectron1_et[n]	= cmgelectron1[i].et;
      cmgelectron1_pt[n]	= cmgelectron1[i].pt;
      cmgelectron1_phi[n]	= cmgelectron1[i].phi;
      cmgelectron1_eta[n]	= cmgelectron1[i].eta;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgmuon.size(); ++i)
    {
      if ( ! cmgmuon[i].selected ) continue;
      cmgmuon_charge[n]	= cmgmuon[i].charge;
      cmgmuon_energy[n]	= cmgmuon[i].energy;
      cmgmuon_et[n]	= cmgmuon[i].et;
      cmgmuon_pt[n]	= cmgmuon[i].pt;
      cmgmuon_phi[n]	= cmgmuon[i].phi;
      cmgmuon_eta[n]	= cmgmuon[i].eta;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgmuon1.size(); ++i)
    {
      if ( ! cmgmuon1[i].selected ) continue;
      cmgmuon1_charge[n]	= cmgmuon1[i].charge;
      cmgmuon1_energy[n]	= cmgmuon1[i].energy;
      cmgmuon1_et[n]	= cmgmuon1[i].et;
      cmgmuon1_pt[n]	= cmgmuon1[i].pt;
      cmgmuon1_phi[n]	= cmgmuon1[i].phi;
      cmgmuon1_eta[n]	= cmgmuon1[i].eta;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgpfjet.size(); ++i)
    {
      if ( ! cmgpfjet[i].selected ) continue;
      cmgpfjet_energy[n]	= cmgpfjet[i].energy;
      cmgpfjet_et[n]	= cmgpfjet[i].et;
      cmgpfjet_pt[n]	= cmgpfjet[i].pt;
      cmgpfjet_phi[n]	= cmgpfjet[i].phi;
      cmgpfjet_eta[n]	= cmgpfjet[i].eta;
      cmgpfjet_jetArea[n]	= cmgpfjet[i].jetArea;
      cmgpfjet_mass[n]	= cmgpfjet[i].mass;
      cmgpfjet_rapidity[n]	= cmgpfjet[i].rapidity;
      cmgpfjet_nConstituents[n]	= cmgpfjet[i].nConstituents;
      cmgpfjet_partonFlavour[n]	= cmgpfjet[i].partonFlavour;
      cmgpfjet_component_0_fraction[n]	= cmgpfjet[i].component_0_fraction;
      cmgpfjet_component_0_number[n]	= cmgpfjet[i].component_0_number;
      cmgpfjet_component_1_fraction[n]	= cmgpfjet[i].component_1_fraction;
      cmgpfjet_component_1_number[n]	= cmgpfjet[i].component_1_number;
      cmgpfjet_component_2_fraction[n]	= cmgpfjet[i].component_2_fraction;
      cmgpfjet_component_2_number[n]	= cmgpfjet[i].component_2_number;
      cmgpfjet_component_3_fraction[n]	= cmgpfjet[i].component_3_fraction;
      cmgpfjet_component_3_number[n]	= cmgpfjet[i].component_3_number;
      cmgpfjet_component_4_fraction[n]	= cmgpfjet[i].component_4_fraction;
      cmgpfjet_component_4_number[n]	= cmgpfjet[i].component_4_number;
      cmgpfjet_component_5_fraction[n]	= cmgpfjet[i].component_5_fraction;
      cmgpfjet_component_5_number[n]	= cmgpfjet[i].component_5_number;
      cmgpfjet_component_6_fraction[n]	= cmgpfjet[i].component_6_fraction;
      cmgpfjet_component_6_number[n]	= cmgpfjet[i].component_6_number;
      cmgpfjet_component_7_fraction[n]	= cmgpfjet[i].component_7_fraction;
      cmgpfjet_component_7_number[n]	= cmgpfjet[i].component_7_number;
      cmgpfjet_trackCountingHighEffBJetTag[n]	= cmgpfjet[i].trackCountingHighEffBJetTag;
      cmgpfjet_trackCountingHighPurBJetTags[n]	= cmgpfjet[i].trackCountingHighPurBJetTags;
      cmgpfjet_jetProbabilityBJetTags[n]	= cmgpfjet[i].jetProbabilityBJetTags;
      cmgpfjet_jetBProbabilityBJetTags[n]	= cmgpfjet[i].jetBProbabilityBJetTags;
      cmgpfjet_combinedSecondaryVertexBJetTags[n]	= cmgpfjet[i].combinedSecondaryVertexBJetTags;
      cmgpfjet_combinedSecondaryVertexMVABJetTags[n]	= cmgpfjet[i].combinedSecondaryVertexMVABJetTags;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgtau.size(); ++i)
    {
      if ( ! cmgtau[i].selected ) continue;
      cmgtau_charge[n]	= cmgtau[i].charge;
      cmgtau_energy[n]	= cmgtau[i].energy;
      cmgtau_et[n]	= cmgtau[i].et;
      cmgtau_pt[n]	= cmgtau[i].pt;
      cmgtau_phi[n]	= cmgtau[i].phi;
      cmgtau_eta[n]	= cmgtau[i].eta;
      cmgtau_byLooseCombinedIsolationDeltaBetaCorr[n]	= cmgtau[i].byLooseCombinedIsolationDeltaBetaCorr;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < cmgtau1.size(); ++i)
    {
      if ( ! cmgtau1[i].selected ) continue;
      cmgtau1_charge[n]	= cmgtau1[i].charge;
      cmgtau1_energy[n]	= cmgtau1[i].energy;
      cmgtau1_et[n]	= cmgtau1[i].et;
      cmgtau1_pt[n]	= cmgtau1[i].pt;
      cmgtau1_phi[n]	= cmgtau1[i].phi;
      cmgtau1_eta[n]	= cmgtau1[i].eta;
      cmgtau1_byTightCombinedIsolationDeltaBetaCorr[n]	= cmgtau1[i].byTightCombinedIsolationDeltaBetaCorr;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < electronhelper.size(); ++i)
    {
      if ( ! electronhelper[i].selected ) continue;
      electronhelper_energy[n]	= electronhelper[i].energy;
      electronhelper_et[n]	= electronhelper[i].et;
      electronhelper_pt[n]	= electronhelper[i].pt;
      electronhelper_phi[n]	= electronhelper[i].phi;
      electronhelper_eta[n]	= electronhelper[i].eta;
      electronhelper_charge[n]	= electronhelper[i].charge;
      electronhelper_eSuperClusterOverP[n]	= electronhelper[i].eSuperClusterOverP;
      electronhelper_deltaEtaSuperClusterTrackAtVtx[n]	= electronhelper[i].deltaEtaSuperClusterTrackAtVtx;
      electronhelper_deltaPhiSuperClusterTrackAtVtx[n]	= electronhelper[i].deltaPhiSuperClusterTrackAtVtx;
      electronhelper_isEB[n]	= electronhelper[i].isEB;
      electronhelper_isEE[n]	= electronhelper[i].isEE;
      electronhelper_scSigmaIEtaIEta[n]	= electronhelper[i].scSigmaIEtaIEta;
      electronhelper_hadronicOverEm[n]	= electronhelper[i].hadronicOverEm;
      electronhelper_superCluster_energy[n]	= electronhelper[i].superCluster_energy;
      electronhelper_gsfTrack_trackerExpectedHitsInner_numberOfHits[n]	= electronhelper[i].gsfTrack_trackerExpectedHitsInner_numberOfHits;
      electronhelper_simpleEleId80relIso[n]	= electronhelper[i].simpleEleId80relIso;
      electronhelper_simpleEleId95relIso[n]	= electronhelper[i].simpleEleId95relIso;
      electronhelper_trackIso[n]	= electronhelper[i].trackIso;
      electronhelper_ecalIso[n]	= electronhelper[i].ecalIso;
      electronhelper_hcalIso[n]	= electronhelper[i].hcalIso;
      electronhelper_caloIso[n]	= electronhelper[i].caloIso;
      electronhelper_dxywrtPV[n]	= electronhelper[i].dxywrtPV;
      electronhelper_dzwrtPV[n]	= electronhelper[i].dzwrtPV;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < genjet.size(); ++i)
    {
      if ( ! genjet[i].selected ) continue;
      genjet_energy[n]	= genjet[i].energy;
      genjet_et[n]	= genjet[i].et;
      genjet_pt[n]	= genjet[i].pt;
      genjet_phi[n]	= genjet[i].phi;
      genjet_eta[n]	= genjet[i].eta;
      genjet_mass[n]	= genjet[i].mass;
      genjet_rapidity[n]	= genjet[i].rapidity;
      genjet_nConstituents[n]	= genjet[i].nConstituents;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < genparticlehelper.size(); ++i)
    {
      if ( ! genparticlehelper[i].selected ) continue;
      genparticlehelper_firstMother[n]	= genparticlehelper[i].firstMother;
      genparticlehelper_lastMother[n]	= genparticlehelper[i].lastMother;
      genparticlehelper_firstDaughter[n]	= genparticlehelper[i].firstDaughter;
      genparticlehelper_lastDaughter[n]	= genparticlehelper[i].lastDaughter;
      genparticlehelper_charge[n]	= genparticlehelper[i].charge;
      genparticlehelper_pdgId[n]	= genparticlehelper[i].pdgId;
      genparticlehelper_status[n]	= genparticlehelper[i].status;
      genparticlehelper_pt[n]	= genparticlehelper[i].pt;
      genparticlehelper_eta[n]	= genparticlehelper[i].eta;
      genparticlehelper_phi[n]	= genparticlehelper[i].phi;
      genparticlehelper_mass[n]	= genparticlehelper[i].mass;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < jethelper.size(); ++i)
    {
      if ( ! jethelper[i].selected ) continue;
      jethelper_energy[n]	= jethelper[i].energy;
      jethelper_uncor_energy[n]	= jethelper[i].uncor_energy;
      jethelper_et[n]	= jethelper[i].et;
      jethelper_uncor_et[n]	= jethelper[i].uncor_et;
      jethelper_pt[n]	= jethelper[i].pt;
      jethelper_uncor_pt[n]	= jethelper[i].uncor_pt;
      jethelper_phi[n]	= jethelper[i].phi;
      jethelper_eta[n]	= jethelper[i].eta;
      jethelper_rapidity[n]	= jethelper[i].rapidity;
      jethelper_mass[n]	= jethelper[i].mass;
      jethelper_jetArea[n]	= jethelper[i].jetArea;
      jethelper_jetCharge03[n]	= jethelper[i].jetCharge03;
      jethelper_jetCharge05[n]	= jethelper[i].jetCharge05;
      jethelper_jetCharge10[n]	= jethelper[i].jetCharge10;
      jethelper_chargedHadronEnergyFraction[n]	= jethelper[i].chargedHadronEnergyFraction;
      jethelper_neutralHadronEnergyFraction[n]	= jethelper[i].neutralHadronEnergyFraction;
      jethelper_chargedEmEnergyFraction[n]	= jethelper[i].chargedEmEnergyFraction;
      jethelper_neutralEmEnergyFraction[n]	= jethelper[i].neutralEmEnergyFraction;
      jethelper_photonEnergyFraction[n]	= jethelper[i].photonEnergyFraction;
      jethelper_muonEnergyFraction[n]	= jethelper[i].muonEnergyFraction;
      jethelper_chargedMultiplicity[n]	= jethelper[i].chargedMultiplicity;
      jethelper_nConstituents[n]	= jethelper[i].nConstituents;
      jethelper_partonFlavour[n]	= jethelper[i].partonFlavour;
      jethelper_trackCountingHighEffBJetTags[n]	= jethelper[i].trackCountingHighEffBJetTags;
      jethelper_trackCountingHighPurBJetTags[n]	= jethelper[i].trackCountingHighPurBJetTags;
      jethelper_jetProbabilityBJetTags[n]	= jethelper[i].jetProbabilityBJetTags;
      jethelper_jetBProbabilityBJetTags[n]	= jethelper[i].jetBProbabilityBJetTags;
      jethelper_combinedSecondaryVertexBJetTags[n]	= jethelper[i].combinedSecondaryVertexBJetTags;
      jethelper_combinedSecondaryVertexMVABJetTags[n]	= jethelper[i].combinedSecondaryVertexMVABJetTags;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < jethelper1.size(); ++i)
    {
      if ( ! jethelper1[i].selected ) continue;
      jethelper1_energy[n]	= jethelper1[i].energy;
      jethelper1_uncor_energy[n]	= jethelper1[i].uncor_energy;
      jethelper1_et[n]	= jethelper1[i].et;
      jethelper1_uncor_et[n]	= jethelper1[i].uncor_et;
      jethelper1_pt[n]	= jethelper1[i].pt;
      jethelper1_uncor_pt[n]	= jethelper1[i].uncor_pt;
      jethelper1_phi[n]	= jethelper1[i].phi;
      jethelper1_eta[n]	= jethelper1[i].eta;
      jethelper1_rapidity[n]	= jethelper1[i].rapidity;
      jethelper1_mass[n]	= jethelper1[i].mass;
      jethelper1_jetArea[n]	= jethelper1[i].jetArea;
      jethelper1_jetCharge03[n]	= jethelper1[i].jetCharge03;
      jethelper1_jetCharge05[n]	= jethelper1[i].jetCharge05;
      jethelper1_jetCharge10[n]	= jethelper1[i].jetCharge10;
      jethelper1_chargedHadronEnergyFraction[n]	= jethelper1[i].chargedHadronEnergyFraction;
      jethelper1_neutralHadronEnergyFraction[n]	= jethelper1[i].neutralHadronEnergyFraction;
      jethelper1_chargedEmEnergyFraction[n]	= jethelper1[i].chargedEmEnergyFraction;
      jethelper1_neutralEmEnergyFraction[n]	= jethelper1[i].neutralEmEnergyFraction;
      jethelper1_photonEnergyFraction[n]	= jethelper1[i].photonEnergyFraction;
      jethelper1_muonEnergyFraction[n]	= jethelper1[i].muonEnergyFraction;
      jethelper1_chargedMultiplicity[n]	= jethelper1[i].chargedMultiplicity;
      jethelper1_nConstituents[n]	= jethelper1[i].nConstituents;
      jethelper1_partonFlavour[n]	= jethelper1[i].partonFlavour;
      jethelper1_trackCountingHighEffBJetTags[n]	= jethelper1[i].trackCountingHighEffBJetTags;
      jethelper1_trackCountingHighPurBJetTags[n]	= jethelper1[i].trackCountingHighPurBJetTags;
      jethelper1_jetProbabilityBJetTags[n]	= jethelper1[i].jetProbabilityBJetTags;
      jethelper1_jetBProbabilityBJetTags[n]	= jethelper1[i].jetBProbabilityBJetTags;
      jethelper1_combinedSecondaryVertexBJetTags[n]	= jethelper1[i].combinedSecondaryVertexBJetTags;
      jethelper1_combinedSecondaryVertexMVABJetTags[n]	= jethelper1[i].combinedSecondaryVertexMVABJetTags;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < jethelper2.size(); ++i)
    {
      if ( ! jethelper2[i].selected ) continue;
      jethelper2_energy[n]	= jethelper2[i].energy;
      jethelper2_uncor_energy[n]	= jethelper2[i].uncor_energy;
      jethelper2_et[n]	= jethelper2[i].et;
      jethelper2_uncor_et[n]	= jethelper2[i].uncor_et;
      jethelper2_pt[n]	= jethelper2[i].pt;
      jethelper2_uncor_pt[n]	= jethelper2[i].uncor_pt;
      jethelper2_phi[n]	= jethelper2[i].phi;
      jethelper2_eta[n]	= jethelper2[i].eta;
      jethelper2_rapidity[n]	= jethelper2[i].rapidity;
      jethelper2_mass[n]	= jethelper2[i].mass;
      jethelper2_jetArea[n]	= jethelper2[i].jetArea;
      jethelper2_jetCharge03[n]	= jethelper2[i].jetCharge03;
      jethelper2_jetCharge05[n]	= jethelper2[i].jetCharge05;
      jethelper2_jetCharge10[n]	= jethelper2[i].jetCharge10;
      jethelper2_nConstituents[n]	= jethelper2[i].nConstituents;
      jethelper2_partonFlavour[n]	= jethelper2[i].partonFlavour;
      jethelper2_trackCountingHighEffBJetTags[n]	= jethelper2[i].trackCountingHighEffBJetTags;
      jethelper2_trackCountingHighPurBJetTags[n]	= jethelper2[i].trackCountingHighPurBJetTags;
      jethelper2_jetProbabilityBJetTags[n]	= jethelper2[i].jetProbabilityBJetTags;
      jethelper2_jetBProbabilityBJetTags[n]	= jethelper2[i].jetBProbabilityBJetTags;
      jethelper2_combinedSecondaryVertexBJetTags[n]	= jethelper2[i].combinedSecondaryVertexBJetTags;
      jethelper2_combinedSecondaryVertexMVABJetTags[n]	= jethelper2[i].combinedSecondaryVertexMVABJetTags;
      jethelper2_numberOfDaughters[n]	= jethelper2[i].numberOfDaughters;
      jethelper2_daughter_0_energy[n]	= jethelper2[i].daughter_0_energy;
      jethelper2_daughter_0_pt[n]	= jethelper2[i].daughter_0_pt;
      jethelper2_daughter_0_eta[n]	= jethelper2[i].daughter_0_eta;
      jethelper2_daughter_0_rapidity[n]	= jethelper2[i].daughter_0_rapidity;
      jethelper2_daughter_0_phi[n]	= jethelper2[i].daughter_0_phi;
      jethelper2_daughter_0_mass[n]	= jethelper2[i].daughter_0_mass;
      jethelper2_daughter_1_energy[n]	= jethelper2[i].daughter_1_energy;
      jethelper2_daughter_1_pt[n]	= jethelper2[i].daughter_1_pt;
      jethelper2_daughter_1_eta[n]	= jethelper2[i].daughter_1_eta;
      jethelper2_daughter_1_rapidity[n]	= jethelper2[i].daughter_1_rapidity;
      jethelper2_daughter_1_phi[n]	= jethelper2[i].daughter_1_phi;
      jethelper2_daughter_1_mass[n]	= jethelper2[i].daughter_1_mass;
      jethelper2_daughter_0_jetCharge03[n]	= jethelper2[i].daughter_0_jetCharge03;
      jethelper2_daughter_0_jetCharge05[n]	= jethelper2[i].daughter_0_jetCharge05;
      jethelper2_daughter_0_jetCharge10[n]	= jethelper2[i].daughter_0_jetCharge10;
      jethelper2_daughter_1_jetCharge03[n]	= jethelper2[i].daughter_1_jetCharge03;
      jethelper2_daughter_1_jetCharge05[n]	= jethelper2[i].daughter_1_jetCharge05;
      jethelper2_daughter_1_jetCharge10[n]	= jethelper2[i].daughter_1_jetCharge10;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < jethelper3.size(); ++i)
    {
      if ( ! jethelper3[i].selected ) continue;
      jethelper3_energy[n]	= jethelper3[i].energy;
      jethelper3_uncor_energy[n]	= jethelper3[i].uncor_energy;
      jethelper3_et[n]	= jethelper3[i].et;
      jethelper3_uncor_et[n]	= jethelper3[i].uncor_et;
      jethelper3_pt[n]	= jethelper3[i].pt;
      jethelper3_uncor_pt[n]	= jethelper3[i].uncor_pt;
      jethelper3_phi[n]	= jethelper3[i].phi;
      jethelper3_eta[n]	= jethelper3[i].eta;
      jethelper3_rapidity[n]	= jethelper3[i].rapidity;
      jethelper3_mass[n]	= jethelper3[i].mass;
      jethelper3_jetArea[n]	= jethelper3[i].jetArea;
      jethelper3_jetCharge03[n]	= jethelper3[i].jetCharge03;
      jethelper3_jetCharge05[n]	= jethelper3[i].jetCharge05;
      jethelper3_jetCharge10[n]	= jethelper3[i].jetCharge10;
      jethelper3_chargedHadronEnergyFraction[n]	= jethelper3[i].chargedHadronEnergyFraction;
      jethelper3_neutralHadronEnergyFraction[n]	= jethelper3[i].neutralHadronEnergyFraction;
      jethelper3_chargedEmEnergyFraction[n]	= jethelper3[i].chargedEmEnergyFraction;
      jethelper3_neutralEmEnergyFraction[n]	= jethelper3[i].neutralEmEnergyFraction;
      jethelper3_photonEnergyFraction[n]	= jethelper3[i].photonEnergyFraction;
      jethelper3_muonEnergyFraction[n]	= jethelper3[i].muonEnergyFraction;
      jethelper3_chargedMultiplicity[n]	= jethelper3[i].chargedMultiplicity;
      jethelper3_nConstituents[n]	= jethelper3[i].nConstituents;
      jethelper3_partonFlavour[n]	= jethelper3[i].partonFlavour;
      jethelper3_trackCountingHighEffBJetTags[n]	= jethelper3[i].trackCountingHighEffBJetTags;
      jethelper3_trackCountingHighPurBJetTags[n]	= jethelper3[i].trackCountingHighPurBJetTags;
      jethelper3_jetProbabilityBJetTags[n]	= jethelper3[i].jetProbabilityBJetTags;
      jethelper3_jetBProbabilityBJetTags[n]	= jethelper3[i].jetBProbabilityBJetTags;
      jethelper3_combinedSecondaryVertexBJetTags[n]	= jethelper3[i].combinedSecondaryVertexBJetTags;
      jethelper3_combinedSecondaryVertexMVABJetTags[n]	= jethelper3[i].combinedSecondaryVertexMVABJetTags;
      jethelper3_qjetsvolatility[n]	= jethelper3[i].qjetsvolatility;
      jethelper3_tau1[n]	= jethelper3[i].tau1;
      jethelper3_tau2[n]	= jethelper3[i].tau2;
      jethelper3_tau3[n]	= jethelper3[i].tau3;
      jethelper3_C2beta17[n]	= jethelper3[i].C2beta17;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < jethelper4.size(); ++i)
    {
      if ( ! jethelper4[i].selected ) continue;
      jethelper4_energy[n]	= jethelper4[i].energy;
      jethelper4_uncor_energy[n]	= jethelper4[i].uncor_energy;
      jethelper4_et[n]	= jethelper4[i].et;
      jethelper4_uncor_et[n]	= jethelper4[i].uncor_et;
      jethelper4_pt[n]	= jethelper4[i].pt;
      jethelper4_uncor_pt[n]	= jethelper4[i].uncor_pt;
      jethelper4_phi[n]	= jethelper4[i].phi;
      jethelper4_eta[n]	= jethelper4[i].eta;
      jethelper4_rapidity[n]	= jethelper4[i].rapidity;
      jethelper4_mass[n]	= jethelper4[i].mass;
      jethelper4_jetArea[n]	= jethelper4[i].jetArea;
      jethelper4_jetCharge03[n]	= jethelper4[i].jetCharge03;
      jethelper4_jetCharge05[n]	= jethelper4[i].jetCharge05;
      jethelper4_jetCharge10[n]	= jethelper4[i].jetCharge10;
      jethelper4_nConstituents[n]	= jethelper4[i].nConstituents;
      jethelper4_partonFlavour[n]	= jethelper4[i].partonFlavour;
      jethelper4_trackCountingHighEffBJetTags[n]	= jethelper4[i].trackCountingHighEffBJetTags;
      jethelper4_trackCountingHighPurBJetTags[n]	= jethelper4[i].trackCountingHighPurBJetTags;
      jethelper4_jetProbabilityBJetTags[n]	= jethelper4[i].jetProbabilityBJetTags;
      jethelper4_jetBProbabilityBJetTags[n]	= jethelper4[i].jetBProbabilityBJetTags;
      jethelper4_combinedSecondaryVertexBJetTags[n]	= jethelper4[i].combinedSecondaryVertexBJetTags;
      jethelper4_combinedSecondaryVertexMVABJetTags[n]	= jethelper4[i].combinedSecondaryVertexMVABJetTags;
      jethelper4_numberOfDaughters[n]	= jethelper4[i].numberOfDaughters;
      jethelper4_daughter_0_energy[n]	= jethelper4[i].daughter_0_energy;
      jethelper4_daughter_0_pt[n]	= jethelper4[i].daughter_0_pt;
      jethelper4_daughter_0_eta[n]	= jethelper4[i].daughter_0_eta;
      jethelper4_daughter_0_rapidity[n]	= jethelper4[i].daughter_0_rapidity;
      jethelper4_daughter_0_phi[n]	= jethelper4[i].daughter_0_phi;
      jethelper4_daughter_0_mass[n]	= jethelper4[i].daughter_0_mass;
      jethelper4_daughter_1_energy[n]	= jethelper4[i].daughter_1_energy;
      jethelper4_daughter_1_pt[n]	= jethelper4[i].daughter_1_pt;
      jethelper4_daughter_1_eta[n]	= jethelper4[i].daughter_1_eta;
      jethelper4_daughter_1_rapidity[n]	= jethelper4[i].daughter_1_rapidity;
      jethelper4_daughter_1_phi[n]	= jethelper4[i].daughter_1_phi;
      jethelper4_daughter_1_mass[n]	= jethelper4[i].daughter_1_mass;
      jethelper4_tau1[n]	= jethelper4[i].tau1;
      jethelper4_tau2[n]	= jethelper4[i].tau2;
      jethelper4_tau3[n]	= jethelper4[i].tau3;
      jethelper4_daughter_0_jetCharge03[n]	= jethelper4[i].daughter_0_jetCharge03;
      jethelper4_daughter_0_jetCharge05[n]	= jethelper4[i].daughter_0_jetCharge05;
      jethelper4_daughter_0_jetCharge10[n]	= jethelper4[i].daughter_0_jetCharge10;
      jethelper4_daughter_1_jetCharge03[n]	= jethelper4[i].daughter_1_jetCharge03;
      jethelper4_daughter_1_jetCharge05[n]	= jethelper4[i].daughter_1_jetCharge05;
      jethelper4_daughter_1_jetCharge10[n]	= jethelper4[i].daughter_1_jetCharge10;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < jethelper5.size(); ++i)
    {
      if ( ! jethelper5[i].selected ) continue;
      jethelper5_energy[n]	= jethelper5[i].energy;
      jethelper5_uncor_energy[n]	= jethelper5[i].uncor_energy;
      jethelper5_et[n]	= jethelper5[i].et;
      jethelper5_uncor_et[n]	= jethelper5[i].uncor_et;
      jethelper5_pt[n]	= jethelper5[i].pt;
      jethelper5_uncor_pt[n]	= jethelper5[i].uncor_pt;
      jethelper5_phi[n]	= jethelper5[i].phi;
      jethelper5_eta[n]	= jethelper5[i].eta;
      jethelper5_rapidity[n]	= jethelper5[i].rapidity;
      jethelper5_mass[n]	= jethelper5[i].mass;
      jethelper5_jetArea[n]	= jethelper5[i].jetArea;
      jethelper5_jetCharge03[n]	= jethelper5[i].jetCharge03;
      jethelper5_jetCharge05[n]	= jethelper5[i].jetCharge05;
      jethelper5_jetCharge10[n]	= jethelper5[i].jetCharge10;
      jethelper5_chargedHadronEnergyFraction[n]	= jethelper5[i].chargedHadronEnergyFraction;
      jethelper5_neutralHadronEnergyFraction[n]	= jethelper5[i].neutralHadronEnergyFraction;
      jethelper5_chargedEmEnergyFraction[n]	= jethelper5[i].chargedEmEnergyFraction;
      jethelper5_neutralEmEnergyFraction[n]	= jethelper5[i].neutralEmEnergyFraction;
      jethelper5_photonEnergyFraction[n]	= jethelper5[i].photonEnergyFraction;
      jethelper5_muonEnergyFraction[n]	= jethelper5[i].muonEnergyFraction;
      jethelper5_chargedMultiplicity[n]	= jethelper5[i].chargedMultiplicity;
      jethelper5_nConstituents[n]	= jethelper5[i].nConstituents;
      jethelper5_partonFlavour[n]	= jethelper5[i].partonFlavour;
      jethelper5_trackCountingHighEffBJetTags[n]	= jethelper5[i].trackCountingHighEffBJetTags;
      jethelper5_trackCountingHighPurBJetTags[n]	= jethelper5[i].trackCountingHighPurBJetTags;
      jethelper5_jetProbabilityBJetTags[n]	= jethelper5[i].jetProbabilityBJetTags;
      jethelper5_jetBProbabilityBJetTags[n]	= jethelper5[i].jetBProbabilityBJetTags;
      jethelper5_combinedSecondaryVertexBJetTags[n]	= jethelper5[i].combinedSecondaryVertexBJetTags;
      jethelper5_combinedSecondaryVertexMVABJetTags[n]	= jethelper5[i].combinedSecondaryVertexMVABJetTags;
      jethelper5_qjetsvolatility[n]	= jethelper5[i].qjetsvolatility;
      jethelper5_tau1[n]	= jethelper5[i].tau1;
      jethelper5_tau2[n]	= jethelper5[i].tau2;
      jethelper5_tau3[n]	= jethelper5[i].tau3;
      jethelper5_C2beta17[n]	= jethelper5[i].C2beta17;
      jethelper5_getNcharged01[n]	= jethelper5[i].getNcharged01;
      jethelper5_getNneutral01[n]	= jethelper5[i].getNneutral01;
      jethelper5_getChargedPt0[n]	= jethelper5[i].getChargedPt0;
      jethelper5_getChargedPt1[n]	= jethelper5[i].getChargedPt1;
      jethelper5_getChargedPt2[n]	= jethelper5[i].getChargedPt2;
      jethelper5_getChargedPt3[n]	= jethelper5[i].getChargedPt3;
      jethelper5_getPt0[n]	= jethelper5[i].getPt0;
      jethelper5_getPt1[n]	= jethelper5[i].getPt1;
      jethelper5_getPt2[n]	= jethelper5[i].getPt2;
      jethelper5_getPt3[n]	= jethelper5[i].getPt3;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < jethelper6.size(); ++i)
    {
      if ( ! jethelper6[i].selected ) continue;
      jethelper6_energy[n]	= jethelper6[i].energy;
      jethelper6_et[n]	= jethelper6[i].et;
      jethelper6_pt[n]	= jethelper6[i].pt;
      jethelper6_phi[n]	= jethelper6[i].phi;
      jethelper6_eta[n]	= jethelper6[i].eta;
      jethelper6_rapidity[n]	= jethelper6[i].rapidity;
      jethelper6_mass[n]	= jethelper6[i].mass;
      jethelper6_jetArea[n]	= jethelper6[i].jetArea;
      jethelper6_genJetCharge03[n]	= jethelper6[i].genJetCharge03;
      jethelper6_genJetCharge05[n]	= jethelper6[i].genJetCharge05;
      jethelper6_genJetCharge10[n]	= jethelper6[i].genJetCharge10;
      jethelper6_nConstituents[n]	= jethelper6[i].nConstituents;
      jethelper6_partonFlavour[n]	= jethelper6[i].partonFlavour;
      jethelper6_genC2beta17[n]	= jethelper6[i].genC2beta17;
      jethelper6_genC2beta17CHS[n]	= jethelper6[i].genC2beta17CHS;
      jethelper6_genTau1[n]	= jethelper6[i].genTau1;
      jethelper6_genTau2[n]	= jethelper6[i].genTau2;
      jethelper6_genTau3[n]	= jethelper6[i].genTau3;
      jethelper6_genTau1Pt2[n]	= jethelper6[i].genTau1Pt2;
      jethelper6_genTau2Pt2[n]	= jethelper6[i].genTau2Pt2;
      jethelper6_genTau1Pt5[n]	= jethelper6[i].genTau1Pt5;
      jethelper6_genTau2Pt5[n]	= jethelper6[i].genTau2Pt5;
      jethelper6_genTau1CHS[n]	= jethelper6[i].genTau1CHS;
      jethelper6_genTau2CHS[n]	= jethelper6[i].genTau2CHS;
      jethelper6_genNCHS[n]	= jethelper6[i].genNCHS;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < jethelper7.size(); ++i)
    {
      if ( ! jethelper7[i].selected ) continue;
      jethelper7_energy[n]	= jethelper7[i].energy;
      jethelper7_et[n]	= jethelper7[i].et;
      jethelper7_pt[n]	= jethelper7[i].pt;
      jethelper7_phi[n]	= jethelper7[i].phi;
      jethelper7_eta[n]	= jethelper7[i].eta;
      jethelper7_rapidity[n]	= jethelper7[i].rapidity;
      jethelper7_mass[n]	= jethelper7[i].mass;
      jethelper7_jetArea[n]	= jethelper7[i].jetArea;
      jethelper7_jetCharge03[n]	= jethelper7[i].jetCharge03;
      jethelper7_jetCharge05[n]	= jethelper7[i].jetCharge05;
      jethelper7_jetCharge10[n]	= jethelper7[i].jetCharge10;
      jethelper7_nConstituents[n]	= jethelper7[i].nConstituents;
      jethelper7_partonFlavour[n]	= jethelper7[i].partonFlavour;
      jethelper7_numberOfDaughters[n]	= jethelper7[i].numberOfDaughters;
      jethelper7_daughter_0_energy[n]	= jethelper7[i].daughter_0_energy;
      jethelper7_daughter_0_pt[n]	= jethelper7[i].daughter_0_pt;
      jethelper7_daughter_0_eta[n]	= jethelper7[i].daughter_0_eta;
      jethelper7_daughter_0_rapidity[n]	= jethelper7[i].daughter_0_rapidity;
      jethelper7_daughter_0_phi[n]	= jethelper7[i].daughter_0_phi;
      jethelper7_daughter_0_mass[n]	= jethelper7[i].daughter_0_mass;
      jethelper7_daughter_1_energy[n]	= jethelper7[i].daughter_1_energy;
      jethelper7_daughter_1_pt[n]	= jethelper7[i].daughter_1_pt;
      jethelper7_daughter_1_eta[n]	= jethelper7[i].daughter_1_eta;
      jethelper7_daughter_1_rapidity[n]	= jethelper7[i].daughter_1_rapidity;
      jethelper7_daughter_1_phi[n]	= jethelper7[i].daughter_1_phi;
      jethelper7_daughter_1_mass[n]	= jethelper7[i].daughter_1_mass;
      jethelper7_genTau1[n]	= jethelper7[i].genTau1;
      jethelper7_genTau2[n]	= jethelper7[i].genTau2;
      jethelper7_genTau3[n]	= jethelper7[i].genTau3;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < leafcandidate.size(); ++i)
    {
      if ( ! leafcandidate[i].selected ) continue;
      leafcandidate_charge[n]	= leafcandidate[i].charge;
      leafcandidate_energy[n]	= leafcandidate[i].energy;
      leafcandidate_et[n]	= leafcandidate[i].et;
      leafcandidate_pt[n]	= leafcandidate[i].pt;
      leafcandidate_phi[n]	= leafcandidate[i].phi;
      leafcandidate_eta[n]	= leafcandidate[i].eta;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < met.size(); ++i)
    {
      if ( ! met[i].selected ) continue;
      met_energy[n]	= met[i].energy;
      met_et[n]	= met[i].et;
      met_pt[n]	= met[i].pt;
      met_phi[n]	= met[i].phi;
      met_sumEt[n]	= met[i].sumEt;
      met_mEtSig[n]	= met[i].mEtSig;
      met_significance[n]	= met[i].significance;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < met1.size(); ++i)
    {
      if ( ! met1[i].selected ) continue;
      met1_energy[n]	= met1[i].energy;
      met1_et[n]	= met1[i].et;
      met1_pt[n]	= met1[i].pt;
      met1_phi[n]	= met1[i].phi;
      met1_sumEt[n]	= met1[i].sumEt;
      met1_mEtSig[n]	= met1[i].mEtSig;
      met1_significance[n]	= met1[i].significance;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < muonhelper.size(); ++i)
    {
      if ( ! muonhelper[i].selected ) continue;
      muonhelper_energy[n]	= muonhelper[i].energy;
      muonhelper_et[n]	= muonhelper[i].et;
      muonhelper_pt[n]	= muonhelper[i].pt;
      muonhelper_phi[n]	= muonhelper[i].phi;
      muonhelper_eta[n]	= muonhelper[i].eta;
      muonhelper_charge[n]	= muonhelper[i].charge;
      muonhelper_TMOneStationTight[n]	= muonhelper[i].TMOneStationTight;
      muonhelper_isGlobalMuon[n]	= muonhelper[i].isGlobalMuon;
      muonhelper_isTrackerMuon[n]	= muonhelper[i].isTrackerMuon;
      muonhelper_isPFMuon[n]	= muonhelper[i].isPFMuon;
      muonhelper_track_hitPattern_trackerLayersWithMeasurement[n]	= muonhelper[i].track_hitPattern_trackerLayersWithMeasurement;
      muonhelper_innerTrack_hitPattern_pixelLayersWithMeasurement[n]	= muonhelper[i].innerTrack_hitPattern_pixelLayersWithMeasurement;
      muonhelper_innerTrack_normalizedChi2[n]	= muonhelper[i].innerTrack_normalizedChi2;
      muonhelper_globalTrack_normalizedChi2[n]	= muonhelper[i].globalTrack_normalizedChi2;
      muonhelper_globalTrack_hitPattern_numberOfValidMuonHits[n]	= muonhelper[i].globalTrack_hitPattern_numberOfValidMuonHits;
      muonhelper_numberOfMatchedStations[n]	= muonhelper[i].numberOfMatchedStations;
      muonhelper_dB[n]	= muonhelper[i].dB;
      muonhelper_innerTrack_hitPattern_numberOfValidPixelHits[n]	= muonhelper[i].innerTrack_hitPattern_numberOfValidPixelHits;
      muonhelper_pfIsolationR04_sumChargedHadronPt[n]	= muonhelper[i].pfIsolationR04_sumChargedHadronPt;
      muonhelper_pfIsolationR04_sumChargedParticlePt[n]	= muonhelper[i].pfIsolationR04_sumChargedParticlePt;
      muonhelper_pfIsolationR04_sumNeutralHadronEt[n]	= muonhelper[i].pfIsolationR04_sumNeutralHadronEt;
      muonhelper_pfIsolationR04_sumPhotonEt[n]	= muonhelper[i].pfIsolationR04_sumPhotonEt;
      muonhelper_pfIsolationR04_sumNeutralHadronEtHighThreshold[n]	= muonhelper[i].pfIsolationR04_sumNeutralHadronEtHighThreshold;
      muonhelper_pfIsolationR04_sumPhotonEtHighThreshold[n]	= muonhelper[i].pfIsolationR04_sumPhotonEtHighThreshold;
      muonhelper_pfIsolationR04_sumPUPt[n]	= muonhelper[i].pfIsolationR04_sumPUPt;
      muonhelper_dxywrtPV[n]	= muonhelper[i].dxywrtPV;
      muonhelper_dzwrtPV[n]	= muonhelper[i].dzwrtPV;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < pileupsummaryinfo.size(); ++i)
    {
      if ( ! pileupsummaryinfo[i].selected ) continue;
      pileupsummaryinfo_getPU_NumInteractions[n]	= pileupsummaryinfo[i].getPU_NumInteractions;
      pileupsummaryinfo_getBunchCrossing[n]	= pileupsummaryinfo[i].getBunchCrossing;
      pileupsummaryinfo_getTrueNumInteractions[n]	= pileupsummaryinfo[i].getTrueNumInteractions;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < tau.size(); ++i)
    {
      if ( ! tau[i].selected ) continue;
      tau_energy[n]	= tau[i].energy;
      tau_et[n]	= tau[i].et;
      tau_pt[n]	= tau[i].pt;
      tau_phi[n]	= tau[i].phi;
      tau_eta[n]	= tau[i].eta;
      tau_byLooseCombinedIsolationDeltaBetaCorr[n]	= tau[i].byLooseCombinedIsolationDeltaBetaCorr;
      tau_byMediumCombinedIsolationDeltaBetaCorr[n]	= tau[i].byMediumCombinedIsolationDeltaBetaCorr;
      tau_trackIso[n]	= tau[i].trackIso;
      tau_ecalIso[n]	= tau[i].ecalIso;
      tau_hcalIso[n]	= tau[i].hcalIso;
      tau_caloIso[n]	= tau[i].caloIso;
      n++;
    }

  n = 0;
  for(unsigned int i=0; i < vertex.size(); ++i)
    {
      if ( ! vertex[i].selected ) continue;
      vertex_isFake[n]	= vertex[i].isFake;
      vertex_chi2[n]	= vertex[i].chi2;
      vertex_ndof[n]	= vertex[i].ndof;
      vertex_x[n]	= vertex[i].x;
      vertex_y[n]	= vertex[i].y;
      vertex_z[n]	= vertex[i].z;
      vertex_xError[n]	= vertex[i].xError;
      vertex_yError[n]	= vertex[i].yError;
      vertex_zError[n]	= vertex[i].zError;
      vertex_position_Rho[n]	= vertex[i].position_Rho;
      n++;
    }
  fillObjectsCalled = false;
}

//-----------------------------------------------------------------------------
// -- Select variables to be read
//-----------------------------------------------------------------------------
void selectVariables(itreestream& stream)
{
  stream.select("recoCaloMET_met.energy", calomet1_energy);
  stream.select("recoCaloMET_met.mEtSig", calomet1_mEtSig);
  stream.select("recoCaloMET_met.phi", calomet1_phi);
  stream.select("recoCaloMET_met.pt", calomet1_pt);
  stream.select("recoCaloMET_met.significance", calomet1_significance);
  stream.select("recoCaloMET_met.sumEt", calomet1_sumEt);
  stream.select("recoCaloMET_corMetGlobalMuons.energy", calomet_energy);
  stream.select("recoCaloMET_corMetGlobalMuons.mEtSig", calomet_mEtSig);
  stream.select("recoCaloMET_corMetGlobalMuons.phi", calomet_phi);
  stream.select("recoCaloMET_corMetGlobalMuons.pt", calomet_pt);
  stream.select("recoCaloMET_corMetGlobalMuons.significance", calomet_significance);
  stream.select("recoCaloMET_corMetGlobalMuons.sumEt", calomet_sumEt);
  stream.select("cmgBaseMET_razorMJMetUp.energy", cmgbasemet1_energy);
  stream.select("cmgBaseMET_razorMJMetUp.et", cmgbasemet1_et);
  stream.select("cmgBaseMET_razorMJMetUp.phi", cmgbasemet1_phi);
  stream.select("cmgBaseMET_razorMJMetUp.pt", cmgbasemet1_pt);
  stream.select("cmgBaseMET_razorMJMetUp.sumEt", cmgbasemet1_sumEt);
  stream.select("cmgBaseMET_cmgPFMET.energy", cmgbasemet2_energy);
  stream.select("cmgBaseMET_cmgPFMET.et", cmgbasemet2_et);
  stream.select("cmgBaseMET_cmgPFMET.phi", cmgbasemet2_phi);
  stream.select("cmgBaseMET_cmgPFMET.pt", cmgbasemet2_pt);
  stream.select("cmgBaseMET_cmgPFMET.sumEt", cmgbasemet2_sumEt);
  stream.select("cmgBaseMET_razorMJMetDown.energy", cmgbasemet_energy);
  stream.select("cmgBaseMET_razorMJMetDown.et", cmgbasemet_et);
  stream.select("cmgBaseMET_razorMJMetDown.phi", cmgbasemet_phi);
  stream.select("cmgBaseMET_razorMJMetDown.pt", cmgbasemet_pt);
  stream.select("cmgBaseMET_razorMJMetDown.sumEt", cmgbasemet_sumEt);
  stream.select("cmgElectron_razorMJElectronTight.charge", cmgelectron1_charge);
  stream.select("cmgElectron_razorMJElectronTight.energy", cmgelectron1_energy);
  stream.select("cmgElectron_razorMJElectronTight.et", cmgelectron1_et);
  stream.select("cmgElectron_razorMJElectronTight.eta", cmgelectron1_eta);
  stream.select("cmgElectron_razorMJElectronTight.phi", cmgelectron1_phi);
  stream.select("cmgElectron_razorMJElectronTight.pt", cmgelectron1_pt);
  stream.select("cmgElectron_razorMJElectronLoose.charge", cmgelectron_charge);
  stream.select("cmgElectron_razorMJElectronLoose.energy", cmgelectron_energy);
  stream.select("cmgElectron_razorMJElectronLoose.et", cmgelectron_et);
  stream.select("cmgElectron_razorMJElectronLoose.eta", cmgelectron_eta);
  stream.select("cmgElectron_razorMJElectronLoose.phi", cmgelectron_phi);
  stream.select("cmgElectron_razorMJElectronLoose.pt", cmgelectron_pt);
  stream.select("cmgMuon_razorMJMuonTight.charge", cmgmuon1_charge);
  stream.select("cmgMuon_razorMJMuonTight.energy", cmgmuon1_energy);
  stream.select("cmgMuon_razorMJMuonTight.et", cmgmuon1_et);
  stream.select("cmgMuon_razorMJMuonTight.eta", cmgmuon1_eta);
  stream.select("cmgMuon_razorMJMuonTight.phi", cmgmuon1_phi);
  stream.select("cmgMuon_razorMJMuonTight.pt", cmgmuon1_pt);
  stream.select("cmgMuon_razorMJMuonLoose.charge", cmgmuon_charge);
  stream.select("cmgMuon_razorMJMuonLoose.energy", cmgmuon_energy);
  stream.select("cmgMuon_razorMJMuonLoose.et", cmgmuon_et);
  stream.select("cmgMuon_razorMJMuonLoose.eta", cmgmuon_eta);
  stream.select("cmgMuon_razorMJMuonLoose.phi", cmgmuon_phi);
  stream.select("cmgMuon_razorMJMuonLoose.pt", cmgmuon_pt);
  stream.select("cmgPFJet_cmgPFJetSelCHS.combinedSecondaryVertexBJetTags", cmgpfjet_combinedSecondaryVertexBJetTags);
  stream.select("cmgPFJet_cmgPFJetSelCHS.combinedSecondaryVertexMVABJetTags", cmgpfjet_combinedSecondaryVertexMVABJetTags);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_0_fraction", cmgpfjet_component_0_fraction);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_0_number", cmgpfjet_component_0_number);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_1_fraction", cmgpfjet_component_1_fraction);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_1_number", cmgpfjet_component_1_number);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_2_fraction", cmgpfjet_component_2_fraction);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_2_number", cmgpfjet_component_2_number);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_3_fraction", cmgpfjet_component_3_fraction);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_3_number", cmgpfjet_component_3_number);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_4_fraction", cmgpfjet_component_4_fraction);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_4_number", cmgpfjet_component_4_number);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_5_fraction", cmgpfjet_component_5_fraction);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_5_number", cmgpfjet_component_5_number);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_6_fraction", cmgpfjet_component_6_fraction);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_6_number", cmgpfjet_component_6_number);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_7_fraction", cmgpfjet_component_7_fraction);
  stream.select("cmgPFJet_cmgPFJetSelCHS.component_7_number", cmgpfjet_component_7_number);
  stream.select("cmgPFJet_cmgPFJetSelCHS.energy", cmgpfjet_energy);
  stream.select("cmgPFJet_cmgPFJetSelCHS.et", cmgpfjet_et);
  stream.select("cmgPFJet_cmgPFJetSelCHS.eta", cmgpfjet_eta);
  stream.select("cmgPFJet_cmgPFJetSelCHS.jetArea", cmgpfjet_jetArea);
  stream.select("cmgPFJet_cmgPFJetSelCHS.jetBProbabilityBJetTags", cmgpfjet_jetBProbabilityBJetTags);
  stream.select("cmgPFJet_cmgPFJetSelCHS.jetProbabilityBJetTags", cmgpfjet_jetProbabilityBJetTags);
  stream.select("cmgPFJet_cmgPFJetSelCHS.mass", cmgpfjet_mass);
  stream.select("cmgPFJet_cmgPFJetSelCHS.nConstituents", cmgpfjet_nConstituents);
  stream.select("cmgPFJet_cmgPFJetSelCHS.partonFlavour", cmgpfjet_partonFlavour);
  stream.select("cmgPFJet_cmgPFJetSelCHS.phi", cmgpfjet_phi);
  stream.select("cmgPFJet_cmgPFJetSelCHS.pt", cmgpfjet_pt);
  stream.select("cmgPFJet_cmgPFJetSelCHS.rapidity", cmgpfjet_rapidity);
  stream.select("cmgPFJet_cmgPFJetSelCHS.trackCountingHighEffBJetTag", cmgpfjet_trackCountingHighEffBJetTag);
  stream.select("cmgPFJet_cmgPFJetSelCHS.trackCountingHighPurBJetTags", cmgpfjet_trackCountingHighPurBJetTags);
  stream.select("cmgTau_razorMJTauTight.byTightCombinedIsolationDeltaBetaCorr", cmgtau1_byTightCombinedIsolationDeltaBetaCorr);
  stream.select("cmgTau_razorMJTauTight.charge", cmgtau1_charge);
  stream.select("cmgTau_razorMJTauTight.energy", cmgtau1_energy);
  stream.select("cmgTau_razorMJTauTight.et", cmgtau1_et);
  stream.select("cmgTau_razorMJTauTight.eta", cmgtau1_eta);
  stream.select("cmgTau_razorMJTauTight.phi", cmgtau1_phi);
  stream.select("cmgTau_razorMJTauTight.pt", cmgtau1_pt);
  stream.select("cmgTau_razorMJTauLoose.byLooseCombinedIsolationDeltaBetaCorr", cmgtau_byLooseCombinedIsolationDeltaBetaCorr);
  stream.select("cmgTau_razorMJTauLoose.charge", cmgtau_charge);
  stream.select("cmgTau_razorMJTauLoose.energy", cmgtau_energy);
  stream.select("cmgTau_razorMJTauLoose.et", cmgtau_et);
  stream.select("cmgTau_razorMJTauLoose.eta", cmgtau_eta);
  stream.select("cmgTau_razorMJTauLoose.phi", cmgtau_phi);
  stream.select("cmgTau_razorMJTauLoose.pt", cmgtau_pt);
  stream.select("patElectronHelper_patElectronsWithTrigger.caloIso", electronhelper_caloIso);
  stream.select("patElectronHelper_patElectronsWithTrigger.charge", electronhelper_charge);
  stream.select("patElectronHelper_patElectronsWithTrigger.deltaEtaSuperClusterTrackAtVtx", electronhelper_deltaEtaSuperClusterTrackAtVtx);
  stream.select("patElectronHelper_patElectronsWithTrigger.deltaPhiSuperClusterTrackAtVtx", electronhelper_deltaPhiSuperClusterTrackAtVtx);
  stream.select("patElectronHelper_patElectronsWithTrigger.dxywrtPV", electronhelper_dxywrtPV);
  stream.select("patElectronHelper_patElectronsWithTrigger.dzwrtPV", electronhelper_dzwrtPV);
  stream.select("patElectronHelper_patElectronsWithTrigger.eSuperClusterOverP", electronhelper_eSuperClusterOverP);
  stream.select("patElectronHelper_patElectronsWithTrigger.ecalIso", electronhelper_ecalIso);
  stream.select("patElectronHelper_patElectronsWithTrigger.energy", electronhelper_energy);
  stream.select("patElectronHelper_patElectronsWithTrigger.et", electronhelper_et);
  stream.select("patElectronHelper_patElectronsWithTrigger.eta", electronhelper_eta);
  stream.select("patElectronHelper_patElectronsWithTrigger.gsfTrack_trackerExpectedHitsInner_numberOfHits", electronhelper_gsfTrack_trackerExpectedHitsInner_numberOfHits);
  stream.select("patElectronHelper_patElectronsWithTrigger.hadronicOverEm", electronhelper_hadronicOverEm);
  stream.select("patElectronHelper_patElectronsWithTrigger.hcalIso", electronhelper_hcalIso);
  stream.select("patElectronHelper_patElectronsWithTrigger.isEB", electronhelper_isEB);
  stream.select("patElectronHelper_patElectronsWithTrigger.isEE", electronhelper_isEE);
  stream.select("patElectronHelper_patElectronsWithTrigger.phi", electronhelper_phi);
  stream.select("patElectronHelper_patElectronsWithTrigger.pt", electronhelper_pt);
  stream.select("patElectronHelper_patElectronsWithTrigger.scSigmaIEtaIEta", electronhelper_scSigmaIEtaIEta);
  stream.select("patElectronHelper_patElectronsWithTrigger.simpleEleId80relIso", electronhelper_simpleEleId80relIso);
  stream.select("patElectronHelper_patElectronsWithTrigger.simpleEleId95relIso", electronhelper_simpleEleId95relIso);
  stream.select("patElectronHelper_patElectronsWithTrigger.superCluster_energy", electronhelper_superCluster_energy);
  stream.select("patElectronHelper_patElectronsWithTrigger.trackIso", electronhelper_trackIso);
  stream.select("edmEventHelper_info.bunchCrossing", eventhelper_bunchCrossing);
  stream.select("edmEventHelper_info.event", eventhelper_event);
  stream.select("edmEventHelper_info.isRealData", eventhelper_isRealData);
  stream.select("edmEventHelper_info.luminosityBlock", eventhelper_luminosityBlock);
  stream.select("edmEventHelper_info.orbitNumber", eventhelper_orbitNumber);
  stream.select("edmEventHelper_info.run", eventhelper_run);
  stream.select("edmEventHelperExtra_info.numberOfPrimaryVertices", eventhelperextra_numberOfPrimaryVertices);
  stream.select("edmEventHelperExtra_info.trackIso", eventhelperextra_trackIso);
  stream.select("edmEventHelperExtra_info.trackIsoLep", eventhelperextra_trackIsoLep);
  stream.select("GenEventInfoProduct_generator.alphaQCD", geneventinfoproduct_alphaQCD);
  stream.select("GenEventInfoProduct_generator.alphaQED", geneventinfoproduct_alphaQED);
  stream.select("GenEventInfoProduct_generator.hasBinningValues", geneventinfoproduct_hasBinningValues);
  stream.select("GenEventInfoProduct_generator.hasPDF", geneventinfoproduct_hasPDF);
  stream.select("GenEventInfoProduct_generator.qScale", geneventinfoproduct_qScale);
  stream.select("GenEventInfoProduct_generator.signalProcessID", geneventinfoproduct_signalProcessID);
  stream.select("GenEventInfoProduct_generator.weight", geneventinfoproduct_weight);
  stream.select("recoGenJet_ak5GenJets.energy", genjet_energy);
  stream.select("recoGenJet_ak5GenJets.et", genjet_et);
  stream.select("recoGenJet_ak5GenJets.eta", genjet_eta);
  stream.select("recoGenJet_ak5GenJets.mass", genjet_mass);
  stream.select("recoGenJet_ak5GenJets.nConstituents", genjet_nConstituents);
  stream.select("recoGenJet_ak5GenJets.phi", genjet_phi);
  stream.select("recoGenJet_ak5GenJets.pt", genjet_pt);
  stream.select("recoGenJet_ak5GenJets.rapidity", genjet_rapidity);
  stream.select("recoGenParticleHelper_genParticles.charge", genparticlehelper_charge);
  stream.select("recoGenParticleHelper_genParticles.eta", genparticlehelper_eta);
  stream.select("recoGenParticleHelper_genParticles.firstDaughter", genparticlehelper_firstDaughter);
  stream.select("recoGenParticleHelper_genParticles.firstMother", genparticlehelper_firstMother);
  stream.select("recoGenParticleHelper_genParticles.lastDaughter", genparticlehelper_lastDaughter);
  stream.select("recoGenParticleHelper_genParticles.lastMother", genparticlehelper_lastMother);
  stream.select("recoGenParticleHelper_genParticles.mass", genparticlehelper_mass);
  stream.select("recoGenParticleHelper_genParticles.pdgId", genparticlehelper_pdgId);
  stream.select("recoGenParticleHelper_genParticles.phi", genparticlehelper_phi);
  stream.select("recoGenParticleHelper_genParticles.pt", genparticlehelper_pt);
  stream.select("recoGenParticleHelper_genParticles.status", genparticlehelper_status);
  stream.select("GenEventInfoProductHelper_generator.id1", geneventinfoproducthelper_id1);
  stream.select("GenEventInfoProductHelper_generator.id2", geneventinfoproducthelper_id2);
  stream.select("GenEventInfoProductHelper_generator.q", geneventinfoproducthelper_q);
  stream.select("GenEventInfoProductHelper_generator.x1", geneventinfoproducthelper_x1);
  stream.select("GenEventInfoProductHelper_generator.x2", geneventinfoproducthelper_x2);
  stream.select("GenRunInfoProduct_generator.crossSection", genruninfoproduct_crossSection);
  stream.select("GenRunInfoProduct_generator.filterEfficiency", genruninfoproduct_filterEfficiency);
  stream.select("patJetHelper_patJetsWithVarCHS.chargedEmEnergyFraction", jethelper1_chargedEmEnergyFraction);
  stream.select("patJetHelper_patJetsWithVarCHS.chargedHadronEnergyFraction", jethelper1_chargedHadronEnergyFraction);
  stream.select("patJetHelper_patJetsWithVarCHS.chargedMultiplicity", jethelper1_chargedMultiplicity);
  stream.select("patJetHelper_patJetsWithVarCHS.combinedSecondaryVertexBJetTags", jethelper1_combinedSecondaryVertexBJetTags);
  stream.select("patJetHelper_patJetsWithVarCHS.combinedSecondaryVertexMVABJetTags", jethelper1_combinedSecondaryVertexMVABJetTags);
  stream.select("patJetHelper_patJetsWithVarCHS.energy", jethelper1_energy);
  stream.select("patJetHelper_patJetsWithVarCHS.et", jethelper1_et);
  stream.select("patJetHelper_patJetsWithVarCHS.eta", jethelper1_eta);
  stream.select("patJetHelper_patJetsWithVarCHS.jetArea", jethelper1_jetArea);
  stream.select("patJetHelper_patJetsWithVarCHS.jetBProbabilityBJetTags", jethelper1_jetBProbabilityBJetTags);
  stream.select("patJetHelper_patJetsWithVarCHS.jetCharge03", jethelper1_jetCharge03);
  stream.select("patJetHelper_patJetsWithVarCHS.jetCharge05", jethelper1_jetCharge05);
  stream.select("patJetHelper_patJetsWithVarCHS.jetCharge10", jethelper1_jetCharge10);
  stream.select("patJetHelper_patJetsWithVarCHS.jetProbabilityBJetTags", jethelper1_jetProbabilityBJetTags);
  stream.select("patJetHelper_patJetsWithVarCHS.mass", jethelper1_mass);
  stream.select("patJetHelper_patJetsWithVarCHS.muonEnergyFraction", jethelper1_muonEnergyFraction);
  stream.select("patJetHelper_patJetsWithVarCHS.nConstituents", jethelper1_nConstituents);
  stream.select("patJetHelper_patJetsWithVarCHS.neutralEmEnergyFraction", jethelper1_neutralEmEnergyFraction);
  stream.select("patJetHelper_patJetsWithVarCHS.neutralHadronEnergyFraction", jethelper1_neutralHadronEnergyFraction);
  stream.select("patJetHelper_patJetsWithVarCHS.partonFlavour", jethelper1_partonFlavour);
  stream.select("patJetHelper_patJetsWithVarCHS.phi", jethelper1_phi);
  stream.select("patJetHelper_patJetsWithVarCHS.photonEnergyFraction", jethelper1_photonEnergyFraction);
  stream.select("patJetHelper_patJetsWithVarCHS.pt", jethelper1_pt);
  stream.select("patJetHelper_patJetsWithVarCHS.rapidity", jethelper1_rapidity);
  stream.select("patJetHelper_patJetsWithVarCHS.trackCountingHighEffBJetTags", jethelper1_trackCountingHighEffBJetTags);
  stream.select("patJetHelper_patJetsWithVarCHS.trackCountingHighPurBJetTags", jethelper1_trackCountingHighPurBJetTags);
  stream.select("patJetHelper_patJetsWithVarCHS.uncor_energy", jethelper1_uncor_energy);
  stream.select("patJetHelper_patJetsWithVarCHS.uncor_et", jethelper1_uncor_et);
  stream.select("patJetHelper_patJetsWithVarCHS.uncor_pt", jethelper1_uncor_pt);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.combinedSecondaryVertexBJetTags", jethelper2_combinedSecondaryVertexBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.combinedSecondaryVertexMVABJetTags", jethelper2_combinedSecondaryVertexMVABJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_0_energy", jethelper2_daughter_0_energy);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_0_eta", jethelper2_daughter_0_eta);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_0_jetCharge03", jethelper2_daughter_0_jetCharge03);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_0_jetCharge05", jethelper2_daughter_0_jetCharge05);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_0_jetCharge10", jethelper2_daughter_0_jetCharge10);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_0_mass", jethelper2_daughter_0_mass);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_0_phi", jethelper2_daughter_0_phi);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_0_pt", jethelper2_daughter_0_pt);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_0_rapidity", jethelper2_daughter_0_rapidity);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_1_energy", jethelper2_daughter_1_energy);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_1_eta", jethelper2_daughter_1_eta);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_1_jetCharge03", jethelper2_daughter_1_jetCharge03);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_1_jetCharge05", jethelper2_daughter_1_jetCharge05);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_1_jetCharge10", jethelper2_daughter_1_jetCharge10);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_1_mass", jethelper2_daughter_1_mass);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_1_phi", jethelper2_daughter_1_phi);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_1_pt", jethelper2_daughter_1_pt);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.daughter_1_rapidity", jethelper2_daughter_1_rapidity);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.energy", jethelper2_energy);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.et", jethelper2_et);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.eta", jethelper2_eta);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.jetArea", jethelper2_jetArea);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.jetBProbabilityBJetTags", jethelper2_jetBProbabilityBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.jetCharge03", jethelper2_jetCharge03);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.jetCharge05", jethelper2_jetCharge05);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.jetCharge10", jethelper2_jetCharge10);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.jetProbabilityBJetTags", jethelper2_jetProbabilityBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.mass", jethelper2_mass);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.nConstituents", jethelper2_nConstituents);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.numberOfDaughters", jethelper2_numberOfDaughters);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.partonFlavour", jethelper2_partonFlavour);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.phi", jethelper2_phi);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.pt", jethelper2_pt);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.rapidity", jethelper2_rapidity);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.trackCountingHighEffBJetTags", jethelper2_trackCountingHighEffBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.trackCountingHighPurBJetTags", jethelper2_trackCountingHighPurBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.uncor_energy", jethelper2_uncor_energy);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.uncor_et", jethelper2_uncor_et);
  stream.select("patJetHelper_selectedPatJetsAK7CHSpruned.uncor_pt", jethelper2_uncor_pt);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.C2beta17", jethelper3_C2beta17);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.chargedEmEnergyFraction", jethelper3_chargedEmEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.chargedHadronEnergyFraction", jethelper3_chargedHadronEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.chargedMultiplicity", jethelper3_chargedMultiplicity);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.combinedSecondaryVertexBJetTags", jethelper3_combinedSecondaryVertexBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.combinedSecondaryVertexMVABJetTags", jethelper3_combinedSecondaryVertexMVABJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.energy", jethelper3_energy);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.et", jethelper3_et);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.eta", jethelper3_eta);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.jetArea", jethelper3_jetArea);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.jetBProbabilityBJetTags", jethelper3_jetBProbabilityBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.jetCharge03", jethelper3_jetCharge03);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.jetCharge05", jethelper3_jetCharge05);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.jetCharge10", jethelper3_jetCharge10);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.jetProbabilityBJetTags", jethelper3_jetProbabilityBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.mass", jethelper3_mass);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.muonEnergyFraction", jethelper3_muonEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.nConstituents", jethelper3_nConstituents);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.neutralEmEnergyFraction", jethelper3_neutralEmEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.neutralHadronEnergyFraction", jethelper3_neutralHadronEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.partonFlavour", jethelper3_partonFlavour);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.phi", jethelper3_phi);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.photonEnergyFraction", jethelper3_photonEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.pt", jethelper3_pt);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.qjetsvolatility", jethelper3_qjetsvolatility);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.rapidity", jethelper3_rapidity);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.tau1", jethelper3_tau1);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.tau2", jethelper3_tau2);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.tau3", jethelper3_tau3);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.trackCountingHighEffBJetTags", jethelper3_trackCountingHighEffBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.trackCountingHighPurBJetTags", jethelper3_trackCountingHighPurBJetTags);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.uncor_energy", jethelper3_uncor_energy);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.uncor_et", jethelper3_uncor_et);
  stream.select("patJetHelper_selectedPatJetsAK7CHSwithQjets.uncor_pt", jethelper3_uncor_pt);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.combinedSecondaryVertexBJetTags", jethelper4_combinedSecondaryVertexBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.combinedSecondaryVertexMVABJetTags", jethelper4_combinedSecondaryVertexMVABJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_0_energy", jethelper4_daughter_0_energy);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_0_eta", jethelper4_daughter_0_eta);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_0_jetCharge03", jethelper4_daughter_0_jetCharge03);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_0_jetCharge05", jethelper4_daughter_0_jetCharge05);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_0_jetCharge10", jethelper4_daughter_0_jetCharge10);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_0_mass", jethelper4_daughter_0_mass);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_0_phi", jethelper4_daughter_0_phi);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_0_pt", jethelper4_daughter_0_pt);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_0_rapidity", jethelper4_daughter_0_rapidity);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_1_energy", jethelper4_daughter_1_energy);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_1_eta", jethelper4_daughter_1_eta);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_1_jetCharge03", jethelper4_daughter_1_jetCharge03);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_1_jetCharge05", jethelper4_daughter_1_jetCharge05);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_1_jetCharge10", jethelper4_daughter_1_jetCharge10);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_1_mass", jethelper4_daughter_1_mass);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_1_phi", jethelper4_daughter_1_phi);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_1_pt", jethelper4_daughter_1_pt);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.daughter_1_rapidity", jethelper4_daughter_1_rapidity);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.energy", jethelper4_energy);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.et", jethelper4_et);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.eta", jethelper4_eta);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.jetArea", jethelper4_jetArea);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.jetBProbabilityBJetTags", jethelper4_jetBProbabilityBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.jetCharge03", jethelper4_jetCharge03);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.jetCharge05", jethelper4_jetCharge05);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.jetCharge10", jethelper4_jetCharge10);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.jetProbabilityBJetTags", jethelper4_jetProbabilityBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.mass", jethelper4_mass);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.nConstituents", jethelper4_nConstituents);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.numberOfDaughters", jethelper4_numberOfDaughters);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.partonFlavour", jethelper4_partonFlavour);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.phi", jethelper4_phi);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.pt", jethelper4_pt);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.rapidity", jethelper4_rapidity);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.tau1", jethelper4_tau1);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.tau2", jethelper4_tau2);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.tau3", jethelper4_tau3);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.trackCountingHighEffBJetTags", jethelper4_trackCountingHighEffBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.trackCountingHighPurBJetTags", jethelper4_trackCountingHighPurBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.uncor_energy", jethelper4_uncor_energy);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.uncor_et", jethelper4_uncor_et);
  stream.select("patJetHelper_selectedPatJetsCA8CHSpruned.uncor_pt", jethelper4_uncor_pt);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.C2beta17", jethelper5_C2beta17);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.chargedEmEnergyFraction", jethelper5_chargedEmEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.chargedHadronEnergyFraction", jethelper5_chargedHadronEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.chargedMultiplicity", jethelper5_chargedMultiplicity);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.combinedSecondaryVertexBJetTags", jethelper5_combinedSecondaryVertexBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.combinedSecondaryVertexMVABJetTags", jethelper5_combinedSecondaryVertexMVABJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.energy", jethelper5_energy);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.et", jethelper5_et);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.eta", jethelper5_eta);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getChargedPt0", jethelper5_getChargedPt0);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getChargedPt1", jethelper5_getChargedPt1);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getChargedPt2", jethelper5_getChargedPt2);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getChargedPt3", jethelper5_getChargedPt3);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getNcharged01", jethelper5_getNcharged01);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getNneutral01", jethelper5_getNneutral01);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getPt0", jethelper5_getPt0);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getPt1", jethelper5_getPt1);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getPt2", jethelper5_getPt2);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.getPt3", jethelper5_getPt3);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.jetArea", jethelper5_jetArea);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.jetBProbabilityBJetTags", jethelper5_jetBProbabilityBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.jetCharge03", jethelper5_jetCharge03);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.jetCharge05", jethelper5_jetCharge05);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.jetCharge10", jethelper5_jetCharge10);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.jetProbabilityBJetTags", jethelper5_jetProbabilityBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.mass", jethelper5_mass);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.muonEnergyFraction", jethelper5_muonEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.nConstituents", jethelper5_nConstituents);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.neutralEmEnergyFraction", jethelper5_neutralEmEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.neutralHadronEnergyFraction", jethelper5_neutralHadronEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.partonFlavour", jethelper5_partonFlavour);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.phi", jethelper5_phi);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.photonEnergyFraction", jethelper5_photonEnergyFraction);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.pt", jethelper5_pt);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.qjetsvolatility", jethelper5_qjetsvolatility);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.rapidity", jethelper5_rapidity);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.tau1", jethelper5_tau1);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.tau2", jethelper5_tau2);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.tau3", jethelper5_tau3);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.trackCountingHighEffBJetTags", jethelper5_trackCountingHighEffBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.trackCountingHighPurBJetTags", jethelper5_trackCountingHighPurBJetTags);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.uncor_energy", jethelper5_uncor_energy);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.uncor_et", jethelper5_uncor_et);
  stream.select("patJetHelper_selectedPatJetsCA8CHSwithQjets.uncor_pt", jethelper5_uncor_pt);
  stream.select("patJetHelper_patGenJetsCA8CHS.energy", jethelper6_energy);
  stream.select("patJetHelper_patGenJetsCA8CHS.et", jethelper6_et);
  stream.select("patJetHelper_patGenJetsCA8CHS.eta", jethelper6_eta);
  stream.select("patJetHelper_patGenJetsCA8CHS.genC2beta17", jethelper6_genC2beta17);
  stream.select("patJetHelper_patGenJetsCA8CHS.genC2beta17CHS", jethelper6_genC2beta17CHS);
  stream.select("patJetHelper_patGenJetsCA8CHS.genJetCharge03", jethelper6_genJetCharge03);
  stream.select("patJetHelper_patGenJetsCA8CHS.genJetCharge05", jethelper6_genJetCharge05);
  stream.select("patJetHelper_patGenJetsCA8CHS.genJetCharge10", jethelper6_genJetCharge10);
  stream.select("patJetHelper_patGenJetsCA8CHS.genNCHS", jethelper6_genNCHS);
  stream.select("patJetHelper_patGenJetsCA8CHS.genTau1", jethelper6_genTau1);
  stream.select("patJetHelper_patGenJetsCA8CHS.genTau1CHS", jethelper6_genTau1CHS);
  stream.select("patJetHelper_patGenJetsCA8CHS.genTau1Pt2", jethelper6_genTau1Pt2);
  stream.select("patJetHelper_patGenJetsCA8CHS.genTau1Pt5", jethelper6_genTau1Pt5);
  stream.select("patJetHelper_patGenJetsCA8CHS.genTau2", jethelper6_genTau2);
  stream.select("patJetHelper_patGenJetsCA8CHS.genTau2CHS", jethelper6_genTau2CHS);
  stream.select("patJetHelper_patGenJetsCA8CHS.genTau2Pt2", jethelper6_genTau2Pt2);
  stream.select("patJetHelper_patGenJetsCA8CHS.genTau2Pt5", jethelper6_genTau2Pt5);
  stream.select("patJetHelper_patGenJetsCA8CHS.genTau3", jethelper6_genTau3);
  stream.select("patJetHelper_patGenJetsCA8CHS.jetArea", jethelper6_jetArea);
  stream.select("patJetHelper_patGenJetsCA8CHS.mass", jethelper6_mass);
  stream.select("patJetHelper_patGenJetsCA8CHS.nConstituents", jethelper6_nConstituents);
  stream.select("patJetHelper_patGenJetsCA8CHS.partonFlavour", jethelper6_partonFlavour);
  stream.select("patJetHelper_patGenJetsCA8CHS.phi", jethelper6_phi);
  stream.select("patJetHelper_patGenJetsCA8CHS.pt", jethelper6_pt);
  stream.select("patJetHelper_patGenJetsCA8CHS.rapidity", jethelper6_rapidity);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_0_energy", jethelper7_daughter_0_energy);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_0_eta", jethelper7_daughter_0_eta);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_0_mass", jethelper7_daughter_0_mass);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_0_phi", jethelper7_daughter_0_phi);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_0_pt", jethelper7_daughter_0_pt);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_0_rapidity", jethelper7_daughter_0_rapidity);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_1_energy", jethelper7_daughter_1_energy);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_1_eta", jethelper7_daughter_1_eta);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_1_mass", jethelper7_daughter_1_mass);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_1_phi", jethelper7_daughter_1_phi);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_1_pt", jethelper7_daughter_1_pt);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.daughter_1_rapidity", jethelper7_daughter_1_rapidity);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.energy", jethelper7_energy);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.et", jethelper7_et);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.eta", jethelper7_eta);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.genTau1", jethelper7_genTau1);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.genTau2", jethelper7_genTau2);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.genTau3", jethelper7_genTau3);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.jetArea", jethelper7_jetArea);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.jetCharge03", jethelper7_jetCharge03);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.jetCharge05", jethelper7_jetCharge05);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.jetCharge10", jethelper7_jetCharge10);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.mass", jethelper7_mass);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.nConstituents", jethelper7_nConstituents);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.numberOfDaughters", jethelper7_numberOfDaughters);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.partonFlavour", jethelper7_partonFlavour);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.phi", jethelper7_phi);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.pt", jethelper7_pt);
  stream.select("patJetHelper_patGenJetsCA8CHSpruned.rapidity", jethelper7_rapidity);
  stream.select("patJetHelper_patJetsWithVar.chargedEmEnergyFraction", jethelper_chargedEmEnergyFraction);
  stream.select("patJetHelper_patJetsWithVar.chargedHadronEnergyFraction", jethelper_chargedHadronEnergyFraction);
  stream.select("patJetHelper_patJetsWithVar.chargedMultiplicity", jethelper_chargedMultiplicity);
  stream.select("patJetHelper_patJetsWithVar.combinedSecondaryVertexBJetTags", jethelper_combinedSecondaryVertexBJetTags);
  stream.select("patJetHelper_patJetsWithVar.combinedSecondaryVertexMVABJetTags", jethelper_combinedSecondaryVertexMVABJetTags);
  stream.select("patJetHelper_patJetsWithVar.energy", jethelper_energy);
  stream.select("patJetHelper_patJetsWithVar.et", jethelper_et);
  stream.select("patJetHelper_patJetsWithVar.eta", jethelper_eta);
  stream.select("patJetHelper_patJetsWithVar.jetArea", jethelper_jetArea);
  stream.select("patJetHelper_patJetsWithVar.jetBProbabilityBJetTags", jethelper_jetBProbabilityBJetTags);
  stream.select("patJetHelper_patJetsWithVar.jetCharge03", jethelper_jetCharge03);
  stream.select("patJetHelper_patJetsWithVar.jetCharge05", jethelper_jetCharge05);
  stream.select("patJetHelper_patJetsWithVar.jetCharge10", jethelper_jetCharge10);
  stream.select("patJetHelper_patJetsWithVar.jetProbabilityBJetTags", jethelper_jetProbabilityBJetTags);
  stream.select("patJetHelper_patJetsWithVar.mass", jethelper_mass);
  stream.select("patJetHelper_patJetsWithVar.muonEnergyFraction", jethelper_muonEnergyFraction);
  stream.select("patJetHelper_patJetsWithVar.nConstituents", jethelper_nConstituents);
  stream.select("patJetHelper_patJetsWithVar.neutralEmEnergyFraction", jethelper_neutralEmEnergyFraction);
  stream.select("patJetHelper_patJetsWithVar.neutralHadronEnergyFraction", jethelper_neutralHadronEnergyFraction);
  stream.select("patJetHelper_patJetsWithVar.partonFlavour", jethelper_partonFlavour);
  stream.select("patJetHelper_patJetsWithVar.phi", jethelper_phi);
  stream.select("patJetHelper_patJetsWithVar.photonEnergyFraction", jethelper_photonEnergyFraction);
  stream.select("patJetHelper_patJetsWithVar.pt", jethelper_pt);
  stream.select("patJetHelper_patJetsWithVar.rapidity", jethelper_rapidity);
  stream.select("patJetHelper_patJetsWithVar.trackCountingHighEffBJetTags", jethelper_trackCountingHighEffBJetTags);
  stream.select("patJetHelper_patJetsWithVar.trackCountingHighPurBJetTags", jethelper_trackCountingHighPurBJetTags);
  stream.select("patJetHelper_patJetsWithVar.uncor_energy", jethelper_uncor_energy);
  stream.select("patJetHelper_patJetsWithVar.uncor_et", jethelper_uncor_et);
  stream.select("patJetHelper_patJetsWithVar.uncor_pt", jethelper_uncor_pt);
  stream.select("recoLeafCandidate_topGenInfo.charge", leafcandidate_charge);
  stream.select("recoLeafCandidate_topGenInfo.energy", leafcandidate_energy);
  stream.select("recoLeafCandidate_topGenInfo.et", leafcandidate_et);
  stream.select("recoLeafCandidate_topGenInfo.eta", leafcandidate_eta);
  stream.select("recoLeafCandidate_topGenInfo.phi", leafcandidate_phi);
  stream.select("recoLeafCandidate_topGenInfo.pt", leafcandidate_pt);
  stream.select("LHEEventProduct_source.hepeup_AQCDUP", lheeventproduct_hepeup_AQCDUP);
  stream.select("LHEEventProduct_source.hepeup_AQEDUP", lheeventproduct_hepeup_AQEDUP);
  stream.select("LHEEventProduct_source.hepeup_IDPRUP", lheeventproduct_hepeup_IDPRUP);
  stream.select("LHEEventProduct_source.hepeup_NUP", lheeventproduct_hepeup_NUP);
  stream.select("LHEEventProduct_source.hepeup_SCALUP", lheeventproduct_hepeup_SCALUP);
  stream.select("LHEEventProduct_source.hepeup_XWGTUP", lheeventproduct_hepeup_XWGTUP);
  stream.select("LHEEventProductHelper_source.mg", lheeventproducthelper_mg);
  stream.select("LHEEventProductHelper_source.mt1", lheeventproducthelper_mt1);
  stream.select("LHEEventProductHelper_source.mz1", lheeventproducthelper_mz1);
  stream.select("patMET_patMETsRaw.energy", met1_energy);
  stream.select("patMET_patMETsRaw.et", met1_et);
  stream.select("patMET_patMETsRaw.mEtSig", met1_mEtSig);
  stream.select("patMET_patMETsRaw.phi", met1_phi);
  stream.select("patMET_patMETsRaw.pt", met1_pt);
  stream.select("patMET_patMETsRaw.significance", met1_significance);
  stream.select("patMET_patMETsRaw.sumEt", met1_sumEt);
  stream.select("patMET_patMETs.energy", met_energy);
  stream.select("patMET_patMETs.et", met_et);
  stream.select("patMET_patMETs.mEtSig", met_mEtSig);
  stream.select("patMET_patMETs.phi", met_phi);
  stream.select("patMET_patMETs.pt", met_pt);
  stream.select("patMET_patMETs.significance", met_significance);
  stream.select("patMET_patMETs.sumEt", met_sumEt);
  stream.select("patMuonHelper_patMuonsWithTrigger.TMOneStationTight", muonhelper_TMOneStationTight);
  stream.select("patMuonHelper_patMuonsWithTrigger.charge", muonhelper_charge);
  stream.select("patMuonHelper_patMuonsWithTrigger.dB", muonhelper_dB);
  stream.select("patMuonHelper_patMuonsWithTrigger.dxywrtPV", muonhelper_dxywrtPV);
  stream.select("patMuonHelper_patMuonsWithTrigger.dzwrtPV", muonhelper_dzwrtPV);
  stream.select("patMuonHelper_patMuonsWithTrigger.energy", muonhelper_energy);
  stream.select("patMuonHelper_patMuonsWithTrigger.et", muonhelper_et);
  stream.select("patMuonHelper_patMuonsWithTrigger.eta", muonhelper_eta);
  stream.select("patMuonHelper_patMuonsWithTrigger.globalTrack_hitPattern_numberOfValidMuonHits", muonhelper_globalTrack_hitPattern_numberOfValidMuonHits);
  stream.select("patMuonHelper_patMuonsWithTrigger.globalTrack_normalizedChi2", muonhelper_globalTrack_normalizedChi2);
  stream.select("patMuonHelper_patMuonsWithTrigger.innerTrack_hitPattern_numberOfValidPixelHits", muonhelper_innerTrack_hitPattern_numberOfValidPixelHits);
  stream.select("patMuonHelper_patMuonsWithTrigger.innerTrack_hitPattern_pixelLayersWithMeasurement", muonhelper_innerTrack_hitPattern_pixelLayersWithMeasurement);
  stream.select("patMuonHelper_patMuonsWithTrigger.innerTrack_normalizedChi2", muonhelper_innerTrack_normalizedChi2);
  stream.select("patMuonHelper_patMuonsWithTrigger.isGlobalMuon", muonhelper_isGlobalMuon);
  stream.select("patMuonHelper_patMuonsWithTrigger.isPFMuon", muonhelper_isPFMuon);
  stream.select("patMuonHelper_patMuonsWithTrigger.isTrackerMuon", muonhelper_isTrackerMuon);
  stream.select("patMuonHelper_patMuonsWithTrigger.numberOfMatchedStations", muonhelper_numberOfMatchedStations);
  stream.select("patMuonHelper_patMuonsWithTrigger.pfIsolationR04_sumChargedHadronPt", muonhelper_pfIsolationR04_sumChargedHadronPt);
  stream.select("patMuonHelper_patMuonsWithTrigger.pfIsolationR04_sumChargedParticlePt", muonhelper_pfIsolationR04_sumChargedParticlePt);
  stream.select("patMuonHelper_patMuonsWithTrigger.pfIsolationR04_sumNeutralHadronEt", muonhelper_pfIsolationR04_sumNeutralHadronEt);
  stream.select("patMuonHelper_patMuonsWithTrigger.pfIsolationR04_sumNeutralHadronEtHighThreshold", muonhelper_pfIsolationR04_sumNeutralHadronEtHighThreshold);
  stream.select("patMuonHelper_patMuonsWithTrigger.pfIsolationR04_sumPUPt", muonhelper_pfIsolationR04_sumPUPt);
  stream.select("patMuonHelper_patMuonsWithTrigger.pfIsolationR04_sumPhotonEt", muonhelper_pfIsolationR04_sumPhotonEt);
  stream.select("patMuonHelper_patMuonsWithTrigger.pfIsolationR04_sumPhotonEtHighThreshold", muonhelper_pfIsolationR04_sumPhotonEtHighThreshold);
  stream.select("patMuonHelper_patMuonsWithTrigger.phi", muonhelper_phi);
  stream.select("patMuonHelper_patMuonsWithTrigger.pt", muonhelper_pt);
  stream.select("patMuonHelper_patMuonsWithTrigger.track_hitPattern_trackerLayersWithMeasurement", muonhelper_track_hitPattern_trackerLayersWithMeasurement);
  stream.select("nrecoCaloMET_corMetGlobalMuons", ncalomet);
  stream.select("nrecoCaloMET_met", ncalomet1);
  stream.select("ncmgBaseMET_razorMJMetDown", ncmgbasemet);
  stream.select("ncmgBaseMET_razorMJMetUp", ncmgbasemet1);
  stream.select("ncmgBaseMET_cmgPFMET", ncmgbasemet2);
  stream.select("ncmgElectron_razorMJElectronLoose", ncmgelectron);
  stream.select("ncmgElectron_razorMJElectronTight", ncmgelectron1);
  stream.select("ncmgMuon_razorMJMuonLoose", ncmgmuon);
  stream.select("ncmgMuon_razorMJMuonTight", ncmgmuon1);
  stream.select("ncmgPFJet_cmgPFJetSelCHS", ncmgpfjet);
  stream.select("ncmgTau_razorMJTauLoose", ncmgtau);
  stream.select("ncmgTau_razorMJTauTight", ncmgtau1);
  stream.select("npatElectronHelper_patElectronsWithTrigger", nelectronhelper);
  stream.select("nrecoGenJet_ak5GenJets", ngenjet);
  stream.select("nrecoGenParticleHelper_genParticles", ngenparticlehelper);
  stream.select("npatJetHelper_patJetsWithVar", njethelper);
  stream.select("npatJetHelper_patJetsWithVarCHS", njethelper1);
  stream.select("npatJetHelper_selectedPatJetsAK7CHSpruned", njethelper2);
  stream.select("npatJetHelper_selectedPatJetsAK7CHSwithQjets", njethelper3);
  stream.select("npatJetHelper_selectedPatJetsCA8CHSpruned", njethelper4);
  stream.select("npatJetHelper_selectedPatJetsCA8CHSwithQjets", njethelper5);
  stream.select("npatJetHelper_patGenJetsCA8CHS", njethelper6);
  stream.select("npatJetHelper_patGenJetsCA8CHSpruned", njethelper7);
  stream.select("nrecoLeafCandidate_topGenInfo", nleafcandidate);
  stream.select("npatMET_patMETs", nmet);
  stream.select("npatMET_patMETsRaw", nmet1);
  stream.select("npatMuonHelper_patMuonsWithTrigger", nmuonhelper);
  stream.select("nPileupSummaryInfo_addPileupInfo", npileupsummaryinfo);
  stream.select("npatTau_selectedPatTaus", ntau);
  stream.select("nrecoVertex_goodOfflinePrimaryVertices", nvertex);
  stream.select("PileupSummaryInfo_addPileupInfo.getBunchCrossing", pileupsummaryinfo_getBunchCrossing);
  stream.select("PileupSummaryInfo_addPileupInfo.getPU_NumInteractions", pileupsummaryinfo_getPU_NumInteractions);
  stream.select("PileupSummaryInfo_addPileupInfo.getTrueNumInteractions", pileupsummaryinfo_getTrueNumInteractions);
  stream.select("sdouble_kt6PFJets_rho.value", sdouble_kt6PFJets_rho_value);
  stream.select("sdouble_vertexWeightSummer12MC53X2012ABCDData.value", sdouble_vertexWeightSummer12MC53X2012ABCDData_value);
  stream.select("sint_hcallasereventfilter2012.value", sint_hcallasereventfilter2012_value);
  stream.select("sint_simpleGenInfo.value", sint_simpleGenInfo_value);
  stream.select("patTau_selectedPatTaus.byLooseCombinedIsolationDeltaBetaCorr", tau_byLooseCombinedIsolationDeltaBetaCorr);
  stream.select("patTau_selectedPatTaus.byMediumCombinedIsolationDeltaBetaCorr", tau_byMediumCombinedIsolationDeltaBetaCorr);
  stream.select("patTau_selectedPatTaus.caloIso", tau_caloIso);
  stream.select("patTau_selectedPatTaus.ecalIso", tau_ecalIso);
  stream.select("patTau_selectedPatTaus.energy", tau_energy);
  stream.select("patTau_selectedPatTaus.et", tau_et);
  stream.select("patTau_selectedPatTaus.eta", tau_eta);
  stream.select("patTau_selectedPatTaus.hcalIso", tau_hcalIso);
  stream.select("patTau_selectedPatTaus.phi", tau_phi);
  stream.select("patTau_selectedPatTaus.pt", tau_pt);
  stream.select("patTau_selectedPatTaus.trackIso", tau_trackIso);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.CSCTightHaloFilterPath", triggerresultshelper1_CSCTightHaloFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.EcalDeadCellTriggerPrimitiveFilterPath", triggerresultshelper1_EcalDeadCellTriggerPrimitiveFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.HBHENoiseFilterPath", triggerresultshelper1_HBHENoiseFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.eeBadScFilterPath", triggerresultshelper1_eeBadScFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.hcalLaserEventFilterPath", triggerresultshelper1_hcalLaserEventFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.metNoiseCleaningPath", triggerresultshelper1_metNoiseCleaningPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.noscrapingFilterPath", triggerresultshelper1_noscrapingFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.primaryVertexFilterPath", triggerresultshelper1_primaryVertexFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.trackingFailureFilterPath", triggerresultshelper1_trackingFailureFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.trkPOGFiltersPath", triggerresultshelper1_trkPOGFiltersPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.totalKinematicsFilterPath", triggerresultshelper2_totalKinematicsFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_PAT.trackIsolationMakerFilterPath", triggerresultshelper3_trackIsolationMakerFilterPath);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v10", triggerresultshelper_HLT_DiPFJetAve320_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v11", triggerresultshelper_HLT_DiPFJetAve320_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v12", triggerresultshelper_HLT_DiPFJetAve320_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v2", triggerresultshelper_HLT_DiPFJetAve320_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v3", triggerresultshelper_HLT_DiPFJetAve320_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v4", triggerresultshelper_HLT_DiPFJetAve320_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v5", triggerresultshelper_HLT_DiPFJetAve320_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v6", triggerresultshelper_HLT_DiPFJetAve320_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v7", triggerresultshelper_HLT_DiPFJetAve320_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v8", triggerresultshelper_HLT_DiPFJetAve320_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve320_v9", triggerresultshelper_HLT_DiPFJetAve320_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v10", triggerresultshelper_HLT_DiPFJetAve400_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v11", triggerresultshelper_HLT_DiPFJetAve400_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v12", triggerresultshelper_HLT_DiPFJetAve400_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v2", triggerresultshelper_HLT_DiPFJetAve400_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v3", triggerresultshelper_HLT_DiPFJetAve400_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v4", triggerresultshelper_HLT_DiPFJetAve400_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v5", triggerresultshelper_HLT_DiPFJetAve400_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v6", triggerresultshelper_HLT_DiPFJetAve400_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v7", triggerresultshelper_HLT_DiPFJetAve400_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v8", triggerresultshelper_HLT_DiPFJetAve400_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_DiPFJetAve400_v9", triggerresultshelper_HLT_DiPFJetAve400_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v11", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v12", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v6", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v7", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v8", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v9", triggerresultshelper_HLT_FatDiPFJetMass750_DR1p1_Deta1p5_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v1", triggerresultshelper_HLT_HT450_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v2", triggerresultshelper_HLT_HT450_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v3", triggerresultshelper_HLT_HT450_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v4", triggerresultshelper_HLT_HT450_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v5", triggerresultshelper_HLT_HT450_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v6", triggerresultshelper_HLT_HT450_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v7", triggerresultshelper_HLT_HT450_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v8", triggerresultshelper_HLT_HT450_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT450_v9", triggerresultshelper_HLT_HT450_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v1", triggerresultshelper_HLT_HT500_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v2", triggerresultshelper_HLT_HT500_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v3", triggerresultshelper_HLT_HT500_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v4", triggerresultshelper_HLT_HT500_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v5", triggerresultshelper_HLT_HT500_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v6", triggerresultshelper_HLT_HT500_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v7", triggerresultshelper_HLT_HT500_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v8", triggerresultshelper_HLT_HT500_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT500_v9", triggerresultshelper_HLT_HT500_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v1", triggerresultshelper_HLT_HT550_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v2", triggerresultshelper_HLT_HT550_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v3", triggerresultshelper_HLT_HT550_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v4", triggerresultshelper_HLT_HT550_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v5", triggerresultshelper_HLT_HT550_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v6", triggerresultshelper_HLT_HT550_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v7", triggerresultshelper_HLT_HT550_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v8", triggerresultshelper_HLT_HT550_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT550_v9", triggerresultshelper_HLT_HT550_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v1", triggerresultshelper_HLT_HT650_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v2", triggerresultshelper_HLT_HT650_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v3", triggerresultshelper_HLT_HT650_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v4", triggerresultshelper_HLT_HT650_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v5", triggerresultshelper_HLT_HT650_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v6", triggerresultshelper_HLT_HT650_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v7", triggerresultshelper_HLT_HT650_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v8", triggerresultshelper_HLT_HT650_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT650_v9", triggerresultshelper_HLT_HT650_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v1", triggerresultshelper_HLT_HT750_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v2", triggerresultshelper_HLT_HT750_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v3", triggerresultshelper_HLT_HT750_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v4", triggerresultshelper_HLT_HT750_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v5", triggerresultshelper_HLT_HT750_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v6", triggerresultshelper_HLT_HT750_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v7", triggerresultshelper_HLT_HT750_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v8", triggerresultshelper_HLT_HT750_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_HT750_v9", triggerresultshelper_HLT_HT750_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT650_v10", triggerresultshelper_HLT_PFHT650_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT650_v11", triggerresultshelper_HLT_PFHT650_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT650_v5", triggerresultshelper_HLT_PFHT650_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT650_v6", triggerresultshelper_HLT_PFHT650_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT650_v7", triggerresultshelper_HLT_PFHT650_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT650_v8", triggerresultshelper_HLT_PFHT650_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT650_v9", triggerresultshelper_HLT_PFHT650_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT700_v10", triggerresultshelper_HLT_PFHT700_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT700_v11", triggerresultshelper_HLT_PFHT700_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT700_v5", triggerresultshelper_HLT_PFHT700_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT700_v6", triggerresultshelper_HLT_PFHT700_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT700_v7", triggerresultshelper_HLT_PFHT700_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT700_v8", triggerresultshelper_HLT_PFHT700_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT700_v9", triggerresultshelper_HLT_PFHT700_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT750_v10", triggerresultshelper_HLT_PFHT750_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT750_v11", triggerresultshelper_HLT_PFHT750_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT750_v5", triggerresultshelper_HLT_PFHT750_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT750_v6", triggerresultshelper_HLT_PFHT750_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT750_v7", triggerresultshelper_HLT_PFHT750_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT750_v8", triggerresultshelper_HLT_PFHT750_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFHT750_v9", triggerresultshelper_HLT_PFHT750_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet320_v10", triggerresultshelper_HLT_PFJet320_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet320_v11", triggerresultshelper_HLT_PFJet320_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet320_v3", triggerresultshelper_HLT_PFJet320_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet320_v4", triggerresultshelper_HLT_PFJet320_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet320_v5", triggerresultshelper_HLT_PFJet320_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet320_v6", triggerresultshelper_HLT_PFJet320_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet320_v7", triggerresultshelper_HLT_PFJet320_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet320_v8", triggerresultshelper_HLT_PFJet320_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet320_v9", triggerresultshelper_HLT_PFJet320_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet400_v10", triggerresultshelper_HLT_PFJet400_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet400_v11", triggerresultshelper_HLT_PFJet400_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet400_v3", triggerresultshelper_HLT_PFJet400_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet400_v4", triggerresultshelper_HLT_PFJet400_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet400_v5", triggerresultshelper_HLT_PFJet400_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet400_v6", triggerresultshelper_HLT_PFJet400_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet400_v7", triggerresultshelper_HLT_PFJet400_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet400_v8", triggerresultshelper_HLT_PFJet400_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFJet400_v9", triggerresultshelper_HLT_PFJet400_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT650_v1", triggerresultshelper_HLT_PFNoPUHT650_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT650_v2", triggerresultshelper_HLT_PFNoPUHT650_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT650_v3", triggerresultshelper_HLT_PFNoPUHT650_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT650_v4", triggerresultshelper_HLT_PFNoPUHT650_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT650_v5", triggerresultshelper_HLT_PFNoPUHT650_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT650_v6", triggerresultshelper_HLT_PFNoPUHT650_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT700_v1", triggerresultshelper_HLT_PFNoPUHT700_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT700_v2", triggerresultshelper_HLT_PFNoPUHT700_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT700_v3", triggerresultshelper_HLT_PFNoPUHT700_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT700_v4", triggerresultshelper_HLT_PFNoPUHT700_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT700_v5", triggerresultshelper_HLT_PFNoPUHT700_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT700_v6", triggerresultshelper_HLT_PFNoPUHT700_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT750_v1", triggerresultshelper_HLT_PFNoPUHT750_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT750_v2", triggerresultshelper_HLT_PFNoPUHT750_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT750_v3", triggerresultshelper_HLT_PFNoPUHT750_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT750_v4", triggerresultshelper_HLT_PFNoPUHT750_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT750_v5", triggerresultshelper_HLT_PFNoPUHT750_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.HLT_PFNoPUHT750_v6", triggerresultshelper_HLT_PFNoPUHT750_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v10", triggerresultshelper_prescaleHLT_DiPFJetAve320_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v11", triggerresultshelper_prescaleHLT_DiPFJetAve320_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v12", triggerresultshelper_prescaleHLT_DiPFJetAve320_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v2", triggerresultshelper_prescaleHLT_DiPFJetAve320_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v3", triggerresultshelper_prescaleHLT_DiPFJetAve320_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v4", triggerresultshelper_prescaleHLT_DiPFJetAve320_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v5", triggerresultshelper_prescaleHLT_DiPFJetAve320_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v6", triggerresultshelper_prescaleHLT_DiPFJetAve320_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v7", triggerresultshelper_prescaleHLT_DiPFJetAve320_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v8", triggerresultshelper_prescaleHLT_DiPFJetAve320_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve320_v9", triggerresultshelper_prescaleHLT_DiPFJetAve320_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v10", triggerresultshelper_prescaleHLT_DiPFJetAve400_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v11", triggerresultshelper_prescaleHLT_DiPFJetAve400_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v12", triggerresultshelper_prescaleHLT_DiPFJetAve400_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v2", triggerresultshelper_prescaleHLT_DiPFJetAve400_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v3", triggerresultshelper_prescaleHLT_DiPFJetAve400_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v4", triggerresultshelper_prescaleHLT_DiPFJetAve400_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v5", triggerresultshelper_prescaleHLT_DiPFJetAve400_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v6", triggerresultshelper_prescaleHLT_DiPFJetAve400_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v7", triggerresultshelper_prescaleHLT_DiPFJetAve400_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v8", triggerresultshelper_prescaleHLT_DiPFJetAve400_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_DiPFJetAve400_v9", triggerresultshelper_prescaleHLT_DiPFJetAve400_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v11", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v12", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v12);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v6", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v7", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v8", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v9", triggerresultshelper_prescaleHLT_FatDiPFJetMass750_DR1p1_Deta1p5_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT450_v1", triggerresultshelper_prescaleHLT_HT450_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT450_v2", triggerresultshelper_prescaleHLT_HT450_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT450_v3", triggerresultshelper_prescaleHLT_HT450_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT450_v4", triggerresultshelper_prescaleHLT_HT450_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT450_v5", triggerresultshelper_prescaleHLT_HT450_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT450_v6", triggerresultshelper_prescaleHLT_HT450_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT450_v7", triggerresultshelper_prescaleHLT_HT450_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT450_v8", triggerresultshelper_prescaleHLT_HT450_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT450_v9", triggerresultshelper_prescaleHLT_HT450_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT500_v1", triggerresultshelper_prescaleHLT_HT500_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT500_v2", triggerresultshelper_prescaleHLT_HT500_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT500_v3", triggerresultshelper_prescaleHLT_HT500_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT500_v4", triggerresultshelper_prescaleHLT_HT500_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT500_v5", triggerresultshelper_prescaleHLT_HT500_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT500_v6", triggerresultshelper_prescaleHLT_HT500_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT500_v7", triggerresultshelper_prescaleHLT_HT500_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT500_v8", triggerresultshelper_prescaleHLT_HT500_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT500_v9", triggerresultshelper_prescaleHLT_HT500_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT550_v1", triggerresultshelper_prescaleHLT_HT550_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT550_v2", triggerresultshelper_prescaleHLT_HT550_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT550_v3", triggerresultshelper_prescaleHLT_HT550_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT550_v4", triggerresultshelper_prescaleHLT_HT550_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT550_v5", triggerresultshelper_prescaleHLT_HT550_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT550_v6", triggerresultshelper_prescaleHLT_HT550_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT550_v7", triggerresultshelper_prescaleHLT_HT550_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT550_v8", triggerresultshelper_prescaleHLT_HT550_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT550_v9", triggerresultshelper_prescaleHLT_HT550_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT650_v1", triggerresultshelper_prescaleHLT_HT650_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT650_v2", triggerresultshelper_prescaleHLT_HT650_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT650_v3", triggerresultshelper_prescaleHLT_HT650_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT650_v4", triggerresultshelper_prescaleHLT_HT650_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT650_v5", triggerresultshelper_prescaleHLT_HT650_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT650_v6", triggerresultshelper_prescaleHLT_HT650_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT650_v7", triggerresultshelper_prescaleHLT_HT650_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT650_v8", triggerresultshelper_prescaleHLT_HT650_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT650_v9", triggerresultshelper_prescaleHLT_HT650_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT750_v1", triggerresultshelper_prescaleHLT_HT750_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT750_v2", triggerresultshelper_prescaleHLT_HT750_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT750_v3", triggerresultshelper_prescaleHLT_HT750_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT750_v4", triggerresultshelper_prescaleHLT_HT750_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT750_v5", triggerresultshelper_prescaleHLT_HT750_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT750_v6", triggerresultshelper_prescaleHLT_HT750_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT750_v7", triggerresultshelper_prescaleHLT_HT750_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT750_v8", triggerresultshelper_prescaleHLT_HT750_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_HT750_v9", triggerresultshelper_prescaleHLT_HT750_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT650_v10", triggerresultshelper_prescaleHLT_PFHT650_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT650_v11", triggerresultshelper_prescaleHLT_PFHT650_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT650_v5", triggerresultshelper_prescaleHLT_PFHT650_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT650_v6", triggerresultshelper_prescaleHLT_PFHT650_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT650_v7", triggerresultshelper_prescaleHLT_PFHT650_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT650_v8", triggerresultshelper_prescaleHLT_PFHT650_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT650_v9", triggerresultshelper_prescaleHLT_PFHT650_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT700_v10", triggerresultshelper_prescaleHLT_PFHT700_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT700_v11", triggerresultshelper_prescaleHLT_PFHT700_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT700_v5", triggerresultshelper_prescaleHLT_PFHT700_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT700_v6", triggerresultshelper_prescaleHLT_PFHT700_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT700_v7", triggerresultshelper_prescaleHLT_PFHT700_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT700_v8", triggerresultshelper_prescaleHLT_PFHT700_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT700_v9", triggerresultshelper_prescaleHLT_PFHT700_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT750_v10", triggerresultshelper_prescaleHLT_PFHT750_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT750_v11", triggerresultshelper_prescaleHLT_PFHT750_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT750_v5", triggerresultshelper_prescaleHLT_PFHT750_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT750_v6", triggerresultshelper_prescaleHLT_PFHT750_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT750_v7", triggerresultshelper_prescaleHLT_PFHT750_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT750_v8", triggerresultshelper_prescaleHLT_PFHT750_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFHT750_v9", triggerresultshelper_prescaleHLT_PFHT750_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet320_v10", triggerresultshelper_prescaleHLT_PFJet320_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet320_v11", triggerresultshelper_prescaleHLT_PFJet320_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet320_v3", triggerresultshelper_prescaleHLT_PFJet320_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet320_v4", triggerresultshelper_prescaleHLT_PFJet320_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet320_v5", triggerresultshelper_prescaleHLT_PFJet320_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet320_v6", triggerresultshelper_prescaleHLT_PFJet320_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet320_v7", triggerresultshelper_prescaleHLT_PFJet320_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet320_v8", triggerresultshelper_prescaleHLT_PFJet320_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet320_v9", triggerresultshelper_prescaleHLT_PFJet320_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet400_v10", triggerresultshelper_prescaleHLT_PFJet400_v10);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet400_v11", triggerresultshelper_prescaleHLT_PFJet400_v11);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet400_v3", triggerresultshelper_prescaleHLT_PFJet400_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet400_v4", triggerresultshelper_prescaleHLT_PFJet400_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet400_v5", triggerresultshelper_prescaleHLT_PFJet400_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet400_v6", triggerresultshelper_prescaleHLT_PFJet400_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet400_v7", triggerresultshelper_prescaleHLT_PFJet400_v7);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet400_v8", triggerresultshelper_prescaleHLT_PFJet400_v8);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFJet400_v9", triggerresultshelper_prescaleHLT_PFJet400_v9);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT650_v1", triggerresultshelper_prescaleHLT_PFNoPUHT650_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT650_v2", triggerresultshelper_prescaleHLT_PFNoPUHT650_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT650_v3", triggerresultshelper_prescaleHLT_PFNoPUHT650_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT650_v4", triggerresultshelper_prescaleHLT_PFNoPUHT650_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT650_v5", triggerresultshelper_prescaleHLT_PFNoPUHT650_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT650_v6", triggerresultshelper_prescaleHLT_PFNoPUHT650_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT700_v1", triggerresultshelper_prescaleHLT_PFNoPUHT700_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT700_v2", triggerresultshelper_prescaleHLT_PFNoPUHT700_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT700_v3", triggerresultshelper_prescaleHLT_PFNoPUHT700_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT700_v4", triggerresultshelper_prescaleHLT_PFNoPUHT700_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT700_v5", triggerresultshelper_prescaleHLT_PFNoPUHT700_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT700_v6", triggerresultshelper_prescaleHLT_PFNoPUHT700_v6);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT750_v1", triggerresultshelper_prescaleHLT_PFNoPUHT750_v1);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT750_v2", triggerresultshelper_prescaleHLT_PFNoPUHT750_v2);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT750_v3", triggerresultshelper_prescaleHLT_PFNoPUHT750_v3);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT750_v4", triggerresultshelper_prescaleHLT_PFNoPUHT750_v4);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT750_v5", triggerresultshelper_prescaleHLT_PFNoPUHT750_v5);
  stream.select("edmTriggerResultsHelper_TriggerResults_HLT.prescaleHLT_PFNoPUHT750_v6", triggerresultshelper_prescaleHLT_PFNoPUHT750_v6);
  stream.select("recoVertex_goodOfflinePrimaryVertices.chi2", vertex_chi2);
  stream.select("recoVertex_goodOfflinePrimaryVertices.isFake", vertex_isFake);
  stream.select("recoVertex_goodOfflinePrimaryVertices.ndof", vertex_ndof);
  stream.select("recoVertex_goodOfflinePrimaryVertices.position_Rho", vertex_position_Rho);
  stream.select("recoVertex_goodOfflinePrimaryVertices.x", vertex_x);
  stream.select("recoVertex_goodOfflinePrimaryVertices.xError", vertex_xError);
  stream.select("recoVertex_goodOfflinePrimaryVertices.y", vertex_y);
  stream.select("recoVertex_goodOfflinePrimaryVertices.yError", vertex_yError);
  stream.select("recoVertex_goodOfflinePrimaryVertices.z", vertex_z);
  stream.select("recoVertex_goodOfflinePrimaryVertices.zError", vertex_zError);

}
//-----------------------------------------------------------------------------
// -- Utilities
//-----------------------------------------------------------------------------
void
error(std::string message)
{
  std::cout << "** error ** " << message << std::endl;
  exit(0);
}

std::string 
strip(std::string line)
{
  int l = line.size();
  if ( l == 0 ) return std::string("");
  int n = 0;
  while (((line[n] == 0)    ||
	  (line[n] == ' ' ) ||
	  (line[n] == '\n') ||
	  (line[n] == '\t')) && n < l) n++;

  int m = l-1;
  while (((line[m] == 0)    ||
	  (line[m] == ' ')  ||
	  (line[m] == '\n') ||
	  (line[m] == '\t')) && m > 0) m--;
  return line.substr(n,m-n+1);
}

std::string
nameonly(std::string filename)
{
  int i = filename.rfind("/");
  int j = filename.rfind(".");
  if ( j < 0 ) j = filename.size();
  return filename.substr(i+1,j-i-1);
}
//-----------------------------------------------------------------------------
struct outputFile
{
  outputFile(std::string filename)
   : filename_(filename),
	 file_(new TFile(filename_.c_str(), "recreate")),
	 tree_(0),
	 b_weight_(0),
	 entry_(0),
	 SAVECOUNT_(50000)
  {
	file_->cd();
	hist_ = new TH1D("counts", "", 1,0,1);
	hist_->SetBit(TH1::kCanRebin);
	hist_->SetStats(0);
  }

  outputFile(std::string filename, itreestream& stream, int savecount=50000) 
   : filename_(filename),
	 file_(new TFile(filename.c_str(), "recreate")),
	 tree_(stream.tree()->CloneTree(0)),
	 b_weight_(tree_->Branch("eventWeight", &weight_, "eventWeight/D")),
	 entry_(0),
	 SAVECOUNT_(savecount)
  {
	std::cout << "events will be skimmed to file "
			  << filename_ << std::endl;
	file_->cd();
	hist_ = new TH1D("counts", "", 1,0,1);
	hist_->SetBit(TH1::kCanRebin);
	hist_->SetStats(0);
  }

  void addEvent(double weight=1)
  {
    if ( tree_ == 0 ) return;
	
    weight_ = weight;
	file_   = tree_->GetCurrentFile();
	file_->cd();
	tree_->Fill();

	entry_++;
	if ( entry_ % SAVECOUNT_ == 0 )
	  tree_->AutoSave("SaveSelf");
  }

  void count(std::string cond, double w=1.0)
  {
    hist_->Fill(cond.c_str(), w);
  }
  
  void close()
  {
  	std::cout << "==> histograms saved to file " << filename_ << std::endl;
    if ( tree_ != 0 )
	  {
	    std::cout << "==> events skimmed to file " << filename_ << std::endl;
	    file_ = tree_->GetCurrentFile();
	  }
	file_->cd();
	//file_->Write("", TObject::kWriteDelete);
	file_->Write();
	file_->ls();
	file_->Close();
  }

  std::string filename_;  
  TFile* file_;
  TTree* tree_;
  TH1D*  hist_;
  TBranch* b_weight_;
  double     weight_;
  int    entry_;
  int    SAVECOUNT_;
};

struct commandLine
{
  std::string progname;
  std::string filelist;
  double xsect;
  double totweight;
  double lumi;
  std::string outputfilename;
  std::string systfilename;
};


void
decodeCommandLine(int argc, char** argv, commandLine& cl)
{
  cl.progname = std::string(argv[0]);

  // 1st (optional) argument
  if ( argc > 1 )
        cl.filelist = std::string(argv[1]);
  else
      	cl.filelist = std::string("filelist.txt");

  // 2nd (optional) argument
  if ( argc > 2 )
        cl.xsect = atof(argv[2]);
  else
      	cl.xsect = -1;

  // 3rd (optional) argument
  if ( argc > 3 )
        cl.totweight = atof(argv[3]);
  else
      	cl.totweight = -1;

  // 4th (optional) argument
  if ( argc > 4 )
        cl.lumi = atof(argv[4]);
  else
      	cl.lumi = -1;

  // 5th (optional) command line argument
  if ( argc > 5 )
        cl.outputfilename = std::string(argv[5]);
  else
        cl.outputfilename = cl.progname + std::string("_histograms");

  // 6th (optional) command line argument
  if ( argc > 6 )
        cl.systfilename = std::string(argv[6]);
  else
        cl.systfilename = std::string("systematics.txt");

  // Make sure extension is ".root"
  std::string name = cl.outputfilename;
  if ( name.substr(name.size()-5, 5) != std::string(".root") )
    cl.outputfilename += std::string(".root");
}

// Read ntuple filenames from file list

std::vector<std::string>
getFilenames(std::string filelist)
{
  std::ifstream stream(filelist.c_str());
  if ( !stream.good() ) error("unable to open file: " + filelist);

  // Get list of ntuple files to be processed

  std::vector<std::string> v;
  std::string filename;
  while ( stream >> filename )
	if ( strip(filename) != "" ) v.push_back(filename);
  return v;
}

#endif
