# Script to get the pileup weights

from ROOT import *

gStyle.SetOptStat(0)
gStyle.SetTextFont(42)
    
fdata = TFile.Open("data_pileup.root")
fmc = TFile.Open("mc_pileup.root")

fout = TFile.Open("pileup_weights.root","RECREATE")
fout.cd()

# Pileup profile in data obtained by:
# pileupCalc.py -i Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt -o data_pileup.root --inputLumiJSON=pileup_latest.txt --calcMode=true --minBiasXsec=69400 --maxPileupBin=100
# inputfiles were taken from the official CMS location

hdata = fdata.Get("pileup")
hdata.Sumw2()
hdata.Rebin(10)
hdata.Scale(1./hdata.Integral())
hdata.GetXaxis().SetTitle("Number of interactions")
hdata.SetLineColor(kBlack)

# Pileup distribution in MC obtained from the variable
# pileupsummaryinfo[0].getTrueNumInteractions

hmc = fmc.Get("h_pileup")
hmc.Sumw2()
hmc.Rebin(10) # something funny going on with the mc histograms
hmc.Scale(1./hmc.Integral())
hmc.GetXaxis().SetTitle("Number of interactions")
hmc.SetLineColor(kRed+2)

# Make the ratio to use for reweighting
hratio = hdata.Clone("pileup_weight")
hratio.Sumw2()
hratio.Divide(hmc)
hratio.Write()
hratio.SetLineColor(kBlue)

hratio2 = hratio.Clone() # to plot the errors
hratio2.SetFillColor(kCyan)

legend = TLegend(0.7,0.7,0.85,0.85)
legend.SetFillColor(0)
legend.SetBorderSize(0)
legend.AddEntry(hdata,"Data","l")
legend.AddEntry(hmc,"MC","l")

# Also make sanity plot
c = TCanvas("c","c",1600,1200)
c.Divide(2,2)
c.cd(1)
hdata.Draw("HIST")
hmc.Draw("HISTsame")
legend.Draw("same")
c.cd(2)
c.cd(2).SetLogy()
hdata.Draw("HIST")
hmc.Draw("HISTsame")
legend.Draw("same")
c.cd(3)
hratio2.Draw("E2")
hratio.Draw("HISTsame")
c.cd(4)
c.cd(4).SetLogy()
hratio2.Draw("E2")
hratio.Draw("HISTsame")
c.cd()
c.SaveAs("compare_pileup_profile.pdf")
#c.SaveAs("compare_pileup_profile.root")

fout.Close()

fdata.Close()
fmc.Close()
