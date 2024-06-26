import ROOT, sys, os, time, re
import numpy as np
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog fileName.root BinNumber")
(opt,args) = parser.parse_args()

ROOT.gROOT.SetStyle("Plain")
#ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.SetBatch(True)

ROOT.gStyle.SetPadTopMargin(0.07);
ROOT.gStyle.SetPadBottomMargin(0.1);
ROOT.gStyle.SetPadLeftMargin(0.15);
ROOT.gStyle.SetPadRightMargin(0.13);

fileName = sys.argv[1]
BinNumber = sys.argv[2]

bin = int(BinNumber)
# bin 3: pt>60 and I_as > 0.05
# bin 25: pt>65 and I_as > 0.175
# bin 28: pt>65 and I_as > 0.3

print("Filename: "+fileName)
input_file = fileName

ProjBin = int(BinNumber)
newFileDir = fileName[0:-5] + "_Bin" + str(ProjBin)


f = ROOT.TFile.Open(input_file)
fileOut = open("SignalBackgroundEff.txt", "a")

isData = False
if ("SingleMuon" in fileName) : isData = True

iDontWannaRunPlots = False

dirs = []
for i in range(0, f.GetListOfKeys().GetEntries()):
  # Remove/modify unnecessary stuff from the name of the plot that was required by SmartHistos to ditinguish plots
  dirname = f.GetListOfKeys().At(i).GetName()
  curr_dir = f.GetDirectory(dirname)
# print("dirname: "+dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      # Match the plot of interest
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = f.GetDirectory(dirname+"/"+keyname)
#                    print("keyname: "+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          if ("__" in keyname2) : continue
          # The plot should be TCanvas
          newname = dirname + "/" + keyname+ "/" + keyname2
#          print("newname: "+newname)
          obj = f.Get(newname)
          
          
          tex2 = ROOT.TLatex(0.13,0.94,"CMS");
          #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
          tex2.SetNDC();
          tex2.SetTextFont(61);
          tex2.SetTextSize(0.0675);
          tex2.SetLineWidth(2);

          #tex3 = ROOT.TLatex(0.27,0.96,"Simulation"); # for square plots
          #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
          tex3 = ROOT.TLatex(0.27,0.94,"Internal");
          tex3.SetNDC();
          tex3.SetTextFont(52);
          tex3.SetTextSize(0.0485);
          tex3.SetLineWidth(2);


          tex4 = ROOT.TLatex()
          if ("BefPreS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"Before pre-selection")
#            if ("BefPreS_Eta" in keyname2) :
#              print("BefPreS number of tracks in plot (" +keyname2 + ") : "+str(obj.Integral()))
          elif ("N1" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"After N-1 selection")
#            if ("N1_Eta" in keyname2) :
#              print("N-1 number of tracks in plot (" +keyname2 + ") : "+str(obj.Integral()))
          elif ("PostPreS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"After pre-selection")
#            if ("PostPreS_Eta" in keyname2) :
#              print("PostPreS number of tracks in plot (" +keyname2 + ") : "+str(obj.Integral()))
          tex4.SetNDC();
          tex4.SetTextFont(52);
          tex4.SetTextSize(0.045);
          tex4.SetLineWidth(2);
          
          codeVersion = fileName[fileName.find("CodeV")+5:fileName.find("CodeV")+9]
          fileVersion = fileName[fileName.find("2018")+5:fileName.find("CodeV")+9]
          tex5 = ROOT.TLatex(0.07,0.03,fileVersion);
          tex5.SetNDC();
          tex5.SetTextFont(52);
          tex5.SetTextSize(0.0185);
          tex5.SetLineWidth(2);
          
          if (keyname2.find("Vs")==-1) :
            axisXTitle = keyname2[keyname2.find("_")+1:]
            axisYTitle = "Yield/bin"
          else :
            axisXTitle = keyname2[keyname2.find("_")+1:keyname2.find("Vs")]
            axisYTitle = keyname2[keyname2.find("Vs")+2:]
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if obj.InheritsFrom("TObject"):
              can = obj
              obj.SetStats(0)
              can = ROOT.TCanvas(newname,newname,800,800)
              # Name of the png to be saved
              name = fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  ".png"
              if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
              if (obj.GetEntries() == 0 ) :
                continue
                
#             when I dont want to plot everything
              if (iDontWannaRunPlots) : continue
                
              if ("Gen" in keyname2 and isData) : continue
#                 print(obj.ClassName())
              if (obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D"):
                obj.SetTitle("")
                if ("ProbQVsProbXY" in keyname2) :
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(3.22),-1)
                  obj.GetYaxis().SetRange(obj.GetYaxis().FindBin(0.0),obj.GetYaxis().FindBin(0.1))
                  obj.Project3D("YZ").Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIh.png")
                  obj.GetXaxis().UnZoom()
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(3.22))
                  obj.GetYaxis().SetRange(obj.GetYaxis().FindBin(0.0),obj.GetYaxis().FindBin(0.1))
                  obj.Project3D("YZ").Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIh.png")
                  obj.GetXaxis().UnZoom()
                  obj.Project3D("YZ").Draw("COLZ")
                  can.SaveAs(name)
                # this maybe should go from bin to bin+1 ?
                if ("IhVsLayer" in keyname2 or "IhVsLayer" in keyname2) :
                  obj.SetMarkerStyle(20)
                  if ("PostPreS_IasPixelIhVsLayer" in keyname2) :
                      ratioOfhighIhPartOlowIh = ROOT.TH1F("RatioOfHighIhOverLowIh",";;Ratio of Ih>5 over Ih<5",7,0.,7.)
                      for x in range(1,obj.GetNbinsX()) :
                        lowIhPart = obj.Project3D("YZ").Integral(x,x,1,obj.Project3D("YZ").GetYaxis().FindBin(5.0))
                        highIhPart = obj.Project3D("YZ").Integral(x,x,obj.Project3D("YZ").GetYaxis().FindBin(5.0),obj.Project3D("YZ").GetNbinsY())
                        ratio = 0
                        if (lowIhPart>0) :
                          ratio = highIhPart/lowIhPart
                        ratioOfhighIhPartOlowIh.SetBinContent(x,ratio)
                      ratioOfhighIhPartOlowIh.Draw("COLZ")
                      ratioOfhighIhPartOlowIh.SetStats(0)
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(1,"BPix L1")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(2,"BPix L2")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(3,"BPix L3")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(4,"BPix L4")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(5,"FPix D1")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(6,"FPix D2")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(7,"FPix D3")
                      ratioOfhighIhPartOlowIh.GetYaxis().SetTitleOffset(1.5)
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(1,"BPix L1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(2,"BPix L2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(3,"BPix L3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(4,"BPix L4")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(5,"FPix D1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(6,"FPix D2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(7,"FPix D3")
                      obj.Project3D("YZ").GetYaxis().SetTitleOffset(0.8)
                      obj.Project3D("YZ").GetYaxis().SetTitle("Ih")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(1,"BPix L1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(2,"BPix L2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(3,"BPix L3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(4,"BPix L4")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(5,"FPix D1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(6,"FPix D2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(7,"FPix D3")
                      obj.Project3D("XZ").GetYaxis().SetTitle("Ias")
                      obj.Project3D("XZ").GetYaxis().SetTitleOffset(0.8)
                  elif ("PostPreS_IasStripIhVsLayer" in keyname2) :
                      ratioOfhighIhPartOlowIh = ROOT.TH1F("RatioOfHighIhOverLowIh",";;Ratio of Ih>5 over Ih<5",23,0.,23.)
                      for x in range(1,obj.GetNbinsX()) :
                        lowIhPart = obj.Project3D("YZ").Integral(x,x,1,obj.Project3D("YZ").GetYaxis().FindBin(5.0))
                        highIhPart = obj.Project3D("YZ").Integral(x,x,obj.Project3D("YZ").GetYaxis().FindBin(5.0),obj.Project3D("YZ").GetNbinsY())
                        ratio = 0
                        if (lowIhPart>0) :
                          ratio = highIhPart/lowIhPart
                        ratioOfhighIhPartOlowIh.SetBinContent(x,ratio)
                      ratioOfhighIhPartOlowIh.Draw("COLZ")
                      ratioOfhighIhPartOlowIh.SetStats(0)
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(1,"TIB L1")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(2,"TIB L2")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(3,"TIB L3")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(4,"TIB L4")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(5,"TOB L1")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(6,"TOB L2")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(7,"TOB L3")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(8,"TOB L4")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(9,"TOB L5")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(10,"TOB L6")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(11,"TID D1")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(12,"TID D2")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(13,"TID D3")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(14,"TEC D1")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(15,"TEC D2")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(16,"TEC D3")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(17,"TEC D4")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(18,"TEC D5")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(19,"TEC D6")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(20,"TEC D7")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(21,"TEC D8")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(22,"TEC D9")
                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(23,"TEC D10")
                      ratioOfhighIhPartOlowIh.GetYaxis().SetTitleOffset(1.5)
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(1,"TIB L1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(2,"TIB L2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(3,"TIB L3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(4,"TIB L4")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(5,"TOB L1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(6,"TOB L2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(7,"TOB L3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(8,"TOB L4")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(9,"TOB L5")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(10,"TOB L6")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(11,"TID D1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(12,"TID D2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(13,"TID D3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(14,"TEC D1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(15,"TEC D2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(16,"TEC D3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(17,"TEC D4")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(18,"TEC D5")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(19,"TEC D6")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(20,"TEC D7")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(21,"TEC D8")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(22,"TEC D9")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(23,"TEC D10")
                      obj.Project3D("YZ").GetYaxis().SetTitleOffset(0.8)
                      obj.Project3D("YZ").GetYaxis().SetTitle("Ih")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(1,"TIB L1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(2,"TIB L2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(3,"TIB L3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(4,"TIB L4")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(5,"TOB L1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(6,"TOB L2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(7,"TOB L3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(8,"TOB L4")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(9,"TOB L5")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(10,"TOB L6")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(11,"TID D1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(12,"TID D2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(13,"TID D3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(14,"TEC D1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(15,"TEC D2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(16,"TEC D3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(17,"TEC D4")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(18,"TEC D5")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(19,"TEC D6")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(20,"TEC D7")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(21,"TEC D8")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(22,"TEC D9")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(23,"TEC D10")
                      obj.Project3D("XZ").GetYaxis().SetTitle("Ias")
                      obj.Project3D("XZ").GetYaxis().SetTitleOffset(0.8)
#                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(1.0))
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_RatioOfLowIasHighIas.png")
                  obj.Project3D("YZ").SetStats(0)
                  obj.Project3D("YZ").Draw("COLZ")
                  obj.Project3D("YZ").SetTitle("")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_allIas.png")
                  obj.Project3D("XZ").SetStats(0)
                  obj.Project3D("XZ").Draw("COLZ")
                  obj.Project3D("XZ").SetTitle("")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_allIh.png")
                  
                  obj.Project3D("YZ").GetYaxis().UnZoom()
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.7),obj.GetXaxis().FindBin(1.0))
                  obj.Project3D("YZ").Draw("COLZ")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIas.png")
                  
                  obj.GetXaxis().UnZoom()
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(0.7))
                  projObj = obj.Project3D("YZ")
                  if (projObj.GetEntries()==0) : continue
                  projObj.Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIas.png")
                  if ("PostPreS_IasPixelIhVsLayer" in keyname2):
                    for i in range(7) :
                      obj.Project3D("YZ").ProjectionY(newname,i,i+1,"e").Draw()
                      can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_PixLayer"+str(i)+".png")
                  elif ("PostPreS_IasStripIhVsLayer" in keyname2):
                    for i in range(23) :
                      obj.Project3D("YZ").ProjectionY(newname,i,i+1,"e").Draw()
                      can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_StripLayer"+str(i)+".png")
                else :
                  obj.SetMarkerStyle(20)
                  obj.GetXaxis().SetRange(bin,bin)
                  obj.Project3D("ZY").Draw("COLZ")
              elif ((obj.ClassName() == "TH2F" or obj.ClassName() == "TH2D") and not (keyname2 == "GenPtVsRecoPt" or "PreS_" in keyname2 or keyname2 == "CutFlowEta" or  keyname2 == "CutFlowPfType"  or "N1_" in keyname2)):
                obj.SetTitle("")
                obj.SetMarkerStyle(20)
                obj.ProjectionY(newname,bin,bin,"e").Draw("COLZ")
              if ((obj.ClassName() == "TH2F") and ("Clu" in keyname2)) :
                profYobj = obj.ProfileY()
                profYobj.GetYaxis().SetTitle(axisXTitle)
                profYobj.GetYaxis().SetTitleOffset(1.3)
                profYobj.GetYaxis().SetLabelSize(0.03)
                profYobj.SetStats(0)
                profYobj.DrawClone("COLZ")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_profileY.png")
                
                if ("CluSpecInCPE" in keyname2) :
                  obj.SetStats(0)
                  obj.GetXaxis().SetTitle("")
                  obj.GetXaxis().SetBinLabel(1,"isOnEdge")
                  obj.GetXaxis().SetBinLabel(2,"hasBadPixels")
                  obj.GetXaxis().SetBinLabel(3,"spansTwoROCs")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  obj.Draw("COLZ L")
                  obj.GetXaxis().SetTitle(axisXTitle)
                obj.DrawClone("COLZ L")
                obj.GetYaxis().SetTitle(axisYTitle)
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(name)
              elif ("GenPtVsRecoPt" in keyname2) :
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.Draw("COLZ")
                print(str(keyname2) + ": " + str(round(obj.GetCorrelationFactor(),2)))
              elif ("IasVs" in keyname2) :
                obj.SetMarkerStyle(20)
                if ("Angle" in keyname2 or "NumSibling" in keyname2) :
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  obj.ProjectionY(newname,obj.GetXaxis().FindBin(0.7),obj.GetNbinsX(),"e").Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIas.png")
                  projObjs = obj.ProjectionY(newname,1,obj.GetXaxis().FindBin(0.7),"e")
                  if (projObjs.GetEntries()==0) : continue
                  projObjs.Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIas.png")
                  obj.ProjectionY(newname,1,obj.GetNbinsX(),"e").Draw("COLZ")
                else :
                  projObject = obj.ProjectionY(newname+"_lowIas",1,obj.GetXaxis().FindBin(0.7),"e")
                  if (projObject.GetEntries()==0) : continue
                  myPie = ROOT.TPie(projObject)
                  myPie.SetLabelFormat("%txt (%perc)")
                  myPie.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIas.png")
                  
                  objProj = obj.ProjectionY(newname+"_highIas",obj.GetXaxis().FindBin(0.7),obj.GetNbinsX(),"e")
                  if (objProj.GetEntries()==0) : continue
                  myPie2 = ROOT.TPie(objProj)
                  myPie2.SetLabelFormat("%txt (%perc)")
                  myPie2.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie2.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIas.png")
                  
                  myPie3 = ROOT.TPie(obj.ProjectionY(newname,1,obj.GetNbinsX(),"e"))
                  myPie3.SetLabelFormat("%txt (%perc)")
                  myPie3.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie3.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  ".png")
              elif ("EoPVs" in keyname2) :
                obj.SetMarkerStyle(20)
                if ("Angle" in keyname2) :
                  obj.ProjectionY(newname+"_lowEoP",1,obj.GetXaxis().FindBin(0.85),"e").Draw("COLZ")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowEoP.png")
#                  obj.ProjectionY(newname+"_highEoP",obj.GetXaxis().FindBin(0.85),obj.GetNbinsX(),"e").Draw("COLZ")
#                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highEoP.png")
                elif ("EoPVsPfType" in keyname2) :
                  obj.SetMarkerColor(1)
                  obj.SetLineColor(1)
                  obj.SetMarkerStyle(20)
                  obj.SetStats(0)
#                  obj.Scale(1/obj.GetMaximum())
                  obj.GetYaxis().SetBinLabel(1,"AllTracks")
                  obj.GetYaxis().SetBinLabel(2,"PFtracks")
                  obj.GetYaxis().SetBinLabel(3,"isElectron")
                  obj.GetYaxis().SetBinLabel(4,"isMuon")
                  obj.GetYaxis().SetBinLabel(5,"isPhoton")
                  obj.GetYaxis().SetBinLabel(6,"isChHadron")
                  obj.GetYaxis().SetBinLabel(7,"isNeutHadron")
                  obj.GetYaxis().SetBinLabel(8,"isUndefined")
                  obj.GetYaxis().SetBinLabel(9,"notPFtrack")
                  obj.GetXaxis().SetTitle("EoP")
                  obj.Draw("COLZ")
                else :
                  projObj = obj.ProjectionY(newname+"_lowEoP",1,obj.GetXaxis().FindBin(0.85),"e")
                  if (projObj.GetEntries() == 0) : continue
                  myPie = ROOT.TPie(projObj)
                  myPie.SetLabelFormat("%txt (%perc)")
                  myPie.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowEoP.png")
                  
                  projObj2 = obj.ProjectionY(newname+"_highEoP",obj.GetXaxis().FindBin(0.85),obj.GetNbinsX(),"e")
                  if (projObj2.GetEntries() == 0) : continue
                  myPie2 = ROOT.TPie(projObj2)
                  myPie2.SetLabelFormat("%txt (%perc)")
                  myPie2.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie2.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highEoP.png")
              elif (obj.ClassName() == "TH2F" and  "ProbQVsIas" in keyname2) :
                obj.GetXaxis().SetTitle(axisXTitle)
                obj.GetYaxis().SetTitle(axisYTitle)
                obj.GetYaxis().SetTitleOffset(1.3)
                obj.GetYaxis().SetLabelSize(0.03)
                obj.SetStats(0)
                obj.DrawClone("COLZ")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                print(str(keyname2) + ": " + str(round(obj.GetCorrelationFactor(),2)))
                can.SaveAs(name)
              elif (keyname2== "CutFlow") :
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
                obj.GetXaxis().SetTitle("")
                obj.GetYaxis().SetTitle("Efficiency")
                obj.GetXaxis().SetBinLabel(1,"Trigger")
                obj.GetXaxis().SetBinLabel(2,"Eta")
                obj.GetXaxis().SetBinLabel(3,"pT")
                obj.GetXaxis().SetBinLabel(4,"NumPixHits")
                obj.GetXaxis().SetBinLabel(5,"ValidFract")
                obj.GetXaxis().SetBinLabel(6,"NumDeDx")
                obj.GetXaxis().SetBinLabel(7,"HighPurity")
                obj.GetXaxis().SetBinLabel(8,"Chi2oDOF")
                obj.GetXaxis().SetBinLabel(9,"EoP")
                obj.GetXaxis().SetBinLabel(10,"dz")
                obj.GetXaxis().SetBinLabel(11,"dxy")
                obj.GetXaxis().SetBinLabel(12,"") #pTerrOverpT
                obj.GetXaxis().SetBinLabel(13,"dRminPfJet")
                obj.GetXaxis().SetBinLabel(14,"MiniIso")
                obj.GetXaxis().SetBinLabel(15,"PFid")
                obj.GetXaxis().SetBinLabel(16,"Ih")
                obj.GetXaxis().SetBinLabel(17,"ProbXY")
                obj.GetXaxis().SetBinLabel(18,"") #ProbQ
                obj.GetXaxis().SetBinLabel(19,"") #MuStat
                obj.GetXaxis().SetBinLabel(20,"") #PhiTOF
                obj.GetXaxis().SetBinLabel(21,"") #EtaTOF
                obj.Draw("COLZ L")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.GetYaxis().SetRangeUser(0.,1.3)
                can.SaveAs(name)
              elif (keyname2== "CutFlowReverse") :
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
                obj.GetXaxis().SetBinLabel(1,"Trigger")
                obj.GetXaxis().SetBinLabel(2,"Eta")
                obj.GetXaxis().SetBinLabel(3,"pT")
                obj.GetXaxis().SetBinLabel(4,"NumPixHits")
                obj.GetXaxis().SetBinLabel(5,"ValidFract")
                obj.GetXaxis().SetBinLabel(6,"NumDeDx")
                obj.GetXaxis().SetBinLabel(7,"HighPurity")
                obj.GetXaxis().SetBinLabel(8,"Chi2oDOF")
                obj.GetXaxis().SetBinLabel(9,"EoP")
                obj.GetXaxis().SetBinLabel(10,"dz")
                obj.GetXaxis().SetBinLabel(11,"dxy")
                obj.GetXaxis().SetBinLabel(12,"")
                obj.GetXaxis().SetBinLabel(13,"dRminPfJet")
                obj.GetXaxis().SetBinLabel(14,"MiniIso")
                obj.GetXaxis().SetBinLabel(15,"PFid")
                obj.GetXaxis().SetBinLabel(16,"Ih")
                obj.GetXaxis().SetBinLabel(17,"ProbXY")
                obj.GetXaxis().SetBinLabel(18,"")
                obj.GetXaxis().SetBinLabel(19,"")
                obj.GetXaxis().SetBinLabel(20,"")
                obj.GetXaxis().SetBinLabel(21,"")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.Draw("COLZ L")
              elif ("_pfType" in keyname2) :
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
                obj.GetXaxis().SetBinLabel(1,"AllTracks")
                obj.GetXaxis().SetBinLabel(2,"PFtracks")
                obj.GetXaxis().SetBinLabel(3,"isElectron")
                obj.GetXaxis().SetBinLabel(4,"isMuon")
                obj.GetXaxis().SetBinLabel(5,"isPhoton")
                obj.GetXaxis().SetBinLabel(6,"isChHadron")
                obj.GetXaxis().SetBinLabel(7,"isNeutHadron")
                obj.GetXaxis().SetBinLabel(8,"isUndefined")
                obj.GetXaxis().SetBinLabel(9,"notPFtrack")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.Draw("COLZ L")
              elif ((keyname2 == "CutFlowEta") or (keyname2 == "CutFlowProbQ") or (keyname2 == "CutFlowPfType") or (keyname2 == "CutFlowProbQ")) :
                obj.SetStats(0)
                obj.GetYaxis().SetBinLabel(1,"Trigger")
                obj.GetYaxis().SetBinLabel(2,"Eta")
                obj.GetYaxis().SetBinLabel(3,"pT")
                obj.GetYaxis().SetBinLabel(4,"NumPixHits")
                obj.GetYaxis().SetBinLabel(5,"ValidFract")
                obj.GetYaxis().SetBinLabel(6,"NumDeDx")
                obj.GetYaxis().SetBinLabel(7,"HighPurity")
                obj.GetYaxis().SetBinLabel(8,"Chi2oDOF")
                obj.GetYaxis().SetBinLabel(9,"EoP")
                obj.GetYaxis().SetBinLabel(10,"dz")
                obj.GetYaxis().SetBinLabel(11,"dxy")
                obj.GetYaxis().SetBinLabel(12,"")
                obj.GetYaxis().SetBinLabel(13,"dRminPfJet")
                obj.GetYaxis().SetBinLabel(14,"MiniIso")
                obj.GetYaxis().SetBinLabel(15,"PFid")
                obj.GetYaxis().SetBinLabel(16,"Ih")
                obj.GetYaxis().SetBinLabel(17,"ProbXY")
                obj.GetYaxis().SetBinLabel(18,"")
                obj.GetYaxis().SetBinLabel(19,"")
                obj.GetYaxis().SetBinLabel(20,"")
                obj.GetYaxis().SetBinLabel(21,"")
                obj.GetYaxis().SetTitle("")
                if (keyname2 == "CutFlowPfType"):
                  obj.Scale(1/obj.GetMaximum())
                  obj.GetXaxis().SetBinLabel(1,"AllTracks")
                  obj.GetXaxis().SetBinLabel(2,"PFtracks")
                  obj.GetXaxis().SetBinLabel(3,"isElectron")
                  obj.GetXaxis().SetBinLabel(4,"isMuon")
                  obj.GetXaxis().SetBinLabel(5,"isPhoton")
                  obj.GetXaxis().SetBinLabel(6,"isChHadron")
                  obj.GetXaxis().SetBinLabel(7,"isNeutHadron")
                  obj.GetXaxis().SetBinLabel(8,"isUndefined")
                  obj.GetXaxis().SetBinLabel(9,"notPFtrack")
                  for x in range(1,obj.GetNbinsX()) :
                    localMax = obj.GetBinContent(x,1)
                    for y in range(1,obj.GetNbinsY()) :
                      value = obj.GetBinContent(x,y)
                      if (value<=0) : continue
                      obj.SetBinContent(x,y, value/localMax)
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.Draw("COLZ")
                can.SaveAs(name)
              else :
                obj.SetMarkerStyle(20)
                obj.SetTitle("")
                obj.GetXaxis().SetTitle(axisXTitle)
                obj.GetYaxis().SetTitle(axisYTitle)
                obj.Draw("COLZ")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(name)
              
              if ("Angle" in keyname2 and obj.ClassName() == "TH2F" ) :
                obj.Draw("COLZ")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_2d.png")
              
              ## ----------------------------------------------------------------------------------------------
              # now let's plot everything in a logy scale
              obj.SetMarkerStyle(20)
              obj.SetMinimum(0.000001)
#              obj.SetMaximum(10000)
              if ("GenID" in keyname2 or obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D") :
                continue
                
              if (keyname2=="CutFlowProbQ" or keyname2=="CutFlowPfType" or keyname2=="CutFlowEta" or "Vs" in keyname2) :
#              or "PostPreS_IasPixelIhVsLayer" in keyname2 or "PostPreS_IasStripIhVsLayer" in keyname2
                obj.SetMaximum(obj.GetMaximum())
                obj.SetMinimum(0.000001)
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SetLogz()

              else :
                obj.SetMaximum(obj.GetMaximum()*100)
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
#                obj.SetMinimum(0.0001)
                can.SetLogy()

                obj.SetTitle("")
              can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_logy.png")
#              if (keyname2 == "BS_ProbQ") :
#                firstBinInt = obj.Integral(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(0.1))
#                totalInt = obj.Integral(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(1.0))
#                efficiency = firstBinInt/float(totalInt)
#                fileOut.write(name+": "+str(efficiency)+"\n")
#                fileOut.close()
              #can.SaveAs(name.replace(".png",".pdf"))
              #can.SaveAs(name.replace(".png",".C"))
              can.Close()

          else:
              print(keyname+"   "+newname + " does not inherit from TObject" )

Mass = f.Get("/analyzer/BaseName/Mass")
Mass_wPred = f.Get("/analyzer/BaseName/Pred_Mass_CB")
if Mass_wPred :
  tex5m = ROOT.TLatex(0.07,0.01,fileVersion)
  tex5m.SetNDC();
  tex5m.SetTextFont(52);
  tex5m.SetTextSize(0.0185);
  tex5m.SetLineWidth(2);
  name = fileName[0:-5] + "_Bin" + str(bin)+ "/"
  if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
  massBins = [10.,50.,100.,200.,300.,500.,1000.,4000.]
  massBinsArray = np.array(massBins)
  Mass_projY_NotRebinned = Mass.ProjectionY("Mass_projY_NotRebinned",ProjBin,ProjBin,"e")
  Mass_wPred_projY_NotRebinned = Mass_wPred.ProjectionY("Mass_wPred_projY_NotRebinned",ProjBin,ProjBin,"e")

  Mass_projY = ROOT.TH1F("Mass_projY" , "Mass_projY" , len(massBinsArray)-1, massBinsArray)
  Mass_wPred_projY = ROOT.TH1F("Mass_wPred_projY" , "Mass_wPred_projY" , len(massBinsArray)-1, massBinsArray)

  print("Mass_projY_NotRebinned.Integral(): ",Mass_projY_NotRebinned.Integral())
  print("Mass_wPred_projY_NotRebinned.Integral(): ",Mass_wPred_projY_NotRebinned.Integral())

  KSvalue = Mass_projY_NotRebinned.KolmogorovTest(Mass_wPred_projY_NotRebinned,"XD")
  print("KS-test: "+str(KSvalue))

  for i in range(1,len(massBinsArray)) :
    voltBin = Mass_projY_NotRebinned.FindBin(massBinsArray[i-1])+1
    currentBin = Mass_projY_NotRebinned.FindBin(massBinsArray[i])
    Mass_projYCont = 0.0
    Mass_wPred_projYCont = 0.0
    Mass_projYCont_err2 = 0.0
    Mass_wPred_projYCont_err2 = 0.0
    for j in range(voltBin,currentBin) :
      Mass_projYCont += Mass_projY_NotRebinned.GetBinContent(j)
      Mass_projYCont_err2 += (Mass_projY_NotRebinned.GetBinError(j) * Mass_projY_NotRebinned.GetBinError(j))
      
      Mass_wPred_projYCont += Mass_wPred_projY_NotRebinned.GetBinContent(j)
      Mass_wPred_projYCont_err2 += (Mass_wPred_projY_NotRebinned.GetBinError(j)*Mass_wPred_projY_NotRebinned.GetBinError(j))
    Mass_projY.SetBinContent(i,Mass_projYCont)
    Mass_projY.SetBinError(i,np.sqrt(Mass_projYCont_err2))
    Mass_wPred_projY.SetBinContent(i,Mass_wPred_projYCont)
    Mass_wPred_projY.SetBinError(i,np.sqrt(Mass_wPred_projYCont_err2))
    

  print("----------------------------------------------")
  KSvalue2 = Mass_projY.KolmogorovTest(Mass_wPred_projY,"XD")
  print("KS-test after rebinning: "+str(KSvalue2))

  Mass_projY.SetMarkerColor(1)
  Mass_projY.SetLineColor(1)
  Mass_projY.SetMarkerStyle(20)
  Mass_projY.SetTitle("")
  Mass_projY.GetXaxis().SetTitleSize(0.05)
  Mass_projY.GetXaxis().SetTitleOffset(1)
  Mass_projY.GetXaxis().SetTitle("Mass [GeV]")
  Mass_projY.GetYaxis().SetTitle("Tracks/bin")
  Mass_projY.GetYaxis().SetTitleSize(0.05)
  Mass_projY.GetYaxis().SetLabelSize(0.03)
  Mass_projY.GetYaxis().SetTitleOffset(1)
  Mass_projY.SetStats(0)
#  Mass_projY.GetYaxis().SetRangeUser(0.001,Mass_projY.GetMaximum())


  Mass_wPred_projY.SetMarkerColor(2)
  Mass_wPred_projY.SetLineColor(2)
  Mass_wPred_projY.SetMarkerStyle(20)
  Mass_wPred_projY.SetTitle("")
  Mass_wPred_projY.GetXaxis().SetTitleSize(0.05)
  Mass_wPred_projY.GetXaxis().SetTitleOffset(1)
  Mass_wPred_projY.GetXaxis().SetTitle("Mass [GeV]")
  Mass_wPred_projY.GetYaxis().SetTitle("Tracks/bin")
  Mass_wPred_projY.GetYaxis().SetTitleSize(0.05)
  Mass_wPred_projY.GetYaxis().SetTitleOffset(1)
  Mass_wPred_projY.GetYaxis().SetLabelSize(0.03)
  Mass_wPred_projY.SetStats(0)
  
  Mass_projY_NotRebinned.SetMarkerColor(1)
  Mass_projY_NotRebinned.SetLineColor(1)
  Mass_projY_NotRebinned.SetMarkerStyle(20)
  Mass_projY_NotRebinned.SetTitle("")
  Mass_projY_NotRebinned.GetXaxis().SetTitleSize(0.05)
  Mass_projY_NotRebinned.GetXaxis().SetTitleOffset(1)
  Mass_projY_NotRebinned.GetXaxis().SetTitle("Mass [GeV]")
  Mass_projY_NotRebinned.GetYaxis().SetTitle("Tracks/bin")
  Mass_projY_NotRebinned.GetYaxis().SetTitleSize(0.05)
  Mass_projY_NotRebinned.GetYaxis().SetTitleOffset(1)
  Mass_projY_NotRebinned.SetStats(0)
#  Mass_projY_NotRebinned.GetYaxis().SetRangeUser(0.001,Mass_projY_NotRebinned.GetMaximum())


  Mass_wPred_projY_NotRebinned.SetMarkerColor(2)
  Mass_wPred_projY_NotRebinned.SetLineColor(2)
  Mass_wPred_projY_NotRebinned.SetMarkerStyle(20)
  Mass_wPred_projY_NotRebinned.SetTitle("")
  Mass_wPred_projY_NotRebinned.GetXaxis().SetTitleSize(0.05)
  Mass_wPred_projY_NotRebinned.GetXaxis().SetTitleOffset(1)
  Mass_wPred_projY_NotRebinned.GetXaxis().SetTitle("Mass [GeV]")
  Mass_wPred_projY_NotRebinned.GetYaxis().SetTitle("Tracks/bin")
  Mass_wPred_projY_NotRebinned.GetYaxis().SetTitleSize(0.05)
  Mass_wPred_projY_NotRebinned.GetYaxis().SetTitleOffset(1)
  Mass_wPred_projY_NotRebinned.SetStats(0)


  print("Mass_projY.Integral(): ",Mass_projY.Integral())
  print("Mass_wPred_projY.Integral(): ",Mass_wPred_projY.Integral())

  legMass =  ROOT.TLegend(.45,.75,.80,.9,"","brNDC")
  legMass.SetTextFont(42)
  legMass.SetTextSize(0.035)
  legMass.SetBorderSize(1);
  legMass.SetLineColor(1);
  legMass.SetLineStyle(1);
  legMass.SetLineWidth(1);
  legMass.SetFillColor(0);
  legMass.SetFillStyle(1001);
  legMass.AddEntry(Mass_wPred_projY,"Prediction","LP")
  legMass.AddEntry(Mass_projY,"Observation","LP")

  tex4 = ROOT.TLatex(0.7,0.93,"K-S test v2: "+str(round(KSvalue2,4)));
  tex4.SetNDC();
  tex4.SetTextFont(52);
  tex4.SetTextSize(0.0485);
  tex4.SetLineWidth(2);

  cMass_projY = ROOT.TCanvas('cMass_projY', 'cMass_projY',800,800)

  rp = ROOT.TRatioPlot(Mass_projY,Mass_wPred_projY, "diffsigerrasym")

  rp.SetH1DrawOpt("P");
  rp.SetH2DrawOpt("P");

  rp.Draw()
  #rp.GetUpperPad().BuildLegend()
  rp.SetLeftMargin(0.13);
  rp.SetRightMargin(0.05);
  rp.SetUpTopMargin(0.1);
  rp.SetLowTopMargin(0.02);
  rp.SetLowBottomMargin(0.35);

  max = Mass_projY.GetMaximum()*1.2
  Mass_projY.SetMaximum(max);
  rp.GetLowerRefGraph().SetMinimum(-4);
  rp.GetLowerRefGraph().SetMaximum(4);
  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
  #rp.GetLowerRefGraph().SetLineColor(0) #0
  rp.GetLowerRefGraph().SetMarkerStyle(20)
  rp.GetLowerRefGraph().SetMarkerSize(1);
  rp.GetLowYaxis().SetNdivisions(505);
  rp.GetLowerRefYaxis().SetTitle("Pull");
  rp.GetLowerRefYaxis().SetTitleSize(0.05);
  rp.GetLowerRefYaxis().SetTitleOffset(1);
  rp.GetLowerRefYaxis().SetLabelSize(0.035);


  rp.GetLowerRefXaxis().SetTitleSize(0.05);
  rp.GetLowerRefXaxis().SetTitleOffset(0.8);
  rp.GetLowerRefXaxis().SetLabelSize(0.035);
  cMass_projY.Modified()
  cMass_projY.Update()
  #Mass_projY.Draw()
  #Mass_wPred_projY.Draw("SAME")
  #rp.Draw("X")

  rp.GetUpperPad().cd();
  legMass.Draw("SAME")
  tex2.Draw("SAME")
  tex3.Draw("SAME")
  tex4.Draw("SAME")
  tex5m.Draw("SAME")
  
  name = newFileDir + "/cMass.png"
  cMass_projY.SaveAs(name)
  
  
  #############################################################################
  cMass_projY_log = ROOT.TCanvas('cMass_projY_log', 'cMass_projY_log',800,800)
  cMass_projY_log.SetLogy()

  rp2 = ROOT.TRatioPlot(Mass_projY,Mass_wPred_projY, "diffsigerrasym")

  rp2.SetH1DrawOpt("P");
  rp2.SetH2DrawOpt("P");

  rp2.Draw()
  #rp2.GetUpperPad().BuildLegend()
  rp2.SetLeftMargin(0.13);
  rp2.SetRightMargin(0.05);
  rp2.SetUpTopMargin(0.1);
  rp2.SetLowTopMargin(0.02);
  rp2.SetLowBottomMargin(0.35);

  max2 = np.maximum(Mass_projY.GetMaximum()*10,10000)
  Mass_projY.SetMaximum(max2);
  Mass_projY.SetMinimum(0.000001);
#  Mass_projY.GetYaxis().SetRangeUser(0.1,100)
  rp2.GetLowerRefGraph().SetMinimum(-4);
  rp2.GetLowerRefGraph().SetMaximum(4);
  #rp2.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
  #rp2.GetLowerRefGraph().SetLineColor(0) #0
  rp2.GetLowerRefGraph().SetMarkerStyle(20)
  rp2.GetLowerRefGraph().SetMarkerSize(1);
  rp2.GetLowYaxis().SetNdivisions(505);
  rp2.GetLowerRefYaxis().SetTitle("Pull");
  rp2.GetLowerRefYaxis().SetTitleSize(0.05);
  rp2.GetLowerRefYaxis().SetTitleOffset(1);
  rp2.GetLowerRefYaxis().SetLabelSize(0.035);


  rp2.GetLowerRefXaxis().SetTitleSize(0.05);
  rp2.GetLowerRefXaxis().SetTitleOffset(0.8);
  rp2.GetLowerRefXaxis().SetLabelSize(0.035);
  cMass_projY_log.Modified()
  cMass_projY_log.Update()
  #Mass_projY.Draw()
  #Mass_wPred_projY.Draw("SAME")
  #rp2.Draw("X")

  rp2.GetUpperPad().cd();
  legMass.Draw("SAME")
  tex2.Draw("SAME")
  tex3.Draw("SAME")
  tex4.Draw("SAME")
  tex5m.Draw("SAME")
  
  name = newFileDir + "/cMass_log.png"
  cMass_projY_log.SaveAs(name)
  
  cMassOrig_projY = ROOT.TCanvas('cMassOrig_projY', 'cMassOrig_projY',800,800)
  rp0 = ROOT.TRatioPlot(Mass_projY_NotRebinned,Mass_wPred_projY_NotRebinned, "diffsigerrasym")
  rp0.Draw()
  rp0.SetH1DrawOpt("P");
  rp0.SetH2DrawOpt("P");
  rp0.SetLeftMargin(0.13);
  rp0.SetRightMargin(0.05);
  rp0.SetUpTopMargin(0.1);
  rp0.SetLowTopMargin(0.02);
  rp0.SetLowBottomMargin(0.35);
  rp0.GetLowerRefGraph().SetMinimum(-4);
  rp0.GetLowerRefGraph().SetMaximum(4);
  #rp0.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
  #rp0.GetLowerRefGraph().SetLineColor(0) #0
  rp0.GetLowerRefGraph().SetMarkerStyle(20)
  rp0.GetLowerRefGraph().SetMarkerSize(1);
  rp0.GetLowYaxis().SetNdivisions(505);
  rp0.GetLowerRefYaxis().SetTitle("Pull");
  rp0.GetLowerRefYaxis().SetTitleSize(0.05);
  rp0.GetLowerRefYaxis().SetTitleOffset(1);
  rp0.GetLowerRefYaxis().SetLabelSize(0.035);


  rp0.GetLowerRefXaxis().SetTitleSize(0.05);
  rp0.GetLowerRefXaxis().SetTitleOffset(0.8);
  rp0.GetLowerRefXaxis().SetLabelSize(0.035);
  cMassOrig_projY.Modified()
  cMassOrig_projY.Update()

  rp0.GetUpperPad().cd();
  legMass.Draw("SAME")
  tex2.Draw("SAME")
  tex3.Draw("SAME")
  tex4.Draw("SAME")
  tex5m.Draw("SAME")
  cMassOrig_projY.SaveAs(newFileDir + "/cMass_NotRebinned.png")
  
  cMassOrig_projY_log = ROOT.TCanvas('cMassOrig_projY_log', 'cMassOrig_projY_log',800,800)
  cMassOrig_projY_log.SetLogy()
  rp0.Draw()
  legMass.Draw("SAME")
  tex2.Draw("SAME")
  tex3.Draw("SAME")
  tex4.Draw("SAME")
  tex5m.Draw("SAME")
  cMassOrig_projY_log.SaveAs(newFileDir + "/cMass_logy_NotRebinned.png")
  


os.system("cp forWebpage/* "+newFileDir+"/.")
print("scp -r "+ newFileDir + " tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV"+codeVersion+"/.")
