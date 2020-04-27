import numpy as np
import ROOT
import pickle

import sys
sys.path.append("../"); #add library path

from JPyPlotRatio import JPyPlotRatio

with open("ipg502.pickle","rb") as f:
	ipg502 = pickle.load(f); #load IP-Glasma curves

f = ROOT.TFile("lhc15o.root","read"); #Load measured vn curves

plot = JPyPlotRatio(panels=(3,3), #number of panels (each with ratio)
	xlabel="Centrality (%)", #create a 3x3 grid plot
	panelLabel={i: "$v_{{{}}}$".format(i+2) for i in range(0,8)}, #label the panels v_n
	panelScaling={1:2.0,2:3.0}, #add scaling to some of the panels (panel index : scale factor)
	ylabel="$v_n$",disableRatio=[0]);

for i in range(0,8):
	gr = f.Get("gr_v{}".format(i+2)); #get TGraphErrors object from .root
	hv2 = plot.AddTGraph(i,gr,fmt="s",color="black"); #add it to panel 'i', square markers and black color

	#Add theory curves
	try:
		x,y,yerr = ipg502[i]["data"];
		hv2_th = plot.Add(i,(x,y,yerr),linestyle="-",edgecolor="black",facecolor="orange",alpha=0.8,plotType="theory"); #Add theory curve to panel 'i', lines and orange color

		plot.Ratio(hv2_th,hv2); #Calculate and plot ratio between data and theory
	except KeyError:
		pass;

plot.Plot();
plot.Save("vn.pdf");
plot.Show();

