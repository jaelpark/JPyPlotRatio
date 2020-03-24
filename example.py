import numpy as np
import ROOT
from JPyPlotRatio import JPyPlotRatio

f = ROOT.TFile("lhc15o.root","read");
v2 = f.Get("gr_v2");
v3 = f.Get("gr_v3");
v4 = f.Get("gr_v4");

#plot = JPyPlotRatio(panels=(1,2),ylim=(0.0,0.15),panelScaling={1:2.0},xlabel="Centrality (%)",ylabel="$v_n$");
#plot = JPyPlotRatio(panels=(1,2),ylim=(0.0,0.15),xlabel="Centrality (%)",ylabel="$v_n$");
#plot = JPyPlotRatio(panels=(1,2),xlabel="Centrality (%)",ylabel="$v_n$");
hv2 = plot.AddTGraph(0,v2,color="black");
hv3 = plot.AddTGraph(0,v3,color="red");
hv4 = plot.AddTGraph(1,v4,color="green");
plot.Ratio(hv2,hv3);

plot.Plot();
plot.Show();

