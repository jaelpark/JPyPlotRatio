
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import matplotlib.container as container

import scipy
from scipy import interpolate

import ROOT
matplotlib.rcParams["axes.linewidth"] = 1.5;

def TGraphErrorsToNumpy(gr):
	n = gr.GetN();
	x,y,xerr,yerr = np.empty(n),np.empty(n),np.empty(n),np.empty(n);

	a = ROOT.Double(0);
	b = ROOT.Double(0);
	for i in range(0,n):
		gr.GetPoint(i,a,b);
		x[i] = float(a);
		y[i] = float(b);
		xerr[i] = gr.GetErrorX(i);
		yerr[i] = gr.GetErrorY(i);

	return x,y,xerr,yerr;

class JPyPlotRatio:
	def __init__(self, panels=(1,1), panelsize=(3,3.375), rowBounds={}, ratioBounds = {}, panelScaling={}, panelLabel = {}, **kwargs):
		self.p,self.ax = plt.subplots(2*panels[0],panels[1]+1,sharex=True,figsize=(panels[1]*panelsize[0],panels[0]*panelsize[1]),gridspec_kw={'width_ratios':[0.0]+panels[1]*[1.0],'height_ratios':panels[0]*[0.7,0.3]});
		self.p.subplots_adjust(wspace=0.0,hspace=0.0);

		self.plots = [];
		self.ratios = [];
		self.usedSet = set(); #set of plot indices where something has been drawn

		self.panelScaling = panelScaling;
		self.panelLabel = panelLabel;
		self.rowBounds = rowBounds;
		self.ratioBounds = ratioBounds;

		try:
			self.ax.flat[0].set_xlim(kwargs['xlim']);
		except KeyError:
			pass;

		try:
			for a in self.ax[-1,1:]:
				a.set_xlabel(kwargs['xlabel'],fontsize=16);
		except KeyError:
			pass;

		try:
			for a in self.ax[:,0][::2]:
				a.set_ylabel(kwargs['ylabel'],fontsize=16);
			for a in self.ax[:,0][1::2]:
				a.set_ylabel("Ratio",fontsize=16);
		except KeyError:
			pass;

		for A in self.ax.flat:
			A.tick_params(which="major",direction="in",length=8.0);
			A.tick_params(which="minor",direction="in",length=2.8);
			#A.xaxis.set_major_locator(plticker.MultipleLocator(10.0));
			#A.xaxis.set_minor_locator(plticker.AutoMinorLocator(5));
			#A.yaxis.set_minor_locator(plticker.AutoMinorLocator(5));
			A.xaxis.set_tick_params(labelsize=13);
			A.yaxis.set_tick_params(labelsize=13);
	
	def Add(self, plotIndex, arrays, label="", **kwargs):
		self.plots.append((plotIndex,arrays,label,kwargs));
		self.usedSet.add(plotIndex);

		return len(self.plots)-1; #handle to the plot, given to the Ratio()
	
	def AddTGraph(self, plotIndex, gr, label="", **kwargs):
		#arrays = TGraphErrorsToNumpy(gr);
		x,y,_,yerr = TGraphErrorsToNumpy(gr);
		return self.Add(plotIndex,(x,y,yerr),label,**kwargs);
	
	def Ratio(self, r1, r2):
		self.ratios.append((r1,r2));
	
	def Plot(self):
		#create a matrix of plot indices and remove the control column
		s = np.shape(self.ax);
		A = np.arange(s[0]*s[1]).reshape(s);
		A = np.delete(A,0,1); #delete control column
		A0 = np.delete(A,2*np.arange(s[1])+1,0);
		a0 = A0.reshape(-1); #plot indices
		a1 = np.delete(A,2*np.arange(s[1]),0).reshape(-1); #ratio indices
		ap = np.arange(s[0]//2*(s[1]-1)).reshape(s[0]//2*(s[1]-1));

		labels = {};

		#plot the data
		for plot in self.plots:
			labels[plot[2]] = self.ax.flat[a0[plot[0]]].errorbar(*plot[1],**plot[3]);
			
		for i,ra0n in enumerate(A0[:,]):
			try:
				self.ax[2*i,0].set_ylim(self.rowBounds[i]);
			except KeyError:
				bounds = (1e6,-1e6);
				for ra0,rap in zip(ra0n,ap[(s[1]-1)*i:]):
					if rap not in self.usedSet:
						continue;
					ylim0 = self.ax.flat[ra0].get_ylim();
					bounds = (min(bounds[0],ylim0[0]),max(bounds[1],ylim0[1]));
				self.ax[2*i,0].set_ylim(bounds);

			try:
				self.ax[2*i+1,0].set_ylim(self.ratioBounds[i]);
			except KeyError:
				self.ax[2*i+1,0].set_ylim([0.5,1.5]);

		#plot the ratios
		for robj in self.ratios:
			x1,y1,yerr1 = self.plots[robj[0]][1];
			x2,y2,yerr2 = self.plots[robj[1]][1];

			sa = max(x1[0],x2[0]);
			sb = min(x1[-1],x2[-1]);
			sx = np.linspace(sa,sb,1000);

			y1d = interpolate.interp1d(x1,y1,bounds_error=False,fill_value="extrapolate")(sx);
			yerr1d = interpolate.interp1d(x1,yerr1,bounds_error=False,fill_value="extrapolate")(sx);
			y2d = interpolate.interp1d(x2,y2,bounds_error=False,fill_value="extrapolate")(sx);
			yerr2d = interpolate.interp1d(x2,yerr2,bounds_error=False,fill_value="extrapolate")(sx);

			ratio = y1d/y2d;
			ratio_err = ratio*np.sqrt((yerr2d/y2d)**2+(yerr1d/y1d)**2);

			plotIndex = self.plots[robj[0]][0];
			self.ax.flat[a1[plotIndex]].fill_between(sx,ratio-ratio_err,ratio+ratio_err,color=self.plots[robj[1]][3]["color"],alpha=0.5);
			self.ax.flat[a1[plotIndex]].plot(sx,ratio,linestyle="-",color=self.plots[robj[1]][3]["color"]);

		#adjust ticks
		for ra0,ra1,rap in zip(a0,a1,ap):
			ij = np.unravel_index(ra0,s);

			ylim1 = self.ax[ij[0],0].get_ylim();

			try:
				(_a1,_b1) = self.panelScaling[rap],0.0;
				self.ax.flat[ra0].set_ylim((ylim1[0]-_b1)/_a1,(ylim1[1]-_b1)/_a1);
				self.ax.flat[ra0].text(0.92,0.92,"$\\mathdefault{{(\\times {a:.1f})}}$"
					.format(a=_a1),horizontalalignment="right",verticalalignment="center",transform=self.ax.flat[ra0].transAxes,size=10);
			except KeyError:
				self.ax.flat[ra0].set_ylim(ylim1);

			try:
				self.ax.flat[ra0].text(0.2,0.92,self.panelLabel[rap],horizontalalignment="right",verticalalignment="center",transform=self.ax.flat[ra0].transAxes,size=16);
			except KeyError:
				pass;

			ylim0 = self.ax.flat[ra0].get_ylim();

			#make sure the ticks are same as in the first column
			ticks = np.asarray(self.ax[ij[0],0].yaxis.get_majorticklocs());
			ticks = ticks[(ticks >= ylim1[0]) & (ticks <= ylim1[1])];
			ticks_scaled = (ticks-ylim1[0])/(ylim1[1]-ylim1[0]);
			self.ax.flat[ra0].set_yticks((ticks_scaled*(ylim0[1]-ylim0[0]))+ylim0[0]);
			self.ax.flat[ra0].set_yticklabels(ticks);

			plt.setp(self.ax.flat[ra0].get_yticklabels(),visible=False);
			plt.setp(self.ax.flat[ra1].get_yticklabels(),visible=False);

			xl = self.ax.flat[ra1].get_xlim();
			xs = xl[1]-xl[0];
			self.ax.flat[ra1].plot([xl[0]+0.05*xs,xl[1]-0.05*xs],[1,1],color="gray",linestyle="--",alpha=0.5);

			ij = np.unravel_index(ra1,s);

			ylim1 = self.ax[ij[0],0].get_ylim();
			self.ax.flat[ra1].set_ylim(ylim1);

		#hide ticks from the control plot
		for a in self.ax[:,0]:
			plt.setp(a.get_xticklabels(),visible=False);

		lines = [labels[p] for p in labels];
		lines = [h[0] if isinstance(h,container.ErrorbarContainer) else h for h in lines];
		self.ax[0,1].legend(lines,labels,frameon=False,prop={'size':10},loc="center",handletextpad=0.1,bbox_to_anchor=(0.52,0.28));
	
	def Save(self, filename):
		self.p.savefig(filename,bbox_inches="tight");

	def Show(self):
		plt.show();

