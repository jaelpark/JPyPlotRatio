
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
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

def SystematicsPatches(x,y,yerr,s,fc="#FF9848",ec="#CC4F1B",alpha=0.5):
	h = 0.5*s;
	return [patches.Rectangle((x[j]-h,y[j]-0.5*yerr[j]),s,yerr[j],facecolor=fc,edgecolor=ec,alpha=alpha,linewidth=0.5) for j in range(len(x))];

class JPyPlotRatio:
	def __init__(self, panels=(1,1), panelsize=(3,3.375), disableRatio=[], rowBounds={}, colBounds={}, ratioBounds={}, ratioIndicator=True, panelScaling={}, panelPrivateScale=[], panelLabel={}, panelLabelLoc=(0.2,0.92), panelLabelSize=16, panelLabelAlign="right", axisLabelSize=16, sharedColLabels=False, legendPanel=0, legendLoc=(0.52,0.28), legendSize=10, **kwargs):
		disableRatio = list(set(disableRatio));
		height_ratios = np.delete(np.array(panels[0]*[0.7,0.3]),2*np.array(disableRatio)+1);
		self.p,self.ax = plt.subplots(2*panels[0]-len(disableRatio),panels[1]+1,sharex='col',figsize=(panels[1]*panelsize[0],np.sum(height_ratios)*panelsize[1]),gridspec_kw={'width_ratios':[0.0]+panels[1]*[1.0],'height_ratios':height_ratios});
		self.p.subplots_adjust(wspace=0.0,hspace=0.0);

		self.plots = [];
		self.systs = [];
		self.ratios = [];
		self.usedSet = set(); #set of plot indices where something has been drawn

		self.panelScaling = panelScaling;
		self.panelPrivateScale = panelPrivateScale;
		self.panelLabel = panelLabel;
		self.panelLabelLoc = panelLabelLoc;
		self.panelLabelSize = panelLabelSize;
		self.panelLabelAlign = panelLabelAlign;
		self.rowBounds = rowBounds;
		self.ratioBounds = ratioBounds;
		self.ratioIndicator = ratioIndicator;
		self.axisLabelSize = axisLabelSize;
		self.legendPanel = legendPanel;
		self.legendLoc = legendLoc;
		self.legendSize = legendSize;

		self.ax = np.atleast_2d(self.ax);
		self.s = np.shape(self.ax);
		self.A = np.arange(self.s[0]*self.s[1]).reshape(self.s);
		self.Ay = self.A[:,0]; #y control column
		self.A = np.delete(self.A,0,1); #delete control column

		panelRowsWithRatio = list(set(range(panels[0]))-set(disableRatio));#[];#[1];#[1,2];
		cr = np.ones(self.s[0]-len(panelRowsWithRatio),dtype=int);
		cr[panelRowsWithRatio] = 2;
		ratioRows = [np.sum(cr[:t+1])-1 for t in panelRowsWithRatio];

		self.A0y = np.delete(self.Ay,ratioRows,0);
		self.A0 = np.delete(self.A,ratioRows,0);
		self.a0 = self.A0.reshape(-1); #plot indices (flat, access with panelIndex)
		noRatioRows = list(set(range(self.s[0]))-set(ratioRows));
		self.A1y = np.delete(self.Ay,noRatioRows,0);
		self.A1 = np.delete(self.A,noRatioRows,0);
		self.a1 = np.ma.array(np.delete(self.A,np.array(ratioRows)-1,0)); #delete all plot rows for which there is a ratio
		self.a1[disableRatio] = np.ma.masked; #mask plot rows for which there is no ratio
		self.a1 = self.a1.reshape(-1);

		#self.A0y = np.delete(self.Ay,2*np.arange(self.s[1])+1,0); #control column for plots
		#self.A0 = np.delete(self.A,2*np.arange(self.s[1])+1,0);
		#self.a0 = self.A0.reshape(-1); #plot indices (flat)
		#self.A1y = np.delete(self.Ay,2*np.arange(self.s[1]),0);
		#self.A1 = np.delete(self.A,2*np.arange(self.s[1]),0);
		#self.a1 = self.A1.reshape(-1); #ratio indices (flat)

		for i,a in enumerate(self.ax[0,1:]):
			try:
				a.set_xlim(colBounds[i]);
			except KeyError:
				pass;

		try:
			if isinstance(kwargs['xlabel'],str):
				if sharedColLabels:
					self.p.text(0.5,0.05,kwargs['xlabel'],size=self.axisLabelSize,horizontalalignment="center");
				else:
					for a in self.ax[-1,1:]:
						a.set_xlabel(kwargs['xlabel'],fontsize=self.axisLabelSize);
			else:
				for i,a in enumerate(self.ax[-1,1:]):
					try:
						a.set_xlabel(kwargs['xlabel'][i],fontsize=self.axisLabelSize);
					except (KeyError,TypeError):
						pass;
		except KeyError:
			pass;

		try:
			if isinstance(kwargs['ylabel'],str):
				for ry in self.A0y:
					self.ax.flat[ry].set_ylabel(kwargs['ylabel'],fontsize=self.axisLabelSize);
			else:
				for i,ry in enumerate(self.A0y):
					try:
						self.ax.flat[ry].set_ylabel(kwargs['ylabel'][i],fontsize=self.axisLabelSize);
					except (KeyError,TypeError):
						pass;

			for ry in self.A1y:
				self.ax.flat[ry].set_ylabel("Ratio",fontsize=self.axisLabelSize);
		except KeyError:
			pass;

		try:
			if isinstance(kwargs['ylabelRight'],str):
				for ry in self.A0[:,-1]:
					self.ax.flat[ry].set_ylabel(kwargs['ylabelRight'],fontsize=self.axisLabelSize);
					self.ax.flat[ry].yaxis.set_label_position("right");
			else:
				for i,ry in enumerate(self.A0[:,-1]):
					try:
						self.ax.flat[ry].set_ylabel(kwargs['ylabelRight'][i],fontsize=self.axisLabelSize);
						self.ax.flat[ry].yaxis.set_label_position("right");
					except (KeyError,TypeError):
						pass;
		except KeyError:
			pass;

		for A in self.ax.flat:
			A.tick_params(which="major",direction="in",length=8.0);
			A.tick_params(which="minor",direction="in",length=2.8);
			#A.xaxis.set_major_locator(plticker.MultipleLocator(1.0));
			#A.xaxis.set_major_locator(plticker.AutoLocator());
			A.xaxis.set_major_locator(plticker.MaxNLocator(6));
			A.xaxis.set_minor_locator(plticker.AutoMinorLocator(5));
			A.yaxis.set_minor_locator(plticker.AutoMinorLocator(5));
			A.xaxis.set_tick_params(labelsize=13);
			A.yaxis.set_tick_params(labelsize=13);
	
	def Add(self, panelIndex, arrays, label="", labelLegendId=0, plotType="data", **kwargs):
		self.plots.append((panelIndex,arrays,label,labelLegendId,plotType,kwargs));
		self.usedSet.add(panelIndex);

		return len(self.plots)-1; #handle to the plot, given to the Ratio()
	
	def AddTGraph(self, panelIndex, gr, label="", labelLegendId=0, plotType="data", scale=1.0, **kwargs):
		#arrays = TGraphErrorsToNumpy(gr);
		x,y,_,yerr = TGraphErrorsToNumpy(gr);
		return self.Add(panelIndex,(x,y*scale,yerr*scale),label,labelLegendId,plotType,**kwargs);
	
	def AddSyst(self, r1, ysys):
		self.systs.append((r1,ysys));
	
	def Ratio(self, r1, r2, **kwargs):
		self.ratios.append((r1,r2,kwargs));
	
	def GetAxes(self, panelIndex):
		return self.ax.flat[self.a0[panelIndex]];
	
	def GetRatioAxes(self, panelIndex):
		return self.ax.flat[self.a1[panelIndex]];
	
	def GetPlot(self):
		return self.p;
	
	def Plot(self):
		#create a matrix of plot indices and remove the control column
		A = self.A;
		A0 = self.A0;
		a0 = self.a0;
		a1 = self.a1;
		#ap = np.arange(s[0]//2*(s[1]-1)).reshape(s[0]//2*(s[1]-1)); #panel indices
		ap = np.arange(self.A.size).reshape(self.A.size);

		labels = {};

		#plot the data
		for plot in self.plots:
			if plot[4] == "data":
				labels[plot[2],plot[3]] = self.ax.flat[a0[plot[0]]].errorbar(*plot[1],**plot[5]);
			elif plot[4] == "theory":
				p1 = self.ax.flat[a0[plot[0]]].fill_between(plot[1][0],plot[1][1]-plot[1][2],plot[1][1]+plot[1][2],**plot[5]);
				labels[plot[2],plot[3]] = (p1,
					#self.ax.flat[a0[plot[0]]].plot(*plot[1][0:2],color=p1.get_edgecolor()[0],linestyle=p1.get_linestyle()[0])[0]);
					self.ax.flat[a0[plot[0]]].plot(*plot[1][0:2],color="black",linestyle=p1.get_linestyle()[0])[0]);
					#self.ax.flat[a0[plot[0]]].plot(*plot[1][0:2],color=matplotlib.colors.colorConverter.to_rgba(p1.get_edgecolor()[0],alpha=1.0),linestyle=p1.get_linestyle()[0])[0]);
			
		for i,(ra0n,ry) in enumerate(zip(A0[:,],self.A0y)):
			try:
				self.ax.flat[ry].set_ylim(self.rowBounds[i]);
			except KeyError:
				bounds = (1e6,-1e6);
				for ra0,rap in zip(ra0n,ap[self.A.shape[1]*i:]):
					if rap not in self.usedSet or rap in self.panelPrivateScale:
						continue;
					ylim0 = self.ax.flat[ra0].get_ylim();
					bounds = (min(bounds[0],ylim0[0]),max(bounds[1],ylim0[1]));
				self.ax.flat[ry].set_ylim(bounds);

		for i,ry in enumerate(self.A1y):
			try:
				self.ax.flat[ry].set_ylim(self.ratioBounds[i]);
			except KeyError:
				self.ax.flat[ry].set_ylim([0.5,1.5]);

		#plot the ratios
		for robj in self.ratios:
			x1,y1,yerr1 = self.plots[robj[0]][1];
			x2,y2,yerr2 = self.plots[robj[1]][1];

			systs1 = list(filter(lambda t: t[0] == robj[0],self.systs));
			if len(systs1) > 0:
				terr = yerr1*yerr1;
				for sys in systs1:
					serr = (sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*yerr1);
					terr += serr*serr;
				yerr1 = np.sqrt(terr);

			systs2 = list(filter(lambda t: t[0] == robj[1],self.systs));
			if len(systs2) > 0:
				terr = yerr2*yerr2;
				for sys in systs2:
					serr = (sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*yerr2);
					terr += serr*serr;
				yerr2 = np.sqrt(terr);

			sa = max(x1[0],x2[0]);
			sb = min(x1[-1],x2[-1]);
			sx = np.linspace(sa,sb,1000);

			y1d = interpolate.interp1d(x1,y1)(sx);
			yerr1d = interpolate.interp1d(x1,yerr1)(sx);
			y2d = interpolate.interp1d(x2,y2)(sx);
			yerr2d = interpolate.interp1d(x2,yerr2)(sx);

			ratio = y1d/y2d;
			ratio_err = ratio*np.sqrt((yerr2d/y2d)**2+(yerr1d/y1d)**2);

			m = ~np.isnan(ratio);
			sx = sx[m];
			
			ratio = ratio[m];
			ratio_err = ratio_err[m];

			panelIndex = self.plots[robj[0]][0];
			if not np.ma.is_masked(a1[panelIndex]):
				if self.plots[robj[0]][4] == "data":
					#if "style" in robj[2] and robj[2]["style"] == "errorbar_fill_syst":
					#	https://alice-publications.web.cern.ch/system/files/draft/5551/2020-03-10-jtpaper_eb.pdf (fig1)
					#	pass;
					ratio1d = interpolate.interp1d(sx,ratio,bounds_error=False,fill_value="extrapolate")(x1);
					ratio_err1d = interpolate.interp1d(sx,ratio_err,bounds_error=False,fill_value="extrapolate")(x1);

					self.ax.flat[a1[panelIndex]].errorbar(x1,ratio1d,2*ratio_err1d,**self.plots[robj[0]][5]);
				elif self.plots[robj[0]][4] == "theory":
					p1 = self.ax.flat[a1[panelIndex]].fill_between(sx,ratio-ratio_err,ratio+ratio_err,**self.plots[robj[0]][5]);
					if robj[2].get("style","default") == "errorbar":
						p1.remove();

						ratio1d = interpolate.interp1d(sx,ratio,bounds_error=False,fill_value="extrapolate")(x1);
						ratio_err1d = interpolate.interp1d(sx,ratio_err,bounds_error=False,fill_value="extrapolate")(x1);

						self.ax.flat[a1[panelIndex]].errorbar(x1,ratio1d,2*ratio_err1d,fmt="s",markerfacecolor=p1.get_facecolor()[0],markeredgecolor=p1.get_edgecolor()[0],color=p1.get_edgecolor()[0],linestyle=p1.get_linestyle()[0]);
					else:
						#self.ax.flat[a1[panelIndex]].plot(sx,ratio,color=p1.get_edgecolor()[0],linestyle=p1.get_linestyle()[0]);
						self.ax.flat[a1[panelIndex]].plot(sx,ratio,color="black",linestyle=p1.get_linestyle()[0]);

		for sys in self.systs:
			x1,y1,yerr1 = self.plots[sys[0]][1];
			ax = self.ax.flat[a0[self.plots[sys[0]][0]]];
			xlim = ax.get_xlim();
			patchWidth = 0.065*(xlim[1]-xlim[0]);
			#for patch in SystematicsPatches(x1,y1,2*(sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*y1),patchWidth,fc="#916f6f",ec="#382a2a"):#fc="#ff5555",ec="black"):
			for patch in SystematicsPatches(x1,y1,2*(sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*y1),patchWidth,fc=self.plots[sys[0]][5]["color"],ec="#382a2a",alpha=0.25):#fc="#ff5555",ec="black"):
				ax.add_patch(patch);

		#adjust ticks
		#TODO: separate loop for ratio
		#for ra0,ra1,rap in zip(a0,a1,ap):
		for ra0,rap in zip(A0.flat,ap):
			try:
				self.ax.flat[ra0].text(*self.panelLabelLoc,self.panelLabel[rap],horizontalalignment=self.panelLabelAlign,verticalalignment="center",transform=self.ax.flat[ra0].transAxes,size=self.panelLabelSize);
			except KeyError:
				pass;

			if rap in self.panelPrivateScale:
				continue;

			ij = np.unravel_index(ra0,self.s);

			ylim1 = self.ax[ij[0],0].get_ylim();

			try:
				(_a1,_b1) = self.panelScaling[rap],0.0;
				self.ax.flat[ra0].set_ylim((ylim1[0]-_b1)/_a1,(ylim1[1]-_b1)/_a1);
				self.ax.flat[ra0].text(0.92,0.92,"$\\mathdefault{{(\\times {a:.1f})}}$"
					.format(a=_a1),horizontalalignment="right",verticalalignment="center",transform=self.ax.flat[ra0].transAxes,size=10);
			except KeyError:
				self.ax.flat[ra0].set_ylim(ylim1);

			ylim0 = self.ax.flat[ra0].get_ylim();

			#make sure the ticks are same as in the first column
			ticks = np.asarray(self.ax[ij[0],0].yaxis.get_majorticklocs());
			ticks = ticks[(ticks >= ylim1[0]) & (ticks <= ylim1[1])];
			ticks_scaled = (ticks-ylim1[0])/(ylim1[1]-ylim1[0]);
			self.ax.flat[ra0].set_yticks((ticks_scaled*(ylim0[1]-ylim0[0]))+ylim0[0]);
			self.ax.flat[ra0].set_yticklabels(ticks);

			plt.setp(self.ax.flat[ra0].get_yticklabels(),visible=False);

		for ra1 in self.A1.flat:
			if self.ratioIndicator:
				xl = self.ax.flat[ra1].get_xlim();
				xs = xl[1]-xl[0];
				self.ax.flat[ra1].plot([xl[0]+0.05*xs,xl[1]-0.05*xs],[1,1],color="gray",linestyle="--",alpha=0.5);

			ij = np.unravel_index(ra1,self.s);

			ylim1 = self.ax[ij[0],0].get_ylim();
			self.ax.flat[ra1].set_ylim(ylim1);

			plt.setp(self.ax.flat[ra1].get_yticklabels(),visible=False);

		#hide ticks from the control plot
		for a in self.ax[:,0]:
			plt.setp(a.get_xticklabels(),visible=False);

		#for ry in self.A0y:
		#	self.ax.flat[ry].yaxis.offsetText.set_visible(True);
		#	offsetStr = self.ax.flat[ry].yaxis.offsetText.get_text();#self.ax.flat[ry].yaxis.get_major_formatter().get_offset();
		#	self.ax.flat[ry].text(0.0,0.9,offsetStr,horizontalalignment="right",verticalalignment="center",transform=self.ax.flat[ry].transAxes,size=13);
		#	#a.text(0.0,0.9,offsetStr,horizontalalignment="right",verticalalignment="center",transform=a.transAxes,size=13);

		if isinstance(self.legendPanel,dict):
			for k in self.legendPanel:
				lines = [labels[p] for p in labels if p[1] == k];
				lines = [h[0] if isinstance(h,container.ErrorbarContainer) else h for h in lines];
				labels1 = [p[0] for p in labels if p[1] == k];
				self.ax.flat[a0[self.legendPanel[k]]].legend(lines,labels1,frameon=False,prop={'size':self.legendSize},loc="center",handletextpad=0.1,bbox_to_anchor=self.legendLoc[k]);
		else:
			lines = [labels[p] for p in labels];
			lines = [h[0] if isinstance(h,container.ErrorbarContainer) else h for h in lines];
			labels1 = [p[0] for p in labels];
			self.ax.flat[a0[self.legendPanel]].legend(lines,labels1,frameon=False,prop={'size':self.legendSize},loc="center",handletextpad=0.1,bbox_to_anchor=self.legendLoc);

		self.p.align_labels(self.ax.flat[self.A0y]);
		self.p.align_labels(self.ax.flat[self.A0[:,-1]]);
		self.p.align_labels(self.ax.flat[self.A1y]);

	def EnableLatex(self, b):
		matplotlib.rcParams["text.usetex"] = b;
		matplotlib.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',    # helvetica font
			r'\usepackage{sansmath}',   # math-font matching helvetica
			r'\sansmath'                # actually tell tex to use it!
			r'\usepackage{siunitx}',    # micro symbols
			r'\sisetup{detect-all}',    # force siunitx to use the fonts
			r'\newcommand{\mathdefault}[1][]{}'
		];
	
	def Save(self, filename):
		self.p.savefig(filename,bbox_inches="tight");

	def Show(self):
		plt.show();

