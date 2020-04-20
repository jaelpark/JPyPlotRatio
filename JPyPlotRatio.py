
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as plticker
import matplotlib.container as container

import scipy
from scipy import interpolate

import ROOT
#matplotlib.font_manager._rebuild()
#matplotlib.rcParams["font.family"] = 'sans-serif';
#matplotlib.rcParams["font.sans-serif"] = 'cm';
#matplotlib.rcParams["font.family"] = 'Helvetica';
#matplotlib.rcParams["font.family"] = 'Droid Sans';
#prop = font_manager.FontProperties(fname='/usr/share/root/fonts/FreeSans.otf');
#print(prop.get_name());
#matplotlib.rcParams["font.family"] = 'Helvetica';
#matplotlib.rcParams["font.cursive"] = 'Helvetica';
#matplotlib.rcParams["font.family"] = 'Helvetica';
#matplotlib.rcParams["font.cursive"] = 'Helvetica'; ###
#matplotlib.rcParams['mathtext.fontset'] = 'custom'; ###
#matplotlib.rcParams["text.latex.preamble"] = r'\newcommand{\mathdefault}[1][]{}';
#matplotlib.rcParams["text.latex.preamble"] = [
#	r'\usepackage{tgheros}',    # helvetica font
#	r'\usepackage{sansmath}',   # math-font matching  helvetica
#	r'\sansmath'                # actually tell tex to use it!
#	r'\usepackage{siunitx}',    # micro symbols
#	r'\sisetup{detect-all}',    # force siunitx to use the fonts
#]
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

def SystematicsPatches(x,y,yerr,s,fc="#FF9848",ec="#CC4F1B"):
	h = 0.5*s;
	return [patches.Rectangle((x[j]-h,y[j]-0.5*yerr[j]),s,yerr[j],facecolor=fc,edgecolor=ec,alpha=0.5,linewidth=0.5) for j in range(len(x))];

class JPyPlotRatio:
	def __init__(self, panels=(1,1), panelsize=(3,3.375), rowBounds={}, colBounds={}, ratioBounds = {}, panelScaling={}, panelLabel={}, panelLabelLoc=(0.2,0.92), panelLabelSize=16, panelLabelAlign="right", axisLabelSize=16, legendLoc=(0.52,0.28), legendSize=10, **kwargs):
		self.p,self.ax = plt.subplots(2*panels[0],panels[1]+1,sharex='col',figsize=(panels[1]*panelsize[0],panels[0]*panelsize[1]),gridspec_kw={'width_ratios':[0.0]+panels[1]*[1.0],'height_ratios':panels[0]*[0.7,0.3]});
		self.p.subplots_adjust(wspace=0.0,hspace=0.0);

		self.plots = [];
		self.systs = [];
		self.ratios = [];
		self.usedSet = set(); #set of plot indices where something has been drawn

		self.panelScaling = panelScaling;
		self.panelLabel = panelLabel;
		self.panelLabelLoc = panelLabelLoc;
		self.panelLabelSize = panelLabelSize;
		self.panelLabelAlign = panelLabelAlign;
		self.rowBounds = rowBounds;
		self.ratioBounds = ratioBounds;
		self.axisLabelSize = axisLabelSize;
		self.legendLoc = legendLoc;
		self.legendSize = legendSize;
		#self.systPatchWidth = systPatchWidth;

		self.s = np.shape(self.ax);
		self.A = np.arange(self.s[0]*self.s[1]).reshape(self.s);
		self.A = np.delete(self.A,0,1); #delete control column
		self.A0 = np.delete(self.A,2*np.arange(self.s[1])+1,0);
		self.a0 = self.A0.reshape(-1); #plot indices

		for i,a in enumerate(self.ax[0,1:]):
			try:
				a.set_xlim(colBounds[i]);
			except KeyError:
				pass;

		try:
			if isinstance(kwargs['xlabel'],str):
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
			for a in self.ax[:,0][::2]:
				a.set_ylabel(kwargs['ylabel'],fontsize=16);
			for a in self.ax[:,0][1::2]:
				a.set_ylabel("Ratio",fontsize=16);
		except KeyError:
			pass;

		for A in self.ax.flat:
			A.tick_params(which="major",direction="in",length=8.0);
			A.tick_params(which="minor",direction="in",length=2.8);
			#A.xaxis.set_major_locator(plticker.MultipleLocator(1.0));
			#A.xaxis.set_minor_locator(plticker.AutoMinorLocator(0.5));
			#A.yaxis.set_minor_locator(plticker.AutoMinorLocator(5));
			A.xaxis.set_tick_params(labelsize=13);
			A.yaxis.set_tick_params(labelsize=13);
	
	def Add(self, panelIndex, arrays, label="", plotType="data", **kwargs):
		self.plots.append((panelIndex,arrays,label,plotType,kwargs));
		self.usedSet.add(panelIndex);

		return len(self.plots)-1; #handle to the plot, given to the Ratio()
	
	def AddTGraph(self, panelIndex, gr, label="", plotType="data", **kwargs):
		#arrays = TGraphErrorsToNumpy(gr);
		x,y,_,yerr = TGraphErrorsToNumpy(gr);
		return self.Add(panelIndex,(x,y,yerr),label,plotType,**kwargs);
	
	def AddSyst(self, r1, ysys):
		self.systs.append((r1,ysys));
	
	def Ratio(self, r1, r2, **kwargs):
		self.ratios.append((r1,r2,kwargs));
	
	def GetAxes(self, panelIndex):
		return self.ax.flat[self.a0[panelIndex]];
	
	def GetPlot(self):
		return self.p;
	
	def Plot(self):
		#create a matrix of plot indices and remove the control column
		s = self.s;
		A = self.A;
		A0 = self.A0;
		a0 = self.a0;
		a1 = np.delete(A,2*np.arange(s[1]),0).reshape(-1); #ratio indices
		ap = np.arange(s[0]//2*(s[1]-1)).reshape(s[0]//2*(s[1]-1));

		labels = {};

		#plot the data
		for plot in self.plots:
			if plot[3] == "data":
				labels[plot[2]] = self.ax.flat[a0[plot[0]]].errorbar(*plot[1],**plot[4]);
			elif plot[3] == "theory":
				p1 = self.ax.flat[a0[plot[0]]].fill_between(plot[1][0],plot[1][1]-plot[1][2],plot[1][1]+plot[1][2],**plot[4]);
				labels[plot[2]] = (p1,
					self.ax.flat[a0[plot[0]]].plot(*plot[1][0:2],color=p1.get_edgecolor()[0],linestyle=p1.get_linestyle()[0])[0]);
			
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

			m = ~np.isnan(ratio);
			sx = sx[m];
			
			ratio = ratio[m];
			ratio_err = ratio_err[m];

			panelIndex = self.plots[robj[0]][0];
			if self.plots[robj[0]][3] == "data":
				ratio1d = interpolate.interp1d(sx,ratio,bounds_error=False,fill_value="extrapolate")(x1);
				ratio_err1d = interpolate.interp1d(sx,ratio_err,bounds_error=False,fill_value="extrapolate")(x1);

				self.ax.flat[a1[panelIndex]].errorbar(x1,ratio1d,ratio_err1d,**self.plots[robj[0]][4]);
			elif self.plots[robj[0]][3] == "theory":
				p1 = self.ax.flat[a1[panelIndex]].fill_between(sx,ratio-ratio_err,ratio+ratio_err,**self.plots[robj[0]][4]);
				if "style" in robj[2] and robj[2]["style"] == "errorbar":
					p1.remove();

					ratio1d = interpolate.interp1d(sx,ratio,bounds_error=False,fill_value="extrapolate")(x1);
					ratio_err1d = interpolate.interp1d(sx,ratio_err,bounds_error=False,fill_value="extrapolate")(x1);

					self.ax.flat[a1[panelIndex]].errorbar(x1,ratio1d,ratio_err1d,fmt="s",markerfacecolor=p1.get_facecolor()[0],markeredgecolor=p1.get_edgecolor()[0],color=p1.get_edgecolor()[0],linestyle=p1.get_linestyle()[0]);
				else:
					self.ax.flat[a1[panelIndex]].plot(sx,ratio,color=p1.get_edgecolor()[0],linestyle=p1.get_linestyle()[0]);

		for sys in self.systs:
			x1,y1,yerr1 = self.plots[sys[0]][1];
			ax = self.ax.flat[a0[self.plots[sys[0]][0]]];
			xlim = ax.get_xlim();
			patchWidth = 0.1*(xlim[1]-xlim[0]);
			for patch in SystematicsPatches(x1,y1,2*(sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*y1),patchWidth,fc="#916f6f",ec="#382a2a"):#fc="#ff5555",ec="black"):
				ax.add_patch(patch);

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
				self.ax.flat[ra0].text(*self.panelLabelLoc,self.panelLabel[rap],horizontalalignment=self.panelLabelAlign,verticalalignment="center",transform=self.ax.flat[ra0].transAxes,size=self.panelLabelSize);
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
		self.ax[0,1].legend(lines,labels,frameon=False,prop={'size':self.legendSize},loc="center",handletextpad=0.1,bbox_to_anchor=self.legendLoc);
	
	def EnableLatex(self, b):
		matplotlib.rcParams["text.usetex"] = b;
		matplotlib.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',    # helvetica font
			r'\usepackage{sansmath}',   # math-font matching  helvetica
			r'\sansmath'                # actually tell tex to use it!
			r'\usepackage{siunitx}',    # micro symbols
			r'\sisetup{detect-all}',    # force siunitx to use the fonts
			r'\newcommand{\mathdefault}[1][]{}'
		];
	
	def Save(self, filename):
		self.p.savefig(filename,bbox_inches="tight");

	def Show(self):
		plt.show();

