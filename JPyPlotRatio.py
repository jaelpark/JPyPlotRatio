
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as plticker
import matplotlib.container as container
from matplotlib.path import Path
from matplotlib.transforms import Affine2D

import scipy
from scipy import interpolate

from collections import namedtuple
from ctypes import c_double
import sys #for checking modules

try:
	import ROOT
except ModuleNotFoundError:
	print("pyROOT not found, disabling ROOT support.");

matplotlib.rcParams["axes.linewidth"] = 1.5;

def TGraphToNumpy(gr):
	n = gr.GetN();
	x,y = [np.empty(n) for i in range(2)];
	a = c_double(0);
	b = c_double(0);
	for i in range(0,n):
		gr.GetPoint(i,a,b);
		x[i] = a.value;
		y[i] = b.value;

	return x,y;

def TGraphErrorsToNumpy(gr):
	n = gr.GetN();
	x,y,xerr,yerr = [np.empty(n) for i in range(4)];

	#a = ROOT.Double(0);
	#b = ROOT.Double(0);
	a = c_double(0);
	b = c_double(0);
	for i in range(0,n):
		gr.GetPoint(i,a,b);
		x[i] = a.value;
		y[i] = b.value;
		xerr[i] = gr.GetErrorX(i);
		yerr[i] = gr.GetErrorY(i);

	return x,y,xerr,yerr;

def TGraphAsymmErrorsToNumpy(gr):
	n = gr.GetN();
	x,y,xerr1,xerr2,yerr1,yerr2 = [np.empty(n) for i in range(6)];

	a = c_double(0);
	b = c_double(0);
	for i in range(0,n):
		gr.GetPoint(i,a,b);
		x[i] = a.value;
		y[i] = b.value;
		xerr1[i] = gr.GetErrorXlow(i);
		xerr2[i] = gr.GetErrorXhigh(i);
		yerr1[i] = gr.GetErrorYlow(i);
		yerr2[i] = gr.GetErrorYhigh(i);

	return x,y,xerr1,xerr2,yerr1,yerr2;

def TH2ToNumpy(h):
	nx = h.GetNbinsX();
	ny = h.GetNbinsY();
	z = np.empty((ny,nx));
	x = np.empty(nx);
	y = np.empty(ny);
	for i in range(0,nx):
		for j in range(0,ny):
			z[j,i] = h.GetBinContent(i+1,j+1);
	t = h.GetXaxis();
	for i in range(0,nx):
		x[i] = t.GetBinCenter(i+1);
	t = h.GetYaxis();
	for i in range(0,ny):
		y[i] = t.GetBinCenter(i+1);

	return z,x,y;

def RatioSamples(a1, a2, mode="ratio", freq=10, ratioRange=(-np.inf,np.inf)):
	x1,y1,yerr1 = a1;
	x2,y2,yerr2 = a2;

	#sa = np.max(np.array([x1[0],x2[0]]));
	#sb = np.min(np.array([x1[-1],x2[-1]]));
	sa = np.max(np.array([np.min(x1),np.min(x2)]));
	sb = np.min(np.array([np.max(x1),np.max(x2)]));
	sx = np.linspace(sa,sb,freq*max(x1.size,x2.size));

	y1d = interpolate.interp1d(x1,y1)(sx);
	yerr1d = interpolate.interp1d(x1,yerr1)(sx);
	y2d = interpolate.interp1d(x2,y2)(sx);
	yerr2d = interpolate.interp1d(x2,yerr2)(sx);

	if mode == "ratio":
		ratio = y1d/y2d;
		ratio_err = np.abs(ratio)*np.sqrt((yerr2d/y2d)**2+(yerr1d/y1d)**2);

	elif mode == "diff":
		ratio = y1d-y2d;
		ratio_err = np.sqrt(yerr1d*yerr1d+yerr2d*yerr2d);
	
	elif mode == "sigma":
		sigma = np.sqrt(yerr1d*yerr1d+yerr2d*yerr2d);
		ratio = np.abs((y1d-y2d)/sigma);
		ratio_err = np.zeros(ratio.size);
	
	elif mode == "direct":
		ratio = y1d/y2d;
		ratio_err = yerr1d/y1d;
	
	elif mode == "ratio_error":
		ratio = yerr1d/yerr2d;
		ratio_err = np.zeros(ratio.size);
	
	elif mode == "ratio_rel_error":
		ratio = (yerr1d*y2d)/(yerr2d*y1d);
		ratio_err = np.zeros(ratio.size);
	else:
		raise ValueError("Invalid ratioType specified '{}'. ratioType must be either 'ratio', 'diff', 'sigma', 'direct', 'ratio_error' or 'ratio_rel_error'.".format(mode));

	m = np.bitwise_and(~np.isnan(ratio),ratioRange[0] <= sx,sx < ratioRange[1]);
	sx = sx[m];
	
	ratio = ratio[m];
	ratio_err = ratio_err[m];

	return sx,ratio,ratio_err;

def SystematicsPatches(x, y, yerr, s, fc="#FF9848", ec="#CC4F1B", alpha=0.5,**kwargs):
	h = 0.5*s;
	return [patches.Rectangle((x[j]-h,y[j]-0.5*yerr[j]),s,yerr[j],facecolor=fc,edgecolor=ec,alpha=alpha,linewidth=0.5,**kwargs) for j in range(x.size)];

def StripAttributes(d, fixed, extras=[]):
	return {k:d[k] for k in d if (k not in fixed and k not in ["xshift","scale","skipAutolim","limitMask","noError","style","xerrLinestyle","yerrLinestyle"] and k not in extras)};

class JPyPlotRatio:
	def __init__(self, panels=(1,1), panelsize=(3,3.375), layoutRatio=0.7, disableRatio=[], rowBounds={}, rowBoundsMax={}, colBounds={}, ratioBounds={}, ratioIndicator=True, ratioType="ratio", ratioSystPlot=False, systLegend=True, panelScaling={}, panelPrivateScale=[], panelScaleLoc=(0.92,0.92), panelPrivateRowBounds={}, panelRatioPrivateScale={}, panelRatioPrivateRowBounds={}, systPatchWidth=0.065, panelLabel={}, panelLabelLoc=(0.2,0.92), panelLabelSize=12, panelLabelAlign="right", axisLabelSize=12, tickLabelSize=12, majorTicks=6, majorTickMultiple=None, logScale=False, sharedColLabels=False, hideLegends=False, legendPanel=0, legendLoc=(0.52,0.28), legendLabelSpacing=matplotlib.rcParams['legend.labelspacing'], legendSize=12, sharex='col', **kwargs):
		disableRatio = list(set(disableRatio));
		disableRatio = np.array(disableRatio,dtype=np.int32);
		if np.any(disableRatio >= panels[0]):
			raise ValueError("disableRatio: one or more indices exceeds the number of rows ({})".format(panels[0]));
		height_ratios = np.delete(np.array(panels[0]*[layoutRatio,1-layoutRatio]),2*disableRatio+1);
		self.p,self.ax = plt.subplots(2*panels[0]-disableRatio.size,panels[1]+1,sharex=sharex,figsize=(panels[1]*panelsize[0],np.sum(height_ratios)*panelsize[1]),gridspec_kw={'width_ratios':[0.0]+panels[1]*[1.0],'height_ratios':height_ratios});
		self.p.subplots_adjust(wspace=0.0,hspace=0.0);

		self.PlotEntry = namedtuple('PlotEntry',['panelIndex','arrays','label','labelLegendId','labelOrder','plotType','xshift','limitMask','kwargs']);

		self.plots = [];
		self.systs = [];
		self.ratios = [];
		self.usedSet = set(); #set of plot indices where something has been drawn

		self.panelScaling = panelScaling;
		self.panelPrivateScale = panelPrivateScale;
		self.panelScaleLoc = panelScaleLoc;
		self.panelPrivateRowBounds = panelPrivateRowBounds;
		self.panelRatioPrivateScale = panelRatioPrivateScale;
		self.panelRatioPrivateRowBounds = panelRatioPrivateRowBounds;
		self.systPatchWidth = systPatchWidth;
		self.panelLabel = panelLabel;
		self.panelLabelLoc = panelLabelLoc;
		self.panelLabelSize = panelLabelSize;
		self.panelLabelAlign = panelLabelAlign;
		self.rowBounds = rowBounds; #{0:rowBounds} if type(rowBounds) == tuple else rowBounds;
		self.rowBoundsMax = rowBoundsMax; #{0:rowBoundsMax} if type(rowBoundsMax) == tuple else rowBoundsMax;
		self.ratioBounds = ratioBounds;
		self.ratioIndicator = ratioIndicator;
		self.ratioType = ratioType;
		self.ratioSystPlot = ratioSystPlot;
		self.systLegend = systLegend;
		self.axisLabelSize = axisLabelSize;
		self.tickLabelSize = tickLabelSize;
		self.majorTicks = majorTicks;
		self.majorTickMultiple = majorTickMultiple;
		self.logScale = logScale;
		self.hideLegends = hideLegends;
		self.legendPanel = legendPanel;
		self.legendLoc = legendLoc;
		self.legendLabelSpacing = legendLabelSpacing;
		self.legendSize = legendSize;

		self.ax = np.atleast_2d(self.ax);
		self.s = np.shape(self.ax);
		self.A = np.arange(self.s[0]*self.s[1]).reshape(self.s);
		self.Ay = self.A[:,0]; #y control column
		self.A = np.delete(self.A,0,1); #delete control column

		panelRowsWithRatio = list(set(range(panels[0]))-set(disableRatio.tolist()));
		cr = np.ones(self.s[0]-len(panelRowsWithRatio),dtype=int);
		cr[panelRowsWithRatio] = 2;
		ratioRows = [np.sum(cr[:t+1])-1 for t in panelRowsWithRatio];

		self.A0y = np.delete(self.Ay,ratioRows,0);
		self.A0 = np.delete(self.A,ratioRows,0);
		self.a0 = self.A0.reshape(-1); #plot indices (flat, access with panelIndex)
		noRatioRows = list(set(range(self.s[0]))-set(ratioRows));
		self.A1y = np.delete(self.Ay,noRatioRows,0);
		self.A1 = np.delete(self.A,noRatioRows,0);
		self.a1 = np.ma.array(np.delete(self.A,np.array(ratioRows,dtype=int)-1,0)); #delete all plot rows for which there is a ratio
		self.a1[disableRatio] = np.ma.masked; #mask plot rows for which there is no ratio
		self.a1 = self.a1.reshape(-1);

		#self.A0y = np.delete(self.Ay,2*np.arange(self.s[1])+1,0); #control column for plots
		#self.A0 = np.delete(self.A,2*np.arange(self.s[1])+1,0);
		#self.a0 = self.A0.reshape(-1); #plot indices (flat)
		#self.A1y = np.delete(self.Ay,2*np.arange(self.s[1]),0);
		#self.A1 = np.delete(self.A,2*np.arange(self.s[1]),0);
		#self.a1 = self.A1.reshape(-1); #ratio indices (flat)

		#create a bar-arrow marker
		verts = [(0,0),(0,-5),(-0.3,-4.5),(0.3,-4.5),(0,-5),(0,-5)]; #arrow down, no horizontal bar
		codes = [Path.MOVETO,Path.LINETO,Path.MOVETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY];
		#---
		self.limitMarkerPath = Path(verts,codes);#.transformed(at.transAxes); #arrow down, no horizontal bar
		self.limitMarkerPathFull = Path([(-1,0),(1,0)]+verts,[Path.MOVETO,Path.LINETO]+codes);
		self.limitMarkerPathFullInverse = self.limitMarkerPathFull.transformed(Affine2D().rotate(np.pi));

		vertsDashed = [(0,t) for t in np.linspace(-0,-5,10)]+[(-0.3,-4.5),(0.3,-4.5),(0,-5),(0,-5)]; #arrow down, no horizontal bar
		codesDashed = [[Path.MOVETO,Path.LINETO][i%2] for i in range(0,10)]+[Path.MOVETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY];
		self.limitMarkerPathDashed = Path(vertsDashed,codesDashed);

		verts = [(0,0),(0,-1.5),(-0.3,-1),(0.3,-1),(0,-1.5),(0,0)];
		#self.limitMarkerLegend = Path([(-0.01,0),(0.01,0)]+verts,[Path.MOVETO,Path.LINETO]+codes);
		self.limitMarkerLegend = Path([(-0.5,0),(0.5,0)]+verts,[Path.MOVETO,Path.LINETO]+codes).transformed(Affine2D().translate(0,0.9));

		for i,a in enumerate(self.ax[0,1:]):
			a.xaxis.set_ticks_position('both');
			try:
				if isinstance(colBounds,dict) or isinstance(colBounds,list):
					a.set_xlim(colBounds[i]);
				else:
					a.set_xlim(colBounds); #tuple
			except KeyError:
				pass;

		for i,a in enumerate(self.ax[0:,-1]):
			a.yaxis.set_ticks_position('both');

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

			ratioDefault = {"ratio":"Ratio","diff":"Diff","sigma":"$\\sigma$","direct":"Ratio ($\\sigma_\\mathrm{num}$)","ratio_error":"$\\sigma_{num}/\\sigma_{den}$"}[self.ratioType];
			ratioLabel = kwargs.get('ylabelRatio',ratioDefault);
			if isinstance(ratioLabel,str):
				for ry in self.A1y:
					self.ax.flat[ry].set_ylabel(ratioLabel,fontsize=self.axisLabelSize);
			else:
				for i,ry in enumerate(self.A1y):
					try:
						self.ax.flat[ry].set_ylabel(ratioLabel[i],fontsize=self.axisLabelSize);
					except (KeyError,TypeError):
						self.ax.flat[ry].set_ylabel(ratioDefault,fontsize=self.axisLabelSize);
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
			if self.majorTickMultiple:
				A.xaxis.set_major_locator(plticker.MultipleLocator(10.0));
			else:
				#A.xaxis.set_major_locator(plticker.AutoLocator());
				A.xaxis.set_major_locator(plticker.MaxNLocator(self.majorTicks));
			A.xaxis.set_minor_locator(plticker.AutoMinorLocator(5));
			A.yaxis.set_minor_locator(plticker.AutoMinorLocator(5));
			#FIXME - minor ticks disappear with log scale
			#locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10));
			#A.yaxis.set_minor_locator(locmin);
			#A.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter());
			A.xaxis.set_tick_params(labelsize=self.tickLabelSize);
			A.yaxis.set_tick_params(labelsize=self.tickLabelSize);

		if self.logScale:
			#for A in self.ax.flat:
			for panelIndex in np.concatenate((self.A0y.reshape(-1),self.a0)):
				self.ax.flat[panelIndex].set_yscale("log");

			#locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=999,subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10));
			#locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=999,subs=(0.2,0.4,0.6,0.8));
			locmin = matplotlib.ticker.LogLocator(base=10.0,numticks=999,subs=(10,));
			for ry in self.A0.flat:
				self.ax.flat[ry].yaxis.set_minor_locator(locmin);
				self.ax.flat[ry].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter());
	
	def Add(self, panelIndex, arrays, label="", labelLegendId=0, labelOrder=0, plotType="data", **kwargs):
		if panelIndex >= self.a0.size:
			raise ValueError("panelIndex exceeds the number of panels.");
		if "ROOT" in sys.modules:
			if isinstance(arrays,ROOT.TF1):
				a = c_double(0);
				b = c_double(0);
				arrays.GetRange(a,b);
				vf = np.vectorize(lambda t: arrays.Eval(t));
				x = np.linspace(a,b,1000);
				arrays = x,vf(x);
			if isinstance(arrays,ROOT.TH1):
				arrays = ROOT.TGraphErrors(arrays);
			#if isinstance(arrays,ROOT.TGraph):
			#	x,y = TGraphToNumpy(arrays);
			#	arrays = (x,y);
			if isinstance(arrays,ROOT.TGraphErrors) or isinstance(arrays,ROOT.TObject):
				x,y,_,yerr = TGraphErrorsToNumpy(arrays);
				arrays = (x,y,yerr) if not kwargs.get("noError",False) else (x,y);
			#elif isinstance(arrays,ROOT.TObject):
			#	raise ValueError("Not a valid plot object ROOT.TObject (label: {})".format(label));

		#set uncertainty to zero if not provided
		if len(arrays) < 3 and "2d" not in plotType:
			arrays = (arrays[0],arrays[1],np.zeros(arrays[1].size));

		#if "scale" in kwargs:
		#	scale = kwargs["scale"];
		#	arrays = (arrays[0],scale*arrays[1],scale*arrays[2]);
		for i,a in enumerate(arrays):
			if not isinstance(a,np.ndarray):
				raise ValueError("Array[{}] is not an np.ndarray, or its conversion failed.".format(i));

		try:
			xshift = kwargs['xshift'];
			del kwargs['xshift'];
		except KeyError:
			xshift = 0.0;

		try:
			limitMask = kwargs['limitMask'];
			del kwargs['limitMask'];
		except KeyError:
			limitMask = None; #np.full(arrays[0].size,False);

		self.plots.append(self.PlotEntry(
			panelIndex=panelIndex,
			arrays=arrays,
			label=label,
			labelLegendId=labelLegendId,
			labelOrder=labelOrder,
			plotType=plotType,
			xshift=xshift,
			limitMask=limitMask,
			kwargs=kwargs));
		self.usedSet.add(panelIndex);

		return len(self.plots)-1; #handle to the plot, given to the Ratio()
	
	#deprecated
	def AddTGraph(self, panelIndex, gr, label="", labelLegendId=0, plotType="data", **kwargs):
		print("WARNING: AddTGraph deprecated: use Add() to plot ROOT.TGraphErrors objects (drop-in replacement)");
		return self.Add(panelIndex,gr,label,labelLegendId,plotType,**kwargs);
	
	#deprecated
	def AddTH1(self, panelIndex, h1, label="", labelLegendId=0, plotType="histogram", **kwargs):
		print("WARNING: AddTH1 deprecated: use Add() to plot ROOT.TH1 objects (drop-in replacement)");
		return self.Add(panelIndex,h1,label,labelLegendId,plotType,**kwargs);
	
	def Add2D(self, panelIndex, arrays, **kwargs):
		#if "ROOT" in sys.modules:
		#	if isinstance(arrays,ROOT.TH2):
		#		xe = arrays.GetXaxis();
		#		ye = arrays.GetYaxis();
		#		arrays,_,_ = TH2ToNumpy(arrays);
		#		kwargs["extent"] = kwargs.get("extent",(xe.GetXmin(),xe.GetXmax(),ye.GetXmin(),ye.GetXmax()));
		if "ROOT" in sys.modules:
			if isinstance(arrays,ROOT.TH2):
				#xe = arrays.GetXaxis();
				#ye = arrays.GetYaxis();
				z,x,y = TH2ToNumpy(arrays);
				arrays = (x,y,z);
				#kwargs["extent"] = kwargs.get("extent",(xe.GetXmin(),xe.GetXmax(),ye.GetXmin(),ye.GetXmax()));
		return self.Add(panelIndex,arrays,plotType="2d",**kwargs);
	
	def AddSyst(self, r1, ysys):
		yofs = 0.0;
		if isinstance(ysys,np.ndarray):
			yofs = np.zeros(ysys.size);
		elif isinstance(ysys,tuple):
			ye1,ye2 = ysys[1],ysys[0];
			ysys = 0.5*(ye1+ye2);
			yofs = 0.5*(ye2-ye1);

		elif "ROOT" in sys.modules:
			if isinstance(ysys,ROOT.TGraphErrors):
				_,_,_,ysys = TGraphErrorsToNumpy(ysys);
				yofs = np.zeros(ysys.size);
			elif isinstance(ysys,ROOT.TGraphAsymmErrors):
				_,_,_,_,ye1,ye2 = TGraphAsymmErrorsToNumpy(ysys);
				#TODO: ye1 and ye2 can be both given to errorbar()
				ysys = 0.5*(ye1+ye2);
				yofs = 0.5*(ye2-ye1);
			elif isinstance(ysys,ROOT.TGraph):
				_,ysys = TGraphToNumpy(ysys);

		if isinstance(ysys,np.ndarray):
			if ysys.size != self.plots[r1].arrays[0].size:
				raise ValueError("Systematics graph number of points does not match with the plot point count");
			ysys *= self.plots[r1].kwargs.get("scale",1.0);

		self.systs.append((r1,ysys,yofs));
	
	def Ratio(self, r1, r2, **kwargs):
		#if r1 == r2:
		#	raise ValueError("Ratio(r1, r2) with r1 == r2");
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

		histograms = self.a0.size*[[]];
		histogramMinY = self.a0.size*[np.inf];#np.full(self.a0.size,np.inf);

		twins = {};

		def labelWithScale(label):
			return "{} ($\\times {:.1f}$)".format(label,scale) if np.abs(scale-1) > 1e-4 else label;

		limitToZeroMasks = {};

		#plot the data
		for plotIndex,plot in enumerate(self.plots):
			#apply scaling here so that ratio calculations won't get affected
			scale = plot.kwargs.get('scale',1.0);
			try:
				x,y,yerr = plot.arrays[0],scale*plot.arrays[1],scale*plot.arrays[2];
			except IndexError:
				x = plot.arrays[0]; #gotta be tuple of one
			if plot.plotType == "data":
				xerr = plot.kwargs.get("xerr",None);
				if plot.kwargs.get("skipAutolim",False):
					try:
						at = twins[plot.panelIndex];
					except KeyError:
						at = self.ax.flat[a0[plot.panelIndex]].twinx();
						twins[plot.panelIndex] = at;
				else:
					at = self.ax.flat[a0[plot.panelIndex]];
				if plot.limitMask is not None:
					systs1 = list(filter(lambda sys: sys[0] == plotIndex,self.systs));
					terr = plot.arrays[2]*plot.arrays[2]; #unscaled yerr
					if len(systs1) > 0:
						for sys in systs1:
							serr = (sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*plot.arrays[2]);
							serr[np.isnan(serr)] = 0;
							terr += serr*serr;
					terr = np.sqrt(terr);
					if isinstance(plot.limitMask,str):
						if plot.limitMask == "yerr0":
							ll = y-yerr > 0;
						elif plot.limitMask == "yerrsys0":
							ll = y-terr > 0;
						elif plot.limitMask == "all":
							ll = np.full(y.size,False);
						else:
							raise ValueError("Invalid limitMask string {}.".format(plot.limitMask));
					else:
						ll = ~plot.limitMask;

					zx,zy,zyerr = x[~ll],y[~ll],scale*terr[~ll]; #y already scaled
					limitToZeroMasks[plotIndex] = ll;

					try:
						zxerr = plot.kwargs["xerr"][~ll];
						xerr = plot.kwargs["xerr"][ll];
					except KeyError:
						zxerr = None;

					limitMarker = (self.limitMarkerPathDashed if plot.kwargs.get("linestyle","none") in ["--",":"] else self.limitMarkerPath);
					targs = {"x":zx+plot.xshift,"y":zy+zyerr,"yerr":None,"xerr":zxerr,"zorder":2,"marker":limitMarker if "xerr" in plot.kwargs else self.limitMarkerPathFull,"markersize":50,"fillstyle":"full","linestyle":"none"};
					at.errorbar(**(targs|StripAttributes(plot.kwargs,targs)));
					x,y,yerr = x[ll],y[ll],yerr[ll];

				targs = {"x":x+plot.xshift,"y":y,"yerr":yerr,"xerr":xerr,"zorder":2*len(self.plots)+plotIndex};
				pr = at.errorbar(**(targs|StripAttributes(plot.kwargs,targs)));
				#pr[-1][0].set_linestyle(":");
				if "yerrLinestyle" in plot.kwargs:
					pr[-1][-1].set_linestyle(plot.kwargs["yerrLinestyle"]);
				if plot.label != "":
					if self.systLegend and any((sys[0] == plotIndex for sys in self.systs)):
						labels[labelWithScale(plot.label),plot.labelLegendId,plot.labelOrder] = (patches.Patch(facecolor=plot.kwargs["color"],edgecolor="black",alpha=0.25),pr[0]);
					else:
						labels[labelWithScale(plot.label),plot.labelLegendId,plot.labelOrder] = pr[0];

			elif plot.plotType == "theory":
				if plot.limitMask is not None:
					ll = ~plot.limitMask;
					zx,zy,zyerr = x[~ll],y[~ll],yerr[~ll];
					limitToZeroMasks[plotIndex] = ll;
					targs = {"x":zx+plot.xshift,"y":zy+zyerr,"yerr":None,"xerr":None,"zorder":2,"marker":self.limitMarkerPathFull,"markersize":50,"fillstyle":"full","linestyle":"none"};
					at.errorbar(**(targs|StripAttributes(plot.kwargs,targs,["fmt","mfc"])));
					x,y,yerr = x[ll],y[ll],yerr[ll];
				p1 = self.ax.flat[a0[plot.panelIndex]].fill_between(x+plot.xshift,y-yerr,y+yerr,zorder=plotIndex,**{k:plot.kwargs[k] for k in plot.kwargs if k not in ["linecolor","skipAutolim","noError","scale"]});
				pr = (p1,
					self.ax.flat[a0[plot.panelIndex]].plot(x+plot.xshift,y,color=plot.kwargs.get("linecolor","black"),linestyle=p1.get_linestyle()[0],zorder=plotIndex)[0]);
				if plot.label != "":
					labels[labelWithScale(plot.label),plot.labelLegendId,plot.labelOrder] = pr;

			elif plot.plotType == "fill_between":
				#In this case, y is the lower limit, and yerr the upper.
				pr = self.ax.flat[a0[plot.panelIndex]].fill_between(x+plot.xshift,y,yerr,zorder=plotIndex,**{k:plot.kwargs[k] for k in plot.kwargs if k not in ["linecolor","skipAutolim","noError","scale"]});
				if plot.label != "":
					labels[labelWithScale(plot.label),plot.labelLegendId,plot.labelOrder] = pr;

			elif plot.plotType == "histogram":
				if plot.label != "":
					labels[labelWithScale(plot.label),plot.labelLegendId,plot.labelOrder] = pr;
				pr = self.ax.flat[a0[plot.panelIndex]].bar(x+plot.xshift,y,x[1]-x[0],yerr=yerr,**plot.kwargs);
				#histogramMinY[plot.panelIndex] = np.minimum(plot.arrays[1],histogramMinY[plot.panelIndex]);
				#try:
				#	for plot1 in histograms[plot.panelIndex]:
				#		mask = plot1.arrays[1] < histogramMinY[plot.panelIndex];
				#		self.ax.flat[a0[plot1.panelIndex]].bar(plot1.arrays[0][mask],plot1.arrays[1][mask],\
				#			(plot1.arrays[0][1]-plot1.arrays[0][0]),yerr=plot1.arrays[2][mask],**plot1.kwargs);
				#	histograms[plot.panelIndex].append(plot);
				#except ValueError:
				#	raise ValueError("Histograms in the same panel must have identical dimensions.");
			#elif plot.plotType == "upperLimit":
			#	pr = self.ax.flat[a0[plot.panelIndex]].errorbar(x+plot.xshift,y+yerr,marker=self.limitMarkerPath if "xerr" in plot.kwargs else self.limitMarkerPathFull,markersize=50,fillstyle="full",linestyle="none",**{k:plot.kwargs[k] for k in plot.kwargs if k not in ["scale","skipAutolim","noError","fillstyle","linestyle","markersize","mfc"]});
			#	if plot.label != "":
			#		labels[labelWithScale(plot.label),plot.labelLegendId,plot.labelOrder] = plt.plot([1],linestyle="-",color=plot.kwargs["color"])[0];#pr;

			elif "2d" in plot.plotType:
				if "cmap" not in plot.kwargs:
					plot.kwargs["cmap"] = "RdBu_r";
				if plot.plotType == "2d" or plot.plotType == "2d_contour":
					if "levels" not in plot.kwargs:
						plot.kwargs["levels"] = 10;
					#In this case, yerr is z
					pr = self.ax.flat[a0[plot.panelIndex]].contour(x,y,yerr,levels=10,norm=matplotlib.colors.LogNorm(1,yerr.max()),colors='k',linewidths=0.2);
					pr = self.ax.flat[a0[plot.panelIndex]].contourf(x,y,yerr,**plot.kwargs);
				elif plot.plotType == "2d_histogram":
					if "aspect" not in plot.kwargs:
						plot.kwargs["aspect"] = "auto"; #auto is generally needed, so that the axes won't get messed up
					elif "aspect" == "equal":
						print("WARNING: \"equal\" aspect not supported. Use \"panelsize\" to explicitly control the aspect.");
					pr = self.ax.flat[a0[plot.panelIndex]].imshow(x,**plot.kwargs);
				self.p.colorbar(pr,ax=self.ax.flat[a0[plot.panelIndex]]);

			elif plot.plotType != "hidden":
				raise ValueError("Invalid plotType specified '{}'. plotType must be 'data', 'theory', 'fill_between', 'histogram', '2d' or 'hidden'.".format(plot.plotType));
			
		for i,(ra0n,ry) in enumerate(zip(A0[:,],self.A0y)):
			try:
				if isinstance(self.rowBounds,dict) or isinstance(self.rowBounds,list):
					self.ax.flat[ry].set_ylim(self.rowBounds[i]);
				else:
					self.ax.flat[ry].set_ylim(self.rowBounds);
			except KeyError:
				bounds = (1e6,-1e6);
				for ra0,rap in zip(ra0n,ap[self.A.shape[1]*i:]):
					if rap not in self.usedSet or rap in self.panelPrivateScale:
						continue;
					ylim0 = self.ax.flat[ra0].get_ylim();
					bounds = (min(bounds[0],ylim0[0]),max(bounds[1],ylim0[1]));
				try:
					rowBoundsMax1 = self.rowBoundsMax[i] if isinstance(self.rowBoundsMax,dict) or isinstance(self.rowBoundsMax,list) else self.rowBoundsMax;
					bounds = (
						max(rowBoundsMax1[0],bounds[0]),
						min(rowBoundsMax1[1],bounds[1]));
				except KeyError:
					pass;
				self.ax.flat[ry].set_ylim(bounds);

		for i,ry in enumerate(self.A1y):
			try:
				if isinstance(self.ratioBounds,dict) or isinstance(self.ratioBounds,list):
					self.ax.flat[ry].set_ylim(self.ratioBounds[i]);
				else:
					self.ax.flat[ry].set_ylim(self.ratioBounds);
			except KeyError:
				self.ax.flat[ry].set_ylim([0.5,1.5]);

		#plot the ratios
		for robj in self.ratios:
			x1,y1,yerr1 = self.plots[robj[0]].arrays;
			x2,y2,yerr2 = self.plots[robj[1]].arrays;

			plotStyle = robj[2].get("style","default");
			ratioRange = robj[2].get("xlim",(-np.inf,np.inf));

			terr1 = np.zeros(yerr1.size);
			terr2 = np.zeros(yerr2.size);

			systs1 = list(filter(lambda sys: sys[0] == robj[0],self.systs));
			if len(systs1) > 0:
				for sys in systs1:
					serr = (sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*y1);
					terr1 += serr*serr;

			systs2 = list(filter(lambda sys: sys[0] == robj[1],self.systs));
			if len(systs2) > 0:
				for sys in systs2:
					serr = (sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*y2);
					terr2 += serr*serr;

			panelIndex = self.plots[robj[0]].panelIndex;

			ratioSystPlot = self.ratioSystPlot if isinstance(self.ratioSystPlot,bool) else (panelIndex in self.ratioSystPlot);
			if not ratioSystPlot and (len(systs1) > 0 or len(systs2) > 0):
				yerr1 = np.sqrt(yerr1*yerr1+terr1);
				yerr2 = np.sqrt(yerr2*yerr2+terr2);

			sx,ratio,ratio_err = RatioSamples((x1,y1,yerr1),(x2,y2,yerr2),self.ratioType,ratioRange=ratioRange);
			if robj[2].get("noError",False):
				ratio_err = np.zeros(ratio.size);
			
			xshift = self.plots[robj[0]].xshift;

			if not np.ma.is_masked(a1[panelIndex]):
				if self.plots[robj[0]].plotType in ["data","upperLimit","hidden"]:
					if plotStyle == "default":
						ratio1d = interpolate.interp1d(sx,ratio,bounds_error=False)(x1);
						ratio_err1d = interpolate.interp1d(sx,ratio_err,bounds_error=False)(x1);

						dparams = self.plots[robj[0]].kwargs.copy();
						dparams.update({k:robj[2][k] for k in robj[2]});
						#for k in ["scale","skipAutolim","noError","style","xlim","xshift","xerr"]:
						#	dparams.pop(k,None);
						#self.ax.flat[a1[panelIndex]].errorbar(x1+xshift,ratio1d,ratio_err1d,**dparams);#**self.plots[robj[0]].kwargs);
						targs = {"x":x1+xshift,"y":ratio1d,"yerr":ratio_err1d};
						self.ax.flat[a1[panelIndex]].errorbar(**(targs|StripAttributes(targs,dparams)));
					else:
						raise ValueError("Invalid plotStyle specified '{}'. plotStyle must be 'default' when plotType is 'data'.".format(plotStyle));

				elif self.plots[robj[0]].plotType == "theory":
					dparams = self.plots[robj[0]].kwargs.copy();
					dparams.update({k:robj[2][k] for k in robj[2]});
					targs = {"x":sx+xshift,"y1":ratio-ratio_err,"y2":ratio+ratio_err};
					p1 = self.ax.flat[a1[panelIndex]].fill_between(**(targs|StripAttributes(dparams,targs,["linecolor"])));

					if plotStyle == "errorbar":
						p1.remove();

						ratio1d = interpolate.interp1d(sx,ratio,bounds_error=False,fill_value="extrapolate")(x1);
						ratio_err1d = interpolate.interp1d(sx,ratio_err,bounds_error=False,fill_value="extrapolate")(x1);

						targs = {"x":x1+xshift,"y":ratio1d,"yerr":ratio_err1d,"fmt":"s","markerfacecolor":p1.get_facecolor()[0],"markeredgecolor":"black","linestyle":p1.get_linestyle()[0]};
						self.ax.flat[a1[panelIndex]].errorbar(**targs);
					elif plotStyle == "default":
						dparams = self.plots[robj[0]].kwargs.copy();
						dparams.update({k:robj[2][k] for k in robj[2]});
						if "color" not in dparams and \
							"linecolor" not in dparams:
							dparams['color'] = self.plots[robj[0]].kwargs.get("linecolor","black");
						if "linestyle" not in dparams:
							dparams['linestyle'] = p1.get_linestyle()[0];
						targs = {"x":sx+xshift,"y":ratio};
						self.ax.flat[a1[panelIndex]].plot(**(targs|StripAttributes(targs,dparams,["facecolor","edgecolor"])));

					else:
						raise ValueError("Invalid plotStyle specified '{}'. plotStyle must be 'default' or 'errorbar' when plotType is 'theory'.".format(plotStyle));

		for sys in self.systs:
			panelIndex = self.plots[sys[0]].panelIndex;
			x1,y1,yerr1 = self.plots[sys[0]].arrays;
			xshift = self.plots[sys[0]].xshift;

			scale = self.plots[sys[0]].kwargs.get('scale',1.0);
			y1 *= scale;
			yerr1 *= scale;

			ax = self.ax.flat[a0[panelIndex]];
			xlim = ax.get_xlim();
			patchWidth = self.systPatchWidth*(xlim[1]-xlim[0]);
			syst = (scale*sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*y1);
			for i,patch in enumerate(SystematicsPatches(x1+xshift,y1+sys[2],2*syst,patchWidth,fc=self.plots[sys[0]].kwargs["color"],ec="black",alpha=0.25,zorder=len(self.plots)+sys[0])):
				if limitToZeroMasks.get(sys[0],np.full(yerr1.size,True))[i]:
					ax.add_patch(patch);

			ratioSystPlot = self.ratioSystPlot if isinstance(self.ratioSystPlot,bool) else (panelIndex in self.ratioSystPlot);
			if ratioSystPlot and not np.ma.is_masked(a1[panelIndex]):
				syst_y1 = syst/y1;
				terr_max = 1.0+syst_y1;
				terr_min = 1.0-syst_y1;
				p1 = self.ax.flat[a1[panelIndex]].fill_between(x1+xshift,terr_min,terr_max,facecolor=self.plots[sys[0]].kwargs["color"],edgecolor="black",alpha=0.25);

		#adjust ticks
		for ra0,rap in zip(A0.flat,ap):
			try:
				self.ax.flat[ra0].text(*self.panelLabelLoc,self.panelLabel[rap],horizontalalignment=self.panelLabelAlign,verticalalignment="center",transform=self.ax.flat[ra0].transAxes,size=self.panelLabelSize);
			except KeyError:
				pass;

			if rap in self.panelPrivateScale:
				try:
					self.ax.flat[ra0].set_ylim(self.panelPrivateRowBounds[rap]);
				except KeyError:
					pass;
				continue;

			ij = np.unravel_index(ra0,self.s);

			ylim1 = self.ax[ij[0],0].get_ylim();

			try:
				(_a1,_b1) = self.panelScaling[rap],0.0;
				self.ax.flat[ra0].set_ylim((ylim1[0]-_b1)/_a1,(ylim1[1]-_b1)/_a1);
				self.ax.flat[ra0].text(*self.panelScaleLoc,"$\\mathdefault{{(\\times {a:.1f})}}$"
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

		for ra1,rap in zip(self.A1.flat,ap):
			if self.ratioIndicator:
				xl = self.ax.flat[ra1].get_xlim();
				xs = xl[1]-xl[0];
				if self.ratioType in ["ratio","direct","ratio_error","ratio_rel_error"]:
					ratioLine = np.array([1,1]);
				else:
					ratioLine = np.array([0,0]);
				self.ax.flat[ra1].plot([xl[0]+0.05*xs,xl[1]-0.05*xs],ratioLine,color="gray",linestyle="--",alpha=0.5);

			if rap in self.panelRatioPrivateScale:
				try:
					self.ax.flat[ra1].set_ylim(self.panelRatioPrivateRowBounds[rap]);
				except KeyError:
					pass;
				continue;

			ij = np.unravel_index(ra1,self.s);

			ylim1 = self.ax[ij[0],0].get_ylim();
			self.ax.flat[ra1].set_ylim(ylim1);

			plt.setp(self.ax.flat[ra1].get_yticklabels(),visible=False);

		#hide ticks from the control plot
		for a in self.ax[:,0]:
			plt.setp(a.get_xticklabels(),visible=False);

		for t in twins:
			self.ax.flat[a0[t]].set_zorder(1);
			self.ax.flat[a0[t]].patch.set_visible(False);
			twins[t].axis("off");
			twins[t].set_yticks([]);
			ylim1 = self.ax.flat[a0[t]].get_ylim();
			twins[t].set_ylim(ylim1);

		#for ry in self.A0y:
		#	self.ax.flat[ry].yaxis.offsetText.set_visible(True);
		#	offsetStr = self.ax.flat[ry].yaxis.offsetText.get_text();#self.ax.flat[ry].yaxis.get_major_formatter().get_offset();
		#	self.ax.flat[ry].text(0.0,0.9,offsetStr,horizontalalignment="right",verticalalignment="center",transform=self.ax.flat[ry].transAxes,size=13);
		#	#a.text(0.0,0.9,offsetStr,horizontalalignment="right",verticalalignment="center",transform=a.transAxes,size=13);

		if not self.hideLegends:
			if isinstance(self.legendPanel,dict):
				l1 = [self.legendPanel[k] for k in self.legendPanel];
				for k in self.legendPanel:
					labelsSorted = sorted(list(labels),key=lambda p: p[2]);
					lines = [labels[p] for p in labelsSorted if p[1] == k];
					lines = [h[0] if isinstance(h,container.ErrorbarContainer) else h for h in lines];
					try:
						labels1 = [p[0] for p in labelsSorted if p[1] == k];
						lbs = self.legendLabelSpacing.get(k,matplotlib.rcParams['legend.labelspacing']) \
							if isinstance(self.legendLabelSpacing,dict) else self.legendLabelSpacing;
						l = self.ax.flat[a0[self.legendPanel[k]]].legend(lines,labels1,frameon=False,labelspacing=lbs,prop={'size':self.legendSize},loc="center",handletextpad=0.25,bbox_to_anchor=self.legendLoc[k]);
					except KeyError:
						raise ValueError("Incompatible input legendPanel and legendLoc: the number of entries must be same.");
					try:
						#hack: add_artist must not be called for the last legend for particular panel
						l1.remove(self.legendPanel[k]);
						l1.index(self.legendPanel[k]);
						self.ax.flat[a0[self.legendPanel[k]]].add_artist(l);
					except:
						pass;
			else:
				#one legend only
				#TODO: create multiple legends with labelLegendId, not legendPanel?
				labelsSorted = sorted(list(labels),key=lambda p: p[2]);
				lines = [labels[p] for p in labelsSorted];
				lines = [h[0] if isinstance(h,container.ErrorbarContainer) else h for h in lines];
				labels1 = [p[0] for p in labelsSorted];
				self.ax.flat[a0[self.legendPanel]].legend(lines,labels1,frameon=False,labelspacing=self.legendLabelSpacing,prop={'size':self.legendSize},loc="center",handletextpad=0.25,bbox_to_anchor=self.legendLoc);

		#for A in [self.ax.flat[0]]:#self.ax.flat[1:]:
		#	locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.1,0.2,0.4,0.6,0.8,1,2,4,6,8,10));
		#	A.yaxis.set_minor_locator(locmin);
		#	A.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter());

		self.p.align_labels(self.ax.flat[self.A0y]);
		self.p.align_labels(self.ax.flat[self.A0[:,-1]]);
		self.p.align_labels(self.ax.flat[self.A1y]);
		self.p.align_ylabels(self.ax[:,0]);

	def EnableLatex(self, b, font="latex", **kwargs):
		matplotlib.rcParams["text.usetex"] = b;
		if font == "latex":
			matplotlib.rcParams['mathtext.fontset'] = 'stix';
			matplotlib.rcParams['font.family'] = 'STIXGeneral';
		else:
			matplotlib.rcParams['text.latex.preamble'] = [
				r'\usepackage{tgheros}',    # helvetica font
				r'\usepackage[EULERGREEK]{sansmath}' if kwargs.get('eulergreek',False) else r'\usepackage{sansmath}',   # math-font matching helvetica
				r'\sansmath'                # actually tell tex to use it!
				r'\usepackage{siunitx}',    # micro symbols
				r'\sisetup{detect-all}',    # force siunitx to use the fonts
				#r'\renewcommand{\mathdefault}[1][]{}'
			];
	
	def Save(self, filename):
		self.p.savefig(filename,bbox_inches="tight");

	def Show(self):
		plt.show();

	def Close(self):
		plt.close();

