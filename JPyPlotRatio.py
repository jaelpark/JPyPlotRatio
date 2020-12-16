
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as plticker
import matplotlib.container as container

import scipy
from scipy import interpolate

from collections import namedtuple
from ctypes import c_double
import sys #for checking modules

try:
	import ROOT
except ModuleNotFoundError:
	print("pyROOT not found, disabling ROOT integration.");

matplotlib.rcParams["axes.linewidth"] = 1.5;

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

	#a = ROOT.Double(0);
	#b = ROOT.Double(0);
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

def SystematicsPatches(x,y,yerr,s,fc="#FF9848",ec="#CC4F1B",alpha=0.5):
	h = 0.5*s;
	return [patches.Rectangle((x[j]-h,y[j]-0.5*yerr[j]),s,yerr[j],facecolor=fc,edgecolor=ec,alpha=alpha,linewidth=0.5) for j in range(x.size)];

class JPyPlotRatio:
	def __init__(self, panels=(1,1), panelsize=(3,3.375), layoutRatio=0.7, disableRatio=[], rowBounds={}, rowBoundsMax={}, colBounds={}, ratioBounds={}, ratioIndicator=True, ratioType="ratio", ratioSystPlot=False, panelScaling={}, panelPrivateScale=[], panelScaleLoc=(0.92,0.92),panelPrivateRowBounds={}, panelRatioPrivateScale={}, panelRatioPrivateRowBounds={}, systPatchWidth=0.065, panelLabel={}, panelLabelLoc=(0.2,0.92), panelLabelSize=16, panelLabelAlign="right", axisLabelSize=16, tickLabelSize=13, sharedColLabels=False, legendPanel=0, legendLoc=(0.52,0.28), legendSize=10, **kwargs):
		disableRatio = list(set(disableRatio));
		height_ratios = np.delete(np.array(panels[0]*[layoutRatio,1-layoutRatio]),2*np.array(disableRatio,dtype=int)+1);
		self.p,self.ax = plt.subplots(2*panels[0]-len(disableRatio),panels[1]+1,sharex='col',figsize=(panels[1]*panelsize[0],np.sum(height_ratios)*panelsize[1]),gridspec_kw={'width_ratios':[0.0]+panels[1]*[1.0],'height_ratios':height_ratios});
		self.p.subplots_adjust(wspace=0.0,hspace=0.0);

		self.PlotEntry = namedtuple('PlotEntry',['panelIndex','arrays','label','labelLegendId','plotType','kwargs']);

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
		self.rowBounds = rowBounds;
		self.rowBoundsMax = rowBoundsMax;
		self.ratioBounds = ratioBounds;
		self.ratioIndicator = ratioIndicator;
		self.ratioType = ratioType;
		self.ratioSystPlot = ratioSystPlot;
		self.axisLabelSize = axisLabelSize;
		self.tickLabelSize = tickLabelSize;
		self.legendPanel = legendPanel;
		self.legendLoc = legendLoc;
		self.legendSize = legendSize;

		self.ax = np.atleast_2d(self.ax);
		self.s = np.shape(self.ax);
		self.A = np.arange(self.s[0]*self.s[1]).reshape(self.s);
		self.Ay = self.A[:,0]; #y control column
		self.A = np.delete(self.A,0,1); #delete control column

		panelRowsWithRatio = list(set(range(panels[0]))-set(disableRatio));
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

			ratioDefault = {"ratio":"Ratio","diff":"Diff"}[self.ratioType];
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
			#A.xaxis.set_major_locator(plticker.MultipleLocator(1.0));
			#A.xaxis.set_major_locator(plticker.AutoLocator());
			A.xaxis.set_major_locator(plticker.MaxNLocator(6));
			A.xaxis.set_minor_locator(plticker.AutoMinorLocator(5));
			A.yaxis.set_minor_locator(plticker.AutoMinorLocator(5));
			A.xaxis.set_tick_params(labelsize=self.tickLabelSize);
			A.yaxis.set_tick_params(labelsize=self.tickLabelSize);
	
	def Add(self, panelIndex, arrays, label="", labelLegendId=0, plotType="data", **kwargs):
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
			if isinstance(arrays,ROOT.TGraphErrors) or isinstance(arrays,ROOT.TObject):
				x,y,_,yerr = TGraphErrorsToNumpy(arrays);
				scale = kwargs.get("scale",1.0);
				arrays = (x,scale*y,scale*yerr);
			#elif isinstance(arrays,ROOT.TObject):
			#	raise ValueError("Not a valid plot object ROOT.TObject (label: {})".format(label));

		#set uncertainty to zero if not provided
		if len(arrays) < 3:
			arrays = (arrays[0],arrays[1],np.zeros(arrays[1].size));

		try:
			arrays = (arrays[0]+kwargs['xshift'],arrays[1],arrays[2]);
			del kwargs['xshift'];
		except KeyError:
			pass;

		self.plots.append(self.PlotEntry(
			panelIndex=panelIndex,
			arrays=arrays,
			label=label,
			labelLegendId=labelLegendId,
			plotType=plotType,
			kwargs=kwargs));
		self.usedSet.add(panelIndex);

		return len(self.plots)-1; #handle to the plot, given to the Ratio()
	
	#deprecated
	def AddTGraph(self, panelIndex, gr, label="", labelLegendId=0, plotType="data", **kwargs):
		print("WARNING: AddTGraph deprecated: use Add() to plot ROOT.TGraphErrors objects (interchangeable)");
		return self.Add(panelIndex,gr,label,labelLegendId,plotType,**kwargs);
	
	#deprecated
	def AddTH1(self, panelIndex, h1, label="", labelLegendId=0, plotType="histogram", **kwargs):
		print("WARNING: AddTH1 deprecated: use Add() to plot ROOT.TH1 objects (interchangeable)");
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
			print(ysys,yofs);

		elif "ROOT" in sys.modules:
			if isinstance(ysys,ROOT.TGraphErrors):
				_,_,_,ysys = TGraphErrorsToNumpy(ysys);
				yofs = np.zeros(ysys.size);
			elif isinstance(ysys,ROOT.TGraphAsymmErrors):
				_,_,_,_,ye1,ye2 = TGraphAsymmErrorsToNumpy(ysys);
				ysys = 0.5*(ye1+ye2);
				yofs = 0.5*(ye2-ye1);

		if isinstance(ysys,np.ndarray):
			if ysys.size != self.plots[r1].arrays[0].size:
				raise ValueError("Systematics graph number of points does not match with the plot point count");
			ysys *= self.plots[r1].kwargs.get("scale",1.0);

		self.systs.append((r1,ysys,yofs));
	
	def Ratio(self, r1, r2, **kwargs):
		if r1 == r2:
			raise ValueError("Ratio(r1, r2) with r1 == r2");
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

		#plot the data
		for plot in self.plots:
			if plot.plotType == "data":
				if plot.kwargs.get("skipAutolim",False):
					try:
						at = twins[plot.panelIndex];
					except KeyError:
						at = self.ax.flat[a0[plot.panelIndex]].twinx();
						twins[plot.panelIndex] = at;
				else:
					at = self.ax.flat[a0[plot.panelIndex]];
				pr = at.errorbar(*plot.arrays,**{k:plot.kwargs[k] for k in plot.kwargs if k not in ["scale","skipAutolim"]});
				#pr = self.ax.flat[a0[plot.panelIndex]].errorbar(*plot.arrays,**{k:plot.kwargs[k] for k in plot.kwargs if k not in ["scale","skipAutolim"]});
				if plot.label != "":
					labels[plot.label,plot.labelLegendId] = pr;
			elif plot.plotType == "theory":
				p1 = self.ax.flat[a0[plot.panelIndex]].fill_between(plot.arrays[0],plot.arrays[1]-plot.arrays[2],plot.arrays[1]+plot.arrays[2],**{k:plot.kwargs[k] for k in plot.kwargs if k not in ["linecolor","skipAutolim"]});
				pr = (p1,
					self.ax.flat[a0[plot.panelIndex]].plot(*plot.arrays[0:2],color=plot.kwargs.get("linecolor","black"),linestyle=p1.get_linestyle()[0])[0]);
				if plot.label != "":
					labels[plot.label,plot.labelLegendId] = pr;

			elif plot.plotType == "fill_between":
				pr = self.ax.flat[a0[plot.panelIndex]].fill_between(plot.arrays[0],plot.arrays[1],plot.arrays[2],**{k:plot.kwargs[k] for k in plot.kwargs if k not in ["linecolor","skipAutolim"]});
				if plot.label != "":
					labels[plot.label,plot.labelLegendId] = pr;

			elif plot.plotType == "histogram":
				if plot.label != "":
					labels[plot.label,plot.labelLegendId] = pr;
				pr = self.ax.flat[a0[plot.panelIndex]].bar(*plot.arrays[0:2],plot.arrays[0][1]-plot.arrays[0][0],yerr=plot.arrays[2],**plot.kwargs);
				histogramMinY[plot.panelIndex] = np.minimum(plot.arrays[1],histogramMinY[plot.panelIndex]);
				try:
					for plot1 in histograms[plot.panelIndex]:
						mask = plot1.arrays[1] < histogramMinY[plot.panelIndex];
						self.ax.flat[a0[plot1.panelIndex]].bar(plot1.arrays[0][mask],plot1.arrays[1][mask],\
							(plot1.arrays[0][1]-plot1.arrays[0][0]),yerr=plot1.arrays[2][mask],**plot1.kwargs);
					histograms[plot.panelIndex].append(plot);
				except ValueError:
					raise ValueError("Histograms in the same panel must have identical dimensions.");
			elif plot.plotType == "2d":
				#set some defaults
				if "cmap" not in plot.kwargs:
					plot.kwargs["cmap"] = "RdBu_r";
				if "levels" not in plot.kwargs:
					plot.kwargs["levels"] = 10;
				#pr = self.ax.flat[a0[plot.panelIndex]].imshow(plot.arrays,aspect="auto",cmap=plot.kwargs.get("cmap","rainbow"),norm=matplotlib.colors.LogNorm(1,plot.arrays.max()),**plot.kwargs);
				#self.p.colorbar(pr,ax=self.ax.flat[a0[plot.panelIndex]]);
				#pr = self.ax.flat[a0[plot.panelIndex]].contour(*plot.arrays,levels=10,norm=matplotlib.colors.LogNorm(1,plot.arrays[2].max()),colors='k',linewidths=0.2);
				pr = self.ax.flat[a0[plot.panelIndex]].contourf(*plot.arrays,**plot.kwargs);
				self.p.colorbar(pr,ax=self.ax.flat[a0[plot.panelIndex]]);
			else:
				raise ValueError("Invalid plotType specified '{}'. plotType must be 'data', 'theory', 'fill_between', 'histogram' or '2d'.".format(plot.plotType));
			
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
				try:
					bounds = (
						max(self.rowBoundsMax[i][0],bounds[0]),
						min(self.rowBoundsMax[i][1],bounds[1]));
				except KeyError:
					pass;
				self.ax.flat[ry].set_ylim(bounds);

		for i,ry in enumerate(self.A1y):
			try:
				self.ax.flat[ry].set_ylim(self.ratioBounds[i]);
			except KeyError:
				self.ax.flat[ry].set_ylim([0.5,1.5]);

		#plot the ratios
		for robj in self.ratios:
			x1,y1,yerr1 = self.plots[robj[0]].arrays;
			x2,y2,yerr2 = self.plots[robj[1]].arrays;

			plotStyle = robj[2].get("style","default");

			terr1 = np.zeros(yerr1.size);
			terr2 = np.zeros(yerr2.size);

			systs1 = list(filter(lambda t: t[0] == robj[0],self.systs));
			if len(systs1) > 0:
				for sys in systs1:
					serr = (sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*y1);
					terr1 += serr*serr;

			systs2 = list(filter(lambda t: t[0] == robj[1],self.systs));
			if len(systs2) > 0:
				for sys in systs2:
					serr = (sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*y2);
					terr2 += serr*serr;

			sa = max(x1[0],x2[0]);
			sb = min(x1[-1],x2[-1]);
			sx = np.linspace(sa,sb,10*max(x1.size,x2.size));

			if not self.ratioSystPlot and (len(systs1) > 0 or len(systs2) > 0):
				yerr1 = np.sqrt(yerr1*yerr1+terr1);
				yerr2 = np.sqrt(yerr2*yerr2+terr2);
			#else:
			#	terr1 = np.sqrt(terr1);
			#	terr2 = np.sqrt(terr2);
			#	#terr1d = interpolate.interp1d(x1,terr1)(sx);
			#	#terr2d = interpolate.interp1d(x2,terr2)(sx);

			y1d = interpolate.interp1d(x1,y1)(sx);
			yerr1d = interpolate.interp1d(x1,yerr1)(sx);
			y2d = interpolate.interp1d(x2,y2)(sx);
			yerr2d = interpolate.interp1d(x2,yerr2)(sx);

			if self.ratioType == "ratio":
				ratio = y1d/y2d;
				ratio_err = ratio*np.sqrt((yerr2d/y2d)**2+(yerr1d/y1d)**2);

			elif self.ratioType == "diff":
				ratio = y1d-y2d;
				ratio_err = np.sqrt(yerr1d*yerr1d+yerr2d*yerr2d);
			else:
				raise ValueError("Invalid ratioType specified '{}'. ratioType must be either 'ratio' or 'diff'.".format(self.ratioType));

			m = ~np.isnan(ratio);
			sx = sx[m];
			
			ratio = ratio[m];
			ratio_err = ratio_err[m];

			panelIndex = self.plots[robj[0]].panelIndex;

			if not np.ma.is_masked(a1[panelIndex]):
				if self.plots[robj[0]].plotType == "data":
					if plotStyle == "default":
						ratio1d = interpolate.interp1d(sx,ratio,bounds_error=False,fill_value="extrapolate")(x1);
						ratio_err1d = interpolate.interp1d(sx,ratio_err,bounds_error=False,fill_value="extrapolate")(x1);

						self.ax.flat[a1[panelIndex]].errorbar(x1,ratio1d,2*ratio_err1d,**self.plots[robj[0]].kwargs);
					else:
						raise ValueError("Invalid plotStyle specified '{}'. plotStyle must be 'default' when plotType is 'data'.".format(plotStyle));

				elif self.plots[robj[0]].plotType == "theory":
					p1 = self.ax.flat[a1[panelIndex]].fill_between(sx,ratio-ratio_err,ratio+ratio_err,**{k:self.plots[robj[0]].kwargs[k] for k in self.plots[robj[0]].kwargs if k not in ["linecolor"]});

					if plotStyle == "errorbar":
						p1.remove();

						ratio1d = interpolate.interp1d(sx,ratio,bounds_error=False,fill_value="extrapolate")(x1);
						ratio_err1d = interpolate.interp1d(sx,ratio_err,bounds_error=False,fill_value="extrapolate")(x1);

						self.ax.flat[a1[panelIndex]].errorbar(x1,ratio1d,2*ratio_err1d,fmt="s",markerfacecolor=p1.get_facecolor()[0],markeredgecolor="black",linestyle=p1.get_linestyle()[0]);
					elif plotStyle == "default":
						self.ax.flat[a1[panelIndex]].plot(sx,ratio,color=self.plots[robj[0]].kwargs.get("linecolor","black"),linestyle=p1.get_linestyle()[0]);
					else:
						raise ValueError("Invalid plotStyle specified '{}'. plotStyle must be 'default' or 'errorbar' when plotType is 'theory'.".format(plotStyle));

		for sys in self.systs:
			panelIndex = self.plots[sys[0]].panelIndex;
			x1,y1,yerr1 = self.plots[sys[0]].arrays;
			ax = self.ax.flat[a0[panelIndex]];
			xlim = ax.get_xlim();
			patchWidth = self.systPatchWidth*(xlim[1]-xlim[0]);
			syst = (sys[1] if isinstance(sys[1],np.ndarray) else sys[1]*y1);
			for patch in SystematicsPatches(x1,y1+sys[2],2*syst,patchWidth,fc=self.plots[sys[0]].kwargs["color"],ec="black",alpha=0.25):
				ax.add_patch(patch);

			if self.ratioSystPlot and not np.ma.is_masked(a1[panelIndex]):
				syst_y1 = syst/y1;
				terr_max = 1.0+syst_y1;
				terr_min = 1.0-syst_y1;
				p1 = self.ax.flat[a1[panelIndex]].fill_between(x1,terr_min,terr_max,facecolor=self.plots[sys[0]].kwargs["color"],edgecolor="black",alpha=0.25);

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
				if self.ratioType == "ratio":
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

		if isinstance(self.legendPanel,dict):
			for k in self.legendPanel:
				lines = [labels[p] for p in labels if p[1] == k];
				lines = [h[0] if isinstance(h,container.ErrorbarContainer) else h for h in lines];
				labels1 = [p[0] for p in labels if p[1] == k];
				self.ax.flat[a0[self.legendPanel[k]]].legend(lines,labels1,frameon=False,prop={'size':self.legendSize},loc="center",handletextpad=0.25,bbox_to_anchor=self.legendLoc[k]);
		else:
			lines = [labels[p] for p in labels];
			lines = [h[0] if isinstance(h,container.ErrorbarContainer) else h for h in lines];
			labels1 = [p[0] for p in labels];
			self.ax.flat[a0[self.legendPanel]].legend(lines,labels1,frameon=False,prop={'size':self.legendSize},loc="center",handletextpad=0.25,bbox_to_anchor=self.legendLoc);

		self.p.align_labels(self.ax.flat[self.A0y]);
		self.p.align_labels(self.ax.flat[self.A0[:,-1]]);
		self.p.align_labels(self.ax.flat[self.A1y]);
		self.p.align_ylabels(self.ax[:,0]);

	def EnableLatex(self, b):
		matplotlib.rcParams["text.usetex"] = b;
		matplotlib.rcParams['text.latex.preamble'] = [
			r'\usepackage{tgheros}',    # helvetica font
			r'\usepackage{sansmath}',   # math-font matching helvetica
			r'\sansmath'                # actually tell tex to use it!
			r'\usepackage{siunitx}',    # micro symbols
			r'\sisetup{detect-all}',    # force siunitx to use the fonts
			r'\renewcommand{\mathdefault}[1][]{}'
		];
	
	def Save(self, filename):
		self.p.savefig(filename,bbox_inches="tight");

	def Show(self):
		plt.show();

