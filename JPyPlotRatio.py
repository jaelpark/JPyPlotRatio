
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import ROOT

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

def loopOn(input):
	if isinstance(input,list):
		for i in input:
			yield i;
	else: yield input;

class JPyPlotRatio:
	def __init__(self, panels=(1,1), panelsize=(3,3.375), rowBounds={}, panelScaling={}, **kwargs):
		self.p,self.ax = plt.subplots(2*panels[0],panels[1]+1,sharex=True,figsize=(panels[1]*panelsize[0],panels[0]*panelsize[1]),gridspec_kw={'width_ratios':[0.0]+panels[1]*[1.0],'height_ratios':panels[0]*[0.7,0.3]});
		self.p.subplots_adjust(wspace=0.0,hspace=0.0);

		self.plots = [];
		self.ratios = [];

		self.panelScaling = panelScaling;
		self.rowBounds = rowBounds;

		try:
			self.ax.flat[0].set_xlim(kwargs['xlim']);
		except KeyError:
			pass;

		#try:
		#	#set y-limits for the rows
		#	for a,lim in zip(self.ax[:,0],loopOn(kwargs['ylim'])):
		#		a.set_ylim(lim);
		#except KeyError:
		#	pass;

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
	
	def Add(self, plotIndex, arrays, **kwargs):
		self.plots.append((plotIndex,arrays,kwargs));

		return len(self.plots)-1; #handle to the plot, given to the Ratio()
	
	def AddTGraph(self, plotIndex, gr, **kwargs):
		arrays = TGraphErrorsToNumpy(gr);
		return self.Add(plotIndex,arrays,**kwargs);
	
	def Ratio(self, r1, r2):
		self.ratios.append((r1,r2));
	
	def Plot(self):
		#create a matrix of plot indices and remove the control column
		s = np.shape(self.ax);
		A = np.arange(s[0]*s[1]).reshape(s);
		A = np.delete(A,0,1); #delete control column
		a0 = np.delete(A,2*np.arange(s[1])+1,0).reshape(-1); #plot indices
		a1 = np.delete(A,2*np.arange(s[1]),0).reshape(-1); #ratio indices
		ap = np.arange(s[0]//2*(s[1]-1)).reshape(s[0]//2*(s[1]-1));

		#plot the data
		for plot in self.plots:
			self.ax.flat[a0[plot[0]]].errorbar(*plot[1],**plot[2]);
			
		for i,a in enumerate(self.ax[:,0]):
			try:
				a.set_ylim(rowBounds[i]);
			except:
				bounds = (1e6,-1e6);
				for aa in self.ax[i,1:]:
					ylim0 = aa.get_ylim();
					bounds = (min(bounds[0],ylim0[0]),max(bounds[1],ylim0[1]));
				a.set_ylim(bounds);

		#plot the ratios
		for robj in self.ratios:
			x1,y1,_,yerr1 = self.plots[robj[0]][1];
			x2,y2,_,yerr2 = self.plots[robj[1]][1];

			ratio = y1/y2;
			ratio_err = ratio*np.sqrt((yerr2/y2)**2+(yerr1/y1)**2);

			plotIndex = self.plots[robj[0]][0];
			self.ax.flat[a1[plotIndex]].errorbar(x1,ratio,ratio_err,**self.plots[robj[1]][2]);

		#adjust ticks
		for ra0,ra1,rap in zip(a0,a1,ap):
			ij = np.unravel_index(ra0,s);

			ylim1 = self.ax[ij[0],0].get_ylim();

			try:
				(_a1,_b1) = self.panelScaling[rap],0.0;
				self.ax.flat[ra0].set_ylim((ylim1[0]-_b1)/_a1,(ylim1[1]-_b1)/_a1);
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
			plt.setp(self.ax.flat[ra1].get_yticklabels(),visible=False);

		#hide ticks from the control plot
		for a in self.ax[:,0]:
			plt.setp(a.get_xticklabels(),visible=False);

		xl = self.ax.flat[a1[0]].get_xlim();
		xs = xl[1]-xl[0];
		self.ax.flat[a1[0]].plot([xl[0]+0.05*xs,xl[1]-0.05*xs],[1,1],color="gray",linestyle="--",alpha=0.5);

	def Show(self):
		plt.show();

