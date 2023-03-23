# JPyPlotRatio
Multipanel plotting class with ratio panels.

## Required packages
- python3 (https://www.python.org)
- ROOT (http://root.cern.ch)
	- with cmake -Dpython3="ON"
	- OR in MAC
		- base=/usr/local/opt/python3/Frameworks/Python.framework/Versions/3.7
		- cmake ../root-6.20.00 -DCMAKE_INSTALL_PREFIX={your build dir} -DPYTHON_EXECUTABLE=${base}/bin/python3 -DPYTHON_INCLUDE_DIR=${base}/Headers -DPYTHON_LIBRARY=${base}/lib/libpython3.5m.dylib -Dgnuinstall=ON -Dpython3=ON -Droofit=ON -Dminuit2=ON
	- cmake -DCMAKE_INSTALL_PREFIX=../root_install ../root_src -DPYTHON_EXECUTABLE=/usr/bin/python3  -Dgnuinstall=ON -Dpython3=ON -Droofit=ON -Dminuit2=ON
	- Recommended via hombrew, brew install root --build-from-source
- https://matplotlib.org
	- pip3 install matplotlib scipy
- For latex style texting via `plot.EnableLatex(True);`
	- http://www.tug.org/mactex/
	- yum install texlive-*

## Usage
Set the PYTHONPATH with the github you checkout like,
PYTHONPATH=: yourdir/GitHub/JPyPlotRatio
Import `JPyPlotRatio` and create JPyPlotRatio instance like
import sys
sys.path.append("JPyPlotRatio");
import JPyPlotRatio


```python
plot = JPyPlotRatio(panels=(1,2),
	... options
```

Returns a JPyPlotRatio class instance.

## Ratio Usages
By default ratioSystPlot is False.
Ratio() combines stat and syst, unless ratioSystPlot is enabled(`ratioSystPlot=True`).
Make sure that `plot.AddSyst(dataPlotId,syst)` is in the code so that Ratio function knows about it.  

In case of `plot.Ratio(model,data)` model/data 
The ratio is based on the statical errors only both from the data and model (if `ratioSystPlot=True`), the systematic errors from the data (relative errors) will be drawn in bands around 1.

## Parameter Description

Parameter | Type | Description
--- | --- | ---
panels | Tuple (h, w) | Multipanel plot dimensions. Default (1,1) (one panel)
panelsize | Tuple (w, h) | Dimensions of one panel
layoutRatio | float | Relative proportions between the plot panel and ratio panel. Default 0.7
disableRatio | list `[rowN1,rowN2,...]` | Row indices for which ratio plot won't be shown 
rowBounds | Dict `{rowN:(ymin,ymax), ...}` or tuple | Dictionary of y-limits for each row. If tuple, apply common to each row.
rowBoundsMax | Dict `{rowN:(ymin,ymax), ...}` or tuple | Dictionary of minimum and maximum y-limits for each row - applied if plot point exceeds these limits. Ignored for a row if `rowBounds` has a setting for it. If tuple, apply common to each row.
colBounds | Dict `{colN:(xmin,xmax), ...}` | Dictionary of x-limits for each column
ratioBounds | Dict `{rowN:(ymin,ymax), ...}` | Dictionary of y-limits for the ratio panels in each row
ratioIndicator | bool | If True, a dashed line will be drawn to indicate a zero difference line
ratioType | str | Ratio panel approach: "ratio": (a/b) with error propagation, "diff": (a-b) with error propagation, "sigma": |a-b|/sigma, "direct": (a/b) with error (Da/Db), "ratio_err": (Da/Db), "ratio_rel_error": (Da*b)/(Db*a)
ratioSystPlot | bool | Plot systematics separately in ratio plot
systLegend | bool | Draw systematic error patches in legend. Default true
panelScaling | Dict `{panelIndex:scale, ...}` | Dictionary of scale factors for the plots in each panel
panelPrivateScale | list `[panelIdex1,panelIndex2,...]` | List of panels that should have their own y-axis scale instead of a shared one
panelPrivateRowBounds | Dict `{panelIndex:(ymin,ymax),...}` | Dictionary of y-limits for each panel included in `panelPrivateScale`
panelRatioPrivateScale | list `[panelIndex1,panelIndex2,...]` | List of panels ratios that should have their own y-axis instead of a shared one
panelRatioPrivateRowBounds | Dict `{panelIndex:(ymin,ymax),...}` | Dictionary of y-limits for each panel ratio included in `panelPrivateScale`
systPatchWidth | float | Fractional width of the systematic uncertainty patches with relation to the panel width. Default 0.065.
xlabel | str or Dict `{colId:str, ...}` | xlabel for all panels (str), or dictionary of xlabels for each column
ylabel | str or Dict `{rowN:str, ...}` | ylabel for all panels
ylabelRight | str or Dict `{rowN:str, ...}` | ylabel for all panels on the right side
ylabelRatio | str or Dict `{rowN:str, ...}` | ylabel for the ratio panels (default "Ratio" or "Diff")
axisLabelSize | int | Axis label text size. Default 16
tickLabelSize | int | Tick labe size. Default 13
majorTicks | int | Maximum number of major ticks on the x-axis. Default 6
majorTickMultiple | int | Multiples of the major ticks. By default not used (None).
logScale | bool | Apply logarithmic scale to each panel.
panelLabel | Dict `{panelIndex:label(str), ...}` | Dictionary of panel labels
panelLabelLoc | Tuple (x, y) | Location for the panel labels in each panel. Default `(0.2,0.92)`
panelLabelSize | int | Panel label text size. Default 16
panelLabelAlign | str | Text alignment for the panel labels, "left", "center", "right" (default).
legendPanel | int or Dict `{legendId:panelIndex, ...}` | Index (indices) of the panel(s) where to pace the legend. Use dictionary if multiple legends is wanted.
legendLoc | Tuple (x, y) or Dict `{legendId:(x,y), ...}` | Legend location(s) in the panel(s)
legendLabelSpacing | float or Dict `{legendId:(x,y), ...}` | Legend label spacing
legendSize | int | Legend text size. Default 10
sharex | str | sharex argument passed on to subplots() Default "col"

Plot curves with `plot.Add(...)`. Input curve should be either as numpy arrays, or ROOT objects. `Add` automatically converts a TGraphErrors and TH1 object to numpy arrays and plots them.

```python
plotIndex = plot.Add(panelIndex, arrays=(x,y,yerr) or gr=ROOT.TObject, label="", labelLegendId=0, plotType="data", **plotParams)
```

Returns the index for the newly added plot (int), which can be used as a reference for ratio plotting.

Parameter | Type | Description
--- | --- | ---
panelIndex | int | Index of the panel to plot into. Panels are indexed from 0 to _n_ in a row-wise order.
arrays | Tuple (x, y, yerr), `np.array` | Tuple of numpy arrays: x-values, y-values, and y sigma values
gr | `ROOT.TObject` | ROOT object to be plotted. Supported objects are `TGraphErrors`, `TH1`, and `TF1`.
label | str | Label to use for the legend. Use same label for plots for which same legend entry is to be used.
labelLegendId | int | Identifier of the legend to which the plot will be labeled
labelOrder | int | Explicit ordering of legend labels
plotType | str | Plotting method and type. See table below for the available options. Default "data"
scale | float | Scale factor for y values of the TGraphErrors object
**plotParams | dict | Supplementary parameters passed to matplotlib `errorbar, `fill_between` or `bar` plotting methods depending on plotType.
noError | bool | If true, the error from ROOT objects is not plotted, without having to explicitly remove it first. Default False
limitMask | `np.array` of bool, `str` | Draw upper limit at `y + yerr` for data points for which the mask element is true. Default None (no upper limits). Value can also be "yerr0" (upper limit when uncertainty below 0), "yerrsys0" (same, but with also systematics combined), or "all" (all points as upper limit)

The following options are available for the `plotType` parameter. Additional styling may be specified by supplying valid arguments that will be passed on to the respective matplotlib methods.

Value | Description
--- | ---
data | Default: draw the graph as points and errorbars using the `errorbar` method.
theory | Draw a theory curve with a colorband (`fill_between` method) using typical HEP plot styling.
fill_between | Draw using `fill_between` method, with user supplying the lower and higher y-limits instead of y and y uncertainty.
histogram | Draw a histogram using `bar` method.
2d_contour | Draw a 2D color and contour map with with `contour` and `contourf` methods.
2d_histogram | Draw a 2D histogram with a colorbar with `imshow` method.

A corresponding method is used to plot 2D graphs and histograms:

```python
plotIndex = plot.Add2D(panelIndex, arrays=np.ndarray or gr=ROOT.TH2, **plotParams)
```

Parameter | Type | Description
--- | --- | ---
panelIndex | int | Index of the panel to plot into. Panels are indexed from 0 to _n_ in a row-wise order.
arrays | `np.ndarray`, `ROOT.TH2` | numpy 2d matrix to be plotted, or `ROOT.TH2` histogram.
**plotParams | dict | Supplementary parameters passed to matplotlib `imshow`.

Ratio plots can be added by calling `plot.Ratio(...)`:

```python
plot.Ratio(r1, r2, style, xlim, noError, **plotParams)
```

Parameter | Type | Description
--- | --- | ---
r1 | int | Index of the numerator plot
r2 | int | Index of the denominator plot. The plot style for the ratio will be inherited from this curve.
style | str | Style for the ratio plot. Specify "errorbar" to always have point plot ratio curves
xlim | tuple | Set the x-limits of the ratio graph. Default `(-np.inf,np.inf)`
noError | bool | Disable error bars for the ratio if true. Default false
**plotParams | dict | By default, the ratio plot will inherit its style properties from the main plot. However, those properties can be overriden by supplying them separately to `Ratio()`.

Add systematic uncertainties with `plot.AddSyst(...)`:

```python
plot.AddSyst(r1, ysys)
```

By default, the systematic error is included in the ratio plot as a quadratic sum. To plot the ratio uncertainty for the systematics separately, specify `ratioSystPlot=True` while constructing JPyPlotRatio. The systematic uncertainty should be a scalar (to be multiplied with the y-data points for the error), numpy array (uncertainty for each y-data point), TGraphErrors (error to be automatically converted to an array), TGraphAsymmErrors (asymmetric error to be converted to an array and displaced).

Parameter | Type | Description
--- | --- | ---
r1 | int | Index of the plot to add systematic uncertainty patches
ysys | float, `np.array`, `ROOT.TObject` | (Relative) systematic uncertainty. If float, the value will be multiplied by the y-values of the curve, else if numpy array, the values from this array will be used directly. ROOT objects will be converted to a numpy array before plotting.

Finally, create the plot by calling `plot.Plot()`. A plot can be saved with `plot.Save(filename)`, and shown in a window wit `plot.Show()`.

## Additional functions

```python
x,y,xerr,yerr = TGraphErrorsToNumpy(gr)
```

Convert `TGraphErrors` to numpy arrays.

```python
x,y,xerr1,xerr2,yerr1,yerr2 = TGraphAsymmErrorsToNumpy(gr)
```

Convert `TGraphAsymmErrors` to numpy arrays.

```python
z,x,y = TH2ToNumpy(h)
```

Convert `TH2` to numpy arrays: 2d-array `z`, and bin arrays `x` and `y`.

```python
x,ratio,ratio_err = RatioSamples((x1,y1,yerr1), (x2,y2,yerr2), mode="ratio", freq=10)
```

Calculate ratio or difference between two (unmatching) graphs of data. The sample frequency (default 10) determines the number of sample points multiplied by the size of the largest array of data.

