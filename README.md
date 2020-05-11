# JPyPlotRatio
Multipanel plotting class with ratio panels.

## Required packages
- python3 (https://www.python.org)
- ROOT (http://root.cern.ch)
	- with cmake  -Dpython3="ON"
	- OR in MAC
		- base=/usr/local/opt/python3/Frameworks/Python.framework/Versions/3.7
		- cmake ../root-6.20.00 -DCMAKE_INSTALL_PREFIX={your build dir} -DPYTHON_EXECUTABLE=${base}/bin/python3 -DPYTHON_INCLUDE_DIR=${base}/Headers -DPYTHON_LIBRARY=${base}/lib/libpython3.5m.dylib -Dgnuinstall=ON -Dpython3=ON -Droofit=ON -Dminuit2=ON
- https://matplotlib.org
	- pip3 install matplotlib
- For latex style texting via plot.EnableLatex(True);
	- http://www.tug.org/mactex/
	- yum install texlive-*

## Usage
Import `JPyPlotRatio` and create JPyPlotRatio instance:

```python
plot = JPyPlotRatio(panels=(1,2),
	... options
```

Returns a JPyPlotRatio class instance.

Parameter | Type | Description
--- | --- | ---
panels | Tuple (h, w) | Multipanel plot dimensions. Default (1,1) (one panel)
panelsize | Tuple (w, h) | Dimensions of one panel
disableRatio | list `[rowN1,rowN2,...]` | Row indices for which ratio plot won't be shown 
rowBounds | Dict `{rowN:(ymin,ymax), ...}` | Dictionary of y-limits for each row
colBounds | Dict `{colN:(xmin,xmax), ...}` | Dictionary of x-limits for each column
ratioBounds | Dict `{rowN:(ymin,ymax), ...}` | Dictionary of y-limits for the ratio panels in each row
ratioIndicator | bool | If True, a dashed line will be drawn to indicate a zero difference line
ratioType | str | Difference type "ratio" (a/b) or "diff" (a-b)
ratioSystPlot | bool | Plot systematics separately in ratio plot
panelScaling | Dict `{panelIndex:scale, ...}` | Dictionary of scale factors for the plots in each panel
panelPrivateScale | list `[panelIdex1,panelIndex2,...] | List of panels that should have their own y-axis scale instead of shared one
panelPrivateRowBounds | Dict `{panelIndex:(ymin,ymax),...}` | Dictionary of y-limits for each panel included in `panelPrivateScale`
xlabel | str or Dict `{colId:str, ...}` | xlabel for all panels (str), or dictionary of xlabels for each column
ylabel | str | ylabel for all panels
axisLabelSize | int | Axis label text size. Default 16
panelLabel | Dict `{panelIndex:label(str), ...}` | Dictionary of panel labels
panelLabelLoc | Tuple (x, y) | Location for the panel labels in each panel. Default `(0.2,0.92)`
panelLabelSize | int | Panel label text size. Default 16
panelLabelAlign | str | Text alignment for the panel labels, "left", "center", "right" (default).
legendPanel | int or Dict `{legendId:panelIndex, ...}` | Index (indices) of the panel(s) where to pace the legend
legendLoc | Tuple (x, y) or Dict `{legendId:(x,y), ...}` | Legend location(s) in the panel(s)
legendSize | int | Legend text size. Default 10

Plot curves with `plot.Add(...)` or `plot.AddTGraph(...)`. `AddTGraph` automatically converts a TGraphErrors object to numpy arrays and plots them.

```python
plotIndex = plot.Add(panelIndex, (x, y, yerr), label, labelLegendId, plotType, **plotParams)
plotIndex = plot.AddTGraph(panelIndex, gr, label, labelLegendId, plotType, scale, **plotParams)
```

Returns the index for the newly added plot (int).

Parameter | Type | Description
--- | --- | ---
panelIndex | int | Index of the panel to plot into. Panels are indexed from 0 to _n_ in a row-wise order.
arrays | Tuple (x, y, yerr), np.array | Tuple of numpy arrays: x-values, y-values, and y sigma values
gr | ROOT.TGraphErrors | TGraphErrors object to be plotted
label | str | Label to use for the legend. Use same label for plots for which same legend entry is to be used.
labelLegendId | int | Identifier of the legend to which the plot will be labeled
plotType | str | "data" (default) or "theory". Data curve will be plotted as points and errorbars, while a theory curve will be shown as a colorband.
scale | float | Scale factor for y values of the TGraphErrors object
**plotParams | dict | Supplementary parameters passed to matplotlib `errorbar` or `fill_between` plotting methods depending on plotType.

Ratio plots can be added by calling `plot.Ratio(...)`:

```python
plot.Ratio(r1, r2, style)
```

Parameter | Type | Description
--- | --- | ---
r1 | int | Index of the numerator plot
r2 | int | Index of the denominator plot. The plot style for the ratio will be inherited from this curve.
style | str | Style for the ratio plot. Specify "errorbar" to always have point plot ratio curves

Add systematic uncertainties with `plot.AddSyst(...)`:

```python
plot.AddSyst(r1, ysys)
```

Parameter | Type | Description
--- | --- | ---
r1 | int | Index of the plot to add systematic uncertainty patches
ysys | float or np.array | (Relative) systematic uncertainty. If float, the value will be multiplied by the y-values of the curve, else if numpy array, the values from this array will be used directly.

Finally, create the plot by calling `plot.Plot()`. A plot can be saved with `plot.Save(filename)`, and shown in a window wit `plot.Show()`.

