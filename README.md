# JPyPlotRatio
Multipanel plotting class with ratio panels.

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
rowBounds | Dict `{rowN:(ymin,ymax), ...}` | Dictionary of y-limits for each row
colBounds | Dict `{colN:(xmin,xmax), ...}` | Dictionary of x-limits for each column
ratioBounds | Dict `{rowN:(ymin,ymax), ...}` | Dictionary of y-limits for the ratio panels in each row
panelScaling | Dict `{panelIndex:scale, ...}` | Dictionary of scale factors for the plots in each panel
panelLabel | Dict `{panelIndex:label(str), ...}` | Dictionary of panel labels
panelLabelLoc | Tuple (x, y) | Location for the panel labels in each panel. Default `(0.2,0.92)`
panelLabelSize | int | Panel label text size. Default 16
panelLabelAlign | str | Text alignment for the panel labels, "left", "center", "right" (default).
legendLoc | Tuple (x, y) | Legend location in the panel
legendSize | int | Legend text size. Default 10
systPatchWidth | float | Width of the systematic uncertainty patches. Default 5.0.

Plot curves with `plot.Add(...)` or `plot.AddTGraph(...)`. `AddTGraph` automatically converts a TGraphErrors object to numpy arrays and plots them.

```python
plotIndex = plot.Add(panelIndex, (x, y, yerr), label, plotType, **plotParams)
plotIndex = plot.AddTGraph(panelIndex, gr, label, plotType, **plotParams)
```

Returns the index for the newly added plot (int).

Parameter | Type | Description
--- | --- | ---
panelIndex | int | Index of the panel to plot into. Panels are indexed from 0 to _n_ in a row-wise order.
arrays | Tuple (x, y, yerr), np.array | Tuple of numpy arrays: x-values, y-values, and y sigma values
gr | ROOT.TGraphErrors | TGraphErrors object to be plotted
label | str | Label to use for the legend. Use same label for plots for which same legend entry is to be used.
plotType | str | "data" (default) or "theory". Data curve will be plotted as points and errorbars, while a theory curve will be shown as a colorband.
**plotParams | dict | Supplementary parameters passed to matplotlib `errorbar` or `fill_between` plotting methods depending on plotType.

Ratio plots can be added by calling `plot.Ratio(...)`:

```python
plot.Ratio(r1, r2)
```

Parameter | Type | Description
--- | --- | ---
r1 | int | Index of the numerator plot
r2 | int | Index of the denominator plot. The plot style for the ratio will be inherited from this curve.

Add systematic uncertainties with `plot.AddSyst(...)`:

```python
plot.AddSyst(r1, ysys)
```

Parameter | Type | Description
--- | --- | ---
r1 | int | Index of the plot to add systematic uncertainty patches
ysys | float or np.array | (Relative) systematic uncertainty. If float, the value will be multiplied by the y-values of the curve, else if numpy array, the values from this array will be used directly.

Finally, create the plot by calling `plot.Plot()`. A plot can be saved with `plot.Save(filename)`, and shown in a window wit `plot.Show()`.

