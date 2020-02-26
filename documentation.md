---
description: 'All functions, inputs and outputs. Some functions borrowed from PyEIS'
---

# Documentation

## mpt\_data

**init**\(self, path, data, cycle='off', mask=\['none','none'\], gph\_width = 6.4, gph\_height = 4.8\)

* inputs: 
  * path must be a string that leads into the folder that has the file in question
  * data must be a list that has the string or strings of the name\(s\) of the file\(s\) in question
  * cycle can be switched on and off if the mpt file has multiple cycles in it. It is off by default
  * mask must be a list of two integers that determine the mask boundary of the frequency in the mpt file
  * gph\_width and gph\_width serve as the graph width and height of the graph plot. This can be changed later
* output:
  * mpt\_data object

**set\_gph\_width**\(self, new\_width\)

* inputs:
  * new\_width takes an integer that will serve as the width for the graph when asked to be plot
* outputs:
  * None

**set\_gph\_height**\(self, new\_height\)

* inputs:
  * new\_height takes an integer that will serve as the height for the graph when asked to be plot
* outputs:
  * None

**set\_new\_gph\_dims**\(self, new\_width,  new\_height\)

* inputs:
  * employs both setting new height and width of the graph at once
* outputs:
  * None

**mpt\_plot**\(self, fitting='off', rr='off', x\_window = 'none', y\_window = 'none'\)

* inputs
  * fitting can either be 'on' or 'off'; this can only work after the fitting function has been performed. The fitting will overlay the Nyvquist plot
  * rr can be 'on' or 'off'; this also only works after fitting function has been performed. The residuals throughout the graph will be plotted below the Nyvquist plot and will show the difference between the predicted values and the actual values
  * x\_window and y\_window takes in a list of two integers which will determine the window of which our graph will look into
* outputs
  * 



