#!/usr/bin/env python
"""
An example of how to use wx or wxagg in an application with the Basemap module
"""

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

from mpl_toolkits.basemap import Basemap

from wx import *

class CanvasFrame(Frame):

   def __init__(self):
       Frame.__init__(self,None,-1,
                                       'CanvasFrame',size=(550,350))

       self.SetBackgroundColour(NamedColor("WHITE"))

       self.figure = Figure()
       
       self.canvas = FigureCanvas(self, -1, self.figure)
       self.ax = self.figure.add_subplot(111)

       self.sizer = BoxSizer(VERTICAL)
       self.sizer.Add(self.canvas, 1, LEFT | TOP | GROW)
       self.SetSizer(self.sizer)
       self.Fit()

       self.add_toolbar()  # comment this out for no toolbar
       
       self.plot_map()
       
   def add_toolbar(self):
       self.toolbar = NavigationToolbar2Wx(self.canvas)
       self.toolbar.Realize()
       if Platform == '__WXMAC__':
          # Mac platform (OSX 10.3, MacPython) does not seem to cope with
          # having a toolbar in a sizer. This work-around gets the buttons
          # back, but at the expense of having the toolbar at the top
          self.SetToolBar(self.toolbar)
       else:
          # On Windows platform, default window size is incorrect, so set
          # toolbar width to figure width.
          tw, th = self.toolbar.GetSizeTuple()
          fw, fh = self.canvas.GetSizeTuple()
          # By adding toolbar in sizer, we are able to put it at the bottom
          # of the frame - so appearance is closer to GTK version.
          # As noted above, doesn't work for Mac.
          self.toolbar.SetSize(Size(fw, th))
          self.sizer.Add(self.toolbar, 0, LEFT | EXPAND)
       # update the axes menu on the toolbar
       self.toolbar.update()
           
   def plot_map(self):
       map = Basemap(ax=self.ax)
       map.drawcoastlines()
       map.drawcountries()
       map.drawmapboundary()
       map.fillcontinents(color='coral', lake_color='aqua')
       map.drawmapboundary(fill_color='aqua')
       self.figure.canvas.draw()

class App(App):

   def OnInit(self):
       'Create the main window and insert the custom frame'
       frame = CanvasFrame()
       frame.Show(True)
       return True

app = App(0)
app.MainLoop()
