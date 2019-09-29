# example for making plots from our differential equation solver
import ROOT as r      # needed to load ROOT libraries
#import numpy as np    # only needed to use numpy arrays, see example below
import sys
from math import sqrt

# graphs are stored with generic names
# each of our 'n' dependent vars are plotted vs the independent variable (t)
name = "Fastball"
tf=r.TFile(name+".root") # open file for read access

tg_x_vs_t=tf.Get("xy0")   # dependent var[0] vs independent var
tg_vx_vs_t=tf.Get("xy1")  # etc...
tg_y_vs_t=tf.Get("xy2")
tg_vy_vs_t=tf.Get("xy3")
tg_z_vs_t=tf.Get("xy4")
tg_vz_vs_t=tf.Get("xy5")

tc=r.TCanvas()

# it will often be intersting to instead plot one dependent var vs another
time=tg_x_vs_t.GetX()    # extract sample times
xval=tg_x_vs_t.GetY()    # extract array of positions along x axis
height=tg_z_vs_t.GetY()  # extract array of positions along y axis
horizontal = tg_y_vs_t.GetY()
vx=tg_vx_vs_t.GetY()     # extract x,y velocity
vy=tg_vy_vs_t.GetY()
nvals=tg_x_vs_t.GetN()

tg_graph = r.TMultiGraph()
tg_graph.SetTitle(name + " y/z vs x;x [ft];y/z [ft]")

tg_y_vs_x=r.TGraph(nvals,xval,horizontal)   # make a new graph of y vs x!
tg_z_vs_x=r.TGraph(nvals,xval,height)

tg_y_vs_x.SetLineStyle(2);
tg_z_vs_x.SetLineStyle(1);


tg_graph.Add(tg_y_vs_x)
tg_graph.Add(tg_z_vs_x)
tg_graph.Draw("a")
tg_graph.GetYaxis().SetRangeUser(-5,2);




tc.SaveAs(name + ".png")

print("Hit return to exit")
sys.stdin.readline()





