#!/usr/bin/python3
# PROJECT: Wavefront Alignments Algorithms 
# LICENCE: MIT License 
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Plot WFA alignment matrices
# USAGE: python3 wfa2png.py -h

import sys
import copy
import glob
import os.path
import argparse
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns

################################################################################
# WFA-file parsing
################################################################################
def wfa_parse_file(filename):
  # Log
  print('[Parsing]',end='',flush=True)
  # Parse
  wfa_info = {}
  with open(filename) as fp:
    line = fp.readline() # Read line
    while line:
      # Parse data
      if line[0] == '#':
        fields = line.split();
        if fields[1] == "Heatmap":
          matrix_name = fields[2]
          matrix_list = []
          matrix_list.append(fp.readline().strip().split(',')) 
          line = fp.readline()
          while line and line[0] != '#':
            matrix_list.append(line.strip().split(',')) 
            line = fp.readline()
          M = np.matrix(matrix_list).astype(int)
          M = np.where(M==-1,np.nan,M)
          wfa_info[matrix_name] = M
          # Parse next line
          if not line: break
          else: continue
        elif fields[1] == "List":  
          if len(fields) <= 3: 
            wfa_info[fields[2]] = None
          else: 
            CIGAR = np.matrix(fields[3].replace(","," ")).astype(int)
            wfa_info[fields[2]] = CIGAR
        else:    
          wfa_info[fields[1]] = fields[2]
      # Read next line
      line = fp.readline()
  # Process some values
  wfa_info["WFAMode"] = wfa_info["WFAMode"][1:-1] 
  wfa_info["Distance"] = wfa_info["Penalties"][1:-1].split(',')[0]
  # Return
  return wfa_info

################################################################################
# WFA-Heatmap plotting
################################################################################
def wfa_plot_ticks(ax,data,plen,tlen):
  ylen = data.shape[0]
  xlen = data.shape[1]
  # X-ticks
  x_ticks = [i for i in range(0,xlen-1,xlen//5)]
  if x_ticks[-1]+5 < xlen-1: x_ticks.append(xlen-1)  
  else: x_ticks[-1] = xlen-1
  x_ratio = tlen/xlen
  x_ticks_labels = [int(i*x_ratio) for i in x_ticks]
  x_ticks_labels[-1] = tlen-1
  # Y-ticks
  y_ticks = [i for i in range(0,ylen-1,ylen//5)]
  if y_ticks[-1]+5 < ylen-1: y_ticks.append(ylen-1)  
  else: y_ticks[-1] = ylen-1
  y_ratio = plen/ylen
  y_ticks_labels = [int(i*y_ratio) for i in y_ticks]
  y_ticks_labels[-1] = plen-1
  # Set ticks
  ax.set_yticks(y_ticks)
  ax.set_yticklabels(y_ticks_labels,fontsize=5)
  ax.set_xticks(x_ticks)
  ax.set_xticklabels(x_ticks_labels,fontsize=5)

def wfa_plot_wavefront(title,ax,data,plen,tlen,
                       xlabel=False,ylabel=True): 
  # Title
  ax.set_title(title,fontsize=10) 
  # Set ticks
  wfa_plot_ticks(ax,data,plen,tlen)
  # Set labels
  if xlabel: ax.set_xlabel('Text',fontsize=8)
  if ylabel: ax.set_ylabel('Pattern',fontsize=8)
  # Create white grid
  ax.grid(which="major",color="black",linestyle='-',linewidth=0.1)
  # Create colorbar
  cmap = copy.copy(plt.cm.jet)
  cmap.set_bad('whitesmoke')
  # Heatmap
  im = ax.imshow(data,cmap=cmap)
  return im

def wfa_plot_behaviour(title,ax,data,plen,tlen,
                       xlabel=False,ylabel=False): 
  # Title
  ax.set_title(title,fontsize=10) 
  # Set ticks
  wfa_plot_ticks(ax,data,plen,tlen)
  # Set labels
  if xlabel: ax.set_xlabel('Text',fontsize=8)
  if ylabel: ax.set_ylabel('Pattern',fontsize=8)
  # Create white grid
  ax.grid(which="major",color="black",linestyle='-',linewidth=0.1)
  # Create colorbar
  cmap_ext = colors.ListedColormap(['dodgerblue','darkred'])
  cmap_ext_bounds = [10,20]
  cmap_ext_norm = colors.BoundaryNorm(cmap_ext_bounds,cmap_ext.N)
  cmap_ext.set_bad('whitesmoke')
  # Heatmap
  im = ax.imshow(data,cmap=cmap_ext)
  return im

def wfa_cigar(ax,wfa_info,plen,tlen,xlen,ylen):
  # Scale & plot
  def wfa_cigar_scale__plot(cigar,x_ratio,y_ratio,marker,color,label):
    x = np.floor(cigar[:,0].astype('float64') * x_ratio)
    y = np.floor(cigar[:,1].astype('float64') * y_ratio)
    ax.scatter([x],[y],marker=marker,color=color,s=0.1,linewidths=0,label=label)
  # Compute dims
  x_ratio = xlen/tlen  
  y_ratio = ylen/plen
  # Fetch CIGAR
  cigar_m = wfa_info["CIGAR-M"]
  cigar_x = wfa_info["CIGAR-X"]
  cigar_i = wfa_info["CIGAR-I"]
  cigar_d = wfa_info["CIGAR-D"]
  # Plot CIGAR
  if (cigar_m is not None): wfa_cigar_scale__plot(cigar_m,x_ratio,y_ratio,',','limegreen','match')
  if (cigar_x is not None): wfa_cigar_scale__plot(cigar_x,x_ratio,y_ratio,'o','red','misms')
  if (cigar_i is not None): wfa_cigar_scale__plot(cigar_i,x_ratio,y_ratio,'>','orange','ins')
  if (cigar_d is not None): wfa_cigar_scale__plot(cigar_d,x_ratio,y_ratio,'v','blue','del')
  ax.legend(loc="upper right",prop={'size': 5},markerscale=5)
 
def wfa_plot(filename,wfa_info,dpi,compact,extended):
  # Log
  print('[Plotting]',end='',flush=True)
  # Parameters
  plen = int(wfa_info["PatternLength"])
  tlen = int(wfa_info["TextLength"])
  ylen = wfa_info["M"].shape[0]
  xlen = wfa_info["M"].shape[1]
  # Create plot 
  if compact: 
    fig, ax1 = plt.subplots(nrows=1,ncols=1,dpi=dpi,sharex=True)
    im1 = wfa_plot_wavefront('M-Wavefront',ax1,wfa_info["M"],plen,tlen,xlabel=True,ylabel=True)
    #if 'A' in wfa_info["WFAMode"]: wfa_cigar(ax1,wfa_info,plen,tlen,xlen,ylen)
  elif extended or wfa_info["Distance"]=="Edit":
    fig, (ax1,ax3) = plt.subplots(nrows=2,ncols=1,dpi=dpi,sharex=True)
    im1 = wfa_plot_wavefront('M-Wavefront',ax1,wfa_info["M"],plen,tlen,ylabel=True)
    im3 = wfa_plot_behaviour('Extend & CIGAR',ax3,wfa_info["Extend"],plen,tlen,xlabel=True,ylabel=True)
    if 'A' in wfa_info["WFAMode"]: wfa_cigar(ax3,wfa_info,plen,tlen,xlen,ylen)
  elif wfa_info["Distance"]=="GapLineal" or wfa_info["Distance"]=="GapAffine":
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,dpi=dpi,sharex=True)
    im1 = wfa_plot_wavefront('M-Wavefront',ax1,wfa_info["M"],plen,tlen,ylabel=True)
    im2 = wfa_plot_wavefront('I1-Wavefront',ax2,wfa_info["I1"],plen,tlen)
    im4 = wfa_plot_wavefront('D1-Wavefront',ax4,wfa_info["D1"],plen,tlen,xlabel=True)
    im3 = wfa_plot_behaviour('Extend & CIGAR',ax3,wfa_info["Extend"],plen,tlen,xlabel=True,ylabel=True)
    if 'A' in wfa_info["WFAMode"]: wfa_cigar(ax3,wfa_info,plen,tlen,xlen,ylen)
  elif wfa_info["Distance"]=="GapAffine2p":    
    fig, ((ax1,ax2,ax5),(ax3,ax4,ax6)) = plt.subplots(nrows=2,ncols=3,dpi=dpi,sharex=True)
    im1 = wfa_plot_wavefront('M-Wavefront',ax1,wfa_info["M"],plen,tlen,ylabel=True)
    im2 = wfa_plot_wavefront('I1-Wavefront',ax2,wfa_info["I1"],plen,tlen)
    im4 = wfa_plot_wavefront('D1-Wavefront',ax4,wfa_info["D1"],plen,tlen,xlabel=True)
    im5 = wfa_plot_wavefront('I2-Wavefront',ax5,wfa_info["I2"],plen,tlen)
    im6 = wfa_plot_wavefront('D2-Wavefront',ax6,wfa_info["D2"],plen,tlen,xlabel=True)
    im3 = wfa_plot_behaviour('Extend & CIGAR',ax3,wfa_info["Extend"],plen,tlen,xlabel=True,ylabel=True)
    if 'A' in wfa_info["WFAMode"]: wfa_cigar(ax3,wfa_info,plen,tlen,xlen,ylen)
  # Color bar
  if compact: 
    p0 = ax1.get_position().get_points().flatten()
    p1 = p0
  elif extended or wfa_info["Distance"]=="Edit":
    p0 = ax1.get_position().get_points().flatten()
    p1 = p0
  elif wfa_info["Distance"]=="GapLineal" or wfa_info["Distance"]=="GapAffine":
    p0 = ax3.get_position().get_points().flatten()
    p1 = ax4.get_position().get_points().flatten()
  elif wfa_info["Distance"]=="GapAffine2p":    
    p0 = ax3.get_position().get_points().flatten()
    p1 = ax6.get_position().get_points().flatten()
  ax_cbar = fig.add_axes([p0[0],0,p1[2]-p0[0],0.025])
  ax_cbar.tick_params(labelsize=5) 
  plt.colorbar(im1,cax=ax_cbar,orientation='horizontal')
  # Title
  file = os.path.basename(filename).replace('.wfa','')
  title = "WFA-Plot(%s) \n %s[%s]" % (file,wfa_info["Penalties"],wfa_info["WFAMode"])
  fig.suptitle(title,fontsize=12)
  plt.subplots_adjust(top=0.85)
  # Plot
  plt.savefig(filename.replace('.wfa','.png'),bbox_inches='tight')

################################################################################
# Main
################################################################################
# Configure arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', action='store', help='Input file')
parser.add_argument('--dpi', type=int, action='store', default=1500, help='Plot resolution') # More than 2000 is hard to handle
parser.add_argument('--compact', action='store_true', default=False, help='Plot M-Wavefront only')
parser.add_argument('--extended', action='store_true', default=True, help='Plot M-Wavefront and extension/CIGAR')
parser.add_argument('--full', action='store_true', default=False, help='Plot all info available')
parser.add_argument('-H', action='store_true', dest="human_readable", default=False)

# Parse arguments
args = parser.parse_args()

# Open input
if args.input:
  input_files = [args.input]
else:
  input_files = glob.glob("*.wfa")
  print('[WFA2png] Searching all *.wfa ( Found %d file%c)' % (len(input_files),'s' if len(input_files)>1 else ' '))
  
# Plot each WFA file
print('[WFA2png] Plotting at %d dpi' % (args.dpi))
idx = 0
for filename in input_files:
  print('[WFA2png] [#%d] Generating \'%s\' ' % (idx,filename),end='',flush=True)
  wfa_info = wfa_parse_file(filename)
  wfa_plot(filename,wfa_info,args.dpi,args.compact,args.extended)
  print('[Done!]',)
  idx += 1



  
  
