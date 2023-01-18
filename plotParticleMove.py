#!/usr/bin/env /usr/local/anaconda3/bin/python3
#
#Import modules
import netCDF4
import csv
import sys
import time
import shutil
import os.path as path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from io import StringIO
from cycler import cycler
from argparse import ArgumentParser
from matplotlib.ticker import FormatStrFormatter

mpl.use('Qt5Agg')
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}
plt.rc('font', **font)
#mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}' #for \text command
plt.rcParams.update({
    #"text.latex.preamble": r'\usepackage{amsmath}',
    #"text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})
# for Palatino and other serif fonts use:
plt.rcParams.update({
    #"text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Palatino"],
})
# It's also possible to use the reduced notation by directly setting font.family:
plt.rcParams.update({
  #"text.usetex": True,
  "font.family": "Helvetica"
})

parser = ArgumentParser(description='Script to plot file with numerical results')
parser.add_argument("-f", "--file", dest='dataFile', default="dataFile.dat", action='store', type=str,
                    help="Read data from FILE", metavar="FILE")
parser.add_argument("-c", "--caseId", dest='CaseId', default="Case Id", action='store', type=str,
                    help="String to be used as Case Identifier", metavar="CASE_ID")
# - #parser.add_argument("-q", "--quiet",
# - #                    action="store_false", dest="verbose", default=True,
# - #                    help="don't print status messages to stdout")
# - #
args = parser.parse_args()

print("sys.float_info.dig = ", sys.float_info.dig)
print("sys.float_info.mant_dig = ", sys.float_info.mant_dig)
epsMachine = np.finfo(float).eps*np.float64(10.0)
print("epsMachine=",epsMachine)
# Definicion de constantes
pi = np.float64(4.0)*np.arctan(1.0)

#Caso Id
#CaseId = sys.argv[1]
CaseId = args.CaseId
caseIdLabel = CaseId.strip().replace(" ", "")

# Lectura de datos entrada
dataFile = args.dataFile
dataFileName = path.basename(dataFile)
dataFileDir  = path.dirname(dataFile)
dataFileNoExt = path.splitext(dataFileName)[0]
print(dataFileNoExt)
nHeaderLines=0;
#if (dataFile == "forceCoeffs.dat"):
#    nHeaderLines=0
#elif (dataFile == "coefficient.dat"):
#    nHeaderLines=12
#else:
#    print("\n WARNING \n File name not recognised: ... ", dataFileName, "\n -------- \n")
#    parser.print_help()

#  # Read initial lines - Unstructured data
#  with open(dataFile) as f:
#      reader = csv.reader(f);
#      if (nHeaderLines != 0):
#          header = [ next(reader) for x in range(nHeaderLines) ]
#  
#  # preProcess file to avoid problems with special characters (, #, ) 
#  bufferData = StringIO(open(dataFile).read().replace('(','').replace(')','').replace("# ",""))
#  
#  with open('checkFile.txt', 'w') as fd:
#    bufferData.seek(0)
#    shutil.copyfileobj(bufferData, fd)
#  
## ----------------------------------------------------------------
# Read structured data
#bufferData.seek(0)
simData = pd.read_csv(dataFileName, delim_whitespace=True, header=None)
simData.to_numpy()
#columns = simData.columns
#simData = simData.iloc[:, :-1]
#simData.columns = columns[1:]
## ----------------------------------------------------------------
## Constructing plots
r=simData.iloc[:,0]
print(r)
analytic=simData.iloc[:,1]
numeric=simData.iloc[:,2]
absError=simData.iloc[:,3]
relError=simData.iloc[:,4]

fig, ax1 = plt.subplots(1,1,figsize=(14.2,7.5))
fig.subplots_adjust(right=0.8)
#lineCd = ax.plot(simData['Time'],simData['Cd'],label=labelText,markevery=2)
#ax1 = simData.plot(x="Time",y="Cd",color=color);

color = 'tab:red'
p1, = ax1.plot(r,analytic,'x',ls='-',color=color,linewidth=1.2,markevery=3,label='Analytic');
ax1.tick_params(axis='y')#,labelcolor=color)
ax1.set_xlabel('Pipe radius, $r(m)$')
ax1.set_ylabel('Flow velocity, $U_z\, (m/s)$')
color = 'tab:blue'
p2, = ax1.plot(r,numeric,'+',ls='-',color=color,linewidth=1.2,markevery=3,label='Numeric');
legend = ax1.legend();
#ax2.tick_params(axis='y', labelcolor=color)
#ax2.set_ylabel('Numerical solution', color=color)  # we already handled the x-label with ax1
# -- #ax1.set_yticks(np.linspace(ax1.get_ybound()[0], ax1.get_ybound()[1], 5))
# -- #ax2.set_yticks(np.linspace(ax2.get_ybound()[0], ax2.get_ybound()[1], 5))
ax1.axis(ymin=0,ymax=1.1*np.max(numeric))
ax1.axis(xmin=1.05*np.min(r),xmax=1.05*np.max(r))
ax1.set_yticks(np.linspace(0,1.1*np.max(numeric),10)) 
ax1.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
ax1.yaxis.set_major_formatter(mticker.EngFormatter(unit='m/s',places=2,sep=' '))
ax1.xaxis.set_major_formatter(mticker.EngFormatter(unit='m',places=2,sep=' '))
#ax2.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
#ax2.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.3f'))
# -- #ax2.yaxis.set_major_formatter(FormatStrFormatter('%.3f'))

#  meanCd = 0.1;
#  stdCd = 0.2;
#  meanCl = 0.3;
#  stdCl = 0.4;
#  meanCm = 0.5;
#  stdCm = 0.6;
#  textstr = '\n'.join((
#      r'$\mu_{Cd}=%.3e$' % (meanCd, ),
#      r'$\sigma_{Cd}=%.3e$' % (stdCd, ),
#      r'$\mu_{Cl}=%.3e$' % (meanCl, ),
#      r'$\sigma_{Cl}=%.3e$' % (stdCl, ),
#      r'$\mu_{CmPitch}=%.3e$' % (meanCm, ),
#      r'$\sigma_{CmPitch}=%.3e$' % (stdCm, ),
#      ))
# # these are matplotlib.patch.Patch properties
# props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
# # place a text box in upper left in axes coords
# ax1.text(0.80, 0.95, textstr, transform=ax1.transAxes, fontsize=14,
#         verticalalignment='top', bbox=props)

plt.grid();
fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.95)
plt.margins(0,0)
plt.title(CaseId)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig(dataFileNoExt+'_'+caseIdLabel+'.png',bbox_inches='tight',pad_inches=0.5)
plt.show()

# ---------  TRASH/WASTE/PAST COMMANDS

# ---------------------------------------------------------------
#simData_downSample_mean = simData.resample('25U', on='Time').mean()
#print(simData_downSample_mean.head(8))
#simData['Cd'].plot(x='Time',alpha=0.5,style='-')
#simData_downSample_mean['Cd'].plot(x='Time',style=':')

#print(simData_downSample_mean)

#xMin = np.float64(simData.loc['xMin']['Value'])
#xMax = np.float64(simData.loc['xMax']['Value'])
#alpha  = np.float64(simData.loc['alpha']['Value'])
#phiXMin= np.float64(simData.loc['phiXMin']['Value'])
#phiXMax= np.float64(simData.loc['phiXMax']['Value'])
#simTime=  np.float64(simData.loc['simTime']['Value'])
#deltaX = np.float64(simData.loc['deltaX']['Value'])
#deltaT = np.float64(simData.loc['deltaT']['Value'])
#nPlots = np.float64(simData.loc['nPlots']['Value'])
