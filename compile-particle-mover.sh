#!/bin/bash

#Files to compile
fortranfiles='create_axis.f90 command_line.f90 write_netcdf-particle.f90 particleMover.f90'

#Compiled output file
outfile='particleMover-out'

#Compiler name
fc=gfortran

#Use nf-config to grab the compile and link flags. Backticks `` run command and grab output
fflags=`nf-config --fflags`
flibs=`nf-config --flibs`
fboth=`nf-config --fflags --flibs`

#Actual compile line. Other flags can be added
#$fc -std=f2008 -c write_netcdf.f90 `nf-config --fflags --flibs`
#$fc -g -std=f2008 -c $fortranfiles `nf-config --fflags --flibs`
$fc -Wall -Wpedantic -g -std=f2008 $fortranfiles -o $outfile $fboth

#Run the executable script:
./$outfile Nx=100 Ny=100 rho=single
