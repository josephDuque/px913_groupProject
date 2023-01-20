#!/bin/bash

#Files to compile
fortranfiles='create_axis.f90 command_line.f90 write_netcdf-particle.f90 particleMover.f90'

#Compiled output file
outfile='particleMover-out'

#Arguments for the executable:
input1=$1
input2=$2
input3=$3

# Verifying correct input of arguments, and checking which inputs
# correspond to which arguments, by comparing substrings:
#------------------Checking array size argument Nx---------------------
if [[ "$input1" == *"Nx="* ]];
then
    nxArg=$input1
elif [[ "$input2" == *"Nx="* ]];
then
    nxArg=$input2
elif [[ "$input3" == *"Nx="* ]];
then
    nxArg=$input3
else
    echo "Please input 'Nx=<value>' to specify array row size"
fi
#--------------------Checking array size argument Ny---------------------
if [[ "$input1" == *"Ny="* ]];
then
    nyArg=$input1
elif [[ "$input2" == *"Ny="* ]];
then
    nyArg=$input2
elif [[ "$input3" == *"Ny="* ]];
then
    nyArg=$input3
else
    echo "Please input 'Ny=<value>' to specify array column size"
fi
#---------------------Checking initialiser argument-------------------
if [[ "$input1" == *"rho="* ]];
then
    initArg=$input1
elif [[ "$input2" == *"rho="* ]];
then
    initArg=$input2
elif [[ "$input3" == *"rho="* ]];
then
    initArg=$input3
else
    echo "Please input 'rho=<'init'>' to specify initalisation state. 3 are possible: null, single or double"
fi

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
./$outfile $nxArg $nyArg $initArg

python3 plotParticleMove.py
