#!/bin/bash

# ${1} data_file
# ${2} mc_file
# ${3} out_label
# ${4} categorisation
# ${5} do_plots

root -l -b -q getScaleFactor.C+\(\"${1}\",\"${2}\",\"${3}\",\"${4}\"\)

if [ ${5} == "True" ] ; then
  root -l -b -q savePlots.C+\(\"${1}\",\"${3}\",\"False\"\)
  root -l -b -q savePlots.C+\(\"${2}\",\"${3}\",\"True\"\)
fi

