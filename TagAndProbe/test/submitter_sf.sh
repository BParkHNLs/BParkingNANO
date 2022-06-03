#!/bin/bash

# ${1} data_file
# ${2} mc_file
# ${3} out_label
# ${4} categorisation

root -l -b -q getScaleFactor.C+\(\"${1}\",\"${2}\",\"${3}\",\"${4}\"\)
