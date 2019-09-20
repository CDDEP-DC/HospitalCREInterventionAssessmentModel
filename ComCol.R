#
# Copyright (C) 2019 Center for Disease Dynamics, Economics and Policy
#
# This file is part of The Hospital CRE Intervention Assessment Model (hCREiAM)
# 
# hCREiAM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# hCREiAM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with hCREiAM  If not, see <https://www.gnu.org/licenses/>.
#
#

################### Environmental Transmission Rate ##########################

enviro = "HET"
minqdp = 0.0001
maxqdp = 0.00046
medqdp = 0.00028

ICOrder = 5

# enviro = "LET"
# minqdp = 0.0003
# maxqdp = 0.00112
# medqdp = 0.00056
# 
# ICOrder = 6

# ################## South Africa Health Center #############################
# 
# # Min Value
# minc = 0.90
# # Max Value
# maxc = 0.9999
# # Median Value
# medc = 0.95
# 
# numpatients = 324
# S  <- c(307/324, 17/324)  # Proportion of susceptible (colonized with susceptible) patient
# 
# minrn.1 = 2
# maxrn.1 = 6
# medrn.1 = 4.04
# 
# minrn.2 = 0.667
# maxrn.2 = 1.333
# medrn.2 = 1
# 
# minrd.1 = 5.62
# maxrd.1 = 9
# medrd.1 = 7.31
# 
# minrd.2 = 1.5
# maxrd.2 = 2.5
# medrd.2 = 2
# 
################## India Health Center #############################

# Min Value
minc = 0.45
# Max Value
maxc = 0.75
# Median Value
medc = 0.60

numpatients = 392
S  <- c(343/392, 49/392)  # Proportion of susceptible (colonized with susceptible) patient

minrn.1 = 3
maxrn.1 = 9
medrn.1 = 6.19

minrn.2 = 1
maxrn.2 = 1.8
medrn.2 = 1.422

minrd.1 = 10
maxrd.1 = 15
medrd.1 = 12.51

minrd.2 = 2.3
maxrd.2 = 3.5
medrd.2 = 2.91

n = 100