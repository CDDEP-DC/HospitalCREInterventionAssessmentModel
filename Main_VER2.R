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

### ENSURE THAT YOU SET Working directory to file location
setwd("/Users/cddepadmin/Documents/IPC HAI CRE Project Code/HospitalCREInterventionAssessmentModel")

rm(list=ls()) # clear
set.seed(1982)
options(stringsAsFactors=F)
if (!("plyr" %in% rownames(installed.packages()))){
  install.packages("plyr", quiet = TRUE)
}
library(plyr)

if (!("readxl" %in% rownames(installed.packages()))){
  install.packages("readxl", quiet = TRUE)
}
library(readxl)

if (!("RColorBrewer" %in% rownames(installed.packages()))){
  install.packages("RColorBrewer", quiet = TRUE)
}
library(RColorBrewer)

if (!("deSolve" %in% rownames(installed.packages()))){
  install.packages("deSolve", quiet = TRUE)
}
library(deSolve)

if (!("data.table" %in% rownames(installed.packages()))){
  install.packages("data.table", type = "binary", quiet = TRUE)
}
library(data.table)

if (!("ggplot2" %in% rownames(installed.packages()))){
  install.packages("ggplot2", quiet = TRUE)
}
library(ggplot2)

if (!("rms" %in% rownames(installed.packages()))){
  install.packages("rms", quiet = TRUE)
}
library(rms)

if (!("triangle" %in% rownames(installed.packages()))){
  install.packages("triangle", quiet = TRUE)
}
library(triangle)

if (!("lhs" %in% rownames(installed.packages()))){
  install.packages("lhs", quiet = TRUE)
}
library(lhs)


## function to convert tibble read in format to comma separated list and call triangleLHS and uniformLHS
triangleLHS_input <- function(n, vals){
	v=as.numeric(unlist(vals))
	triangleLHS(n, v[1], v[2], v[3])    
}

triangleLHS <- function(n, min, max, mode = (min + max)/2){
    p <- randomLHS(n, 1)[,1]
    qtriangle(p = p, a = min, b = max, c = mode)
}

uniformLHS_input <- function(n, vals){
  v=as.numeric(unlist(vals))
  uniformLHS(n, v[1], v[2])    
}

uniformLHS <- function(n, min, max){
  p <- randomLHS(n, 1)[,1]
  qunif(p = p, a = min, b = max)
}

#Source the model code
source("Model.R")
source("prcc.R")

## Defining the start time of running the model to help organize files
runstartday = format(Sys.time(), "%Y%m%d_%H%M")

## Output files are saved to this folder, "~/hCREiAM Data (Start Time)"
dir.create(paste("hCREiAM Data",runstartday,""))

## Importing the Healthcare Facility specific parameters

HCFparamvals <- read_excel("hCREiAM_IPC_Parameters.xlsx",sheet="HealthCareFacility")
hpv = HCFparamvals[1:36,]

## The two following variables change the duration of the simulation and are defined
## at the end of the HealthCareFacility spreadsheet

## The total number of partitions in the Latin Hypercube Sample
n = as.numeric(HCFparamvals[which(HCFparamvals$Parameter=="Number of Runs"),2])

## The total number of days a single iteration of the simulation will run for
days = as.numeric(HCFparamvals[which(HCFparamvals$Parameter=="Days"),2])

## Importing cost data from the excel spreadsheet
CSTparamvals <- read_excel("hCREiAM_IPC_Parameters.xlsx",sheet="Cost")

## Each entry in the following for loop corresponds with a different healthcare facility
## The parameters ddefining the HCF are defined first

for(envnum in 7:10) {

  hpvrun = matrix(NA, nrow(hpv), 5)
  colnames(hpvrun) <- colnames(hpv[1:5])

   for (g in 1:nrow(hpv)){
    if(isTRUE(hpv[g,envnum]==1)) {
      hpvrun[g,1:5] <- unlist(hpv[g,1:5])
    }
  }
  hpvrun <- as.data.frame(hpvrun[rowSums(is.na(hpvrun)) != ncol(hpvrun),])

  rn.1 = triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="rn.1"),2:4])
  rn.2 = triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="rn.2"),2:4])
  rd.1 = triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="rd.1"),2:4])
  rd.2 = triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="rd.2"),2:4])
  
  Qnp.lhs = triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="qnp"),2:4])
  Qdp.lhs = triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="qdp"),2:4])
  
  ICOrder = triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="ICOrder"),2:4])
  
  numbeds = triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="numbeds"),2:4])  # Total number of patients in each ward
  
  # NUMBER OF WARDS - The first number in each parenthesis below is the
  # value for ward 1, the second number will be the value for ward 2.  All
  # parentheses MUST have the same number of values!
  
  S.1 = mean(triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="S.1"),2:4]))
  S.2 = mean(triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="S.2"),2:4]))
  S  <- c(S.1, S.2)     # Proportion of patients susceptible to pathogen 
  
  Cu <- c(0, 0)         # Proportion of patients colonized with resistant pathogen and it's undetected
  Cd <- c(0, 0)         # Proportion of patients colonized with resistant pathogen and it's detected
  I  <- c(0, 0)         # Proportion of patients infected with resistant pathogen
  XR <- c(0, 0)         # Proportion of patients successfully prophylaxed against resistant pathogen
  
  Y  <- c(0, 0)         # Environmental contamination on surfaces (excluding water)
  Tx <- c(0, 0)         # Environmental contamination on textiles
  M  <- c(0, 0)         # Environmental contamination on medical textiles
  L  <- c(0, 0)         # Environmental contamination on mobile medial textiles
  W  <- c(0, 0)         # Environmental contamination on water sources
  
  Cd1<- c(0, 0)         # Proportion of patients colonized with resistant pathogen and it's detected NOT COHORTED
  Cd2<- c(0, 0)         # Proportion of patients colonized with resistant pathogen and it's detected COHORTED
  I1 <- c(0, 0)         # Proportion of patients infected with resistant pathogen NOT COHORTED
  I2 <- c(0, 0)         # Proportion of patients infected with resistant pathogen COHORTED
  
  # The following parameters correspond with outputs that are important to capture
  
  CH <- c(0, 0)         # Number of acummulated colonized patients created inside the healthcare facility
  CI <- c(0, 0)         # Number of accumalated colonized patients imported into the hospital
  
  IH <- c(0, 0)         # Number of acummulated infected patients created inside the healthcare facility
  II <- c(0, 0)         # Number of accumalated infected patients imported into the hospital
  
  Ie <- c(0, 0)         # Number of acummulated infected patients through the environment
  Ic <- c(0, 0)         # Number of accumalated infected patients through contacts
  
  Deaths <- c(0, 0)     # Total number of deaths
  
  # Intial state puts all of these values into a data frame
  InitialState <- as.data.frame(cbind(S, XR, Cu, Cd, I, Y, Tx, M, L, W, CH, CI, IH, II, Ie, Ic, Deaths))
  InStateCHRT  <- as.data.frame(cbind(S, XR, Cu, Cd1, Cd2, I1, I2, Y, Tx, M, L, W, CH, CI, IH, II, Ie, Ic, Deaths))
  
  patients.1 <- numbeds*sum(InitialState[1,1:5])  # Total number of patients in the General Ward
  patients.2 <- numbeds*sum(InitialState[2,1:5])  # Total number of patients in the ICU
  
  #######################################################
  # LHS 
  #######################################################
  
  # Need this for all reps by group
  aS.1  <- triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="aS.1"),2:4]) 
  aS.2  <- triangleLHS_input(n, hpvrun[which(hpvrun$Parameter=="aS.2"),2:4])
  aC.1  <- (1-aS.1)
  aC.2  <- (1-aS.2)
  aI.1  <- rep(0, n)
  aI.2  <- rep(0, n)
  aXR.1 <- rep(0, n)
  aXR.2 <- rep(0, n)
  
  # Read in the baseline parameters -- Move parameters spreadsheet name and associated parameters to the top to initialize
  baselineparamvals <- read_excel("hCREiAM_IPC_Parameters.xlsx",sheet="Baseline")
  bpv = baselineparamvals
  
  # Set the values based on read in values
  omega.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="omega.lhs"),2:4])
  rhoO.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="rhoO.lhs"),2:4])
  pi.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="pi.lhs"),2:4])
  phi.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="phi.lhs"),2:4])
  rhoR.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="rhoR.lhs"),2:4])
  epsilon.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="epsilon.lhs"),2:4])
  psi.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="psi.lhs"),2:4])
  theta.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="theta.lhs"),2:4])
  tau.1 = triangleLHS_input(n, bpv[which(bpv$Parameter=="tau.1"),2:4])
  tau.2 = triangleLHS_input(n, bpv[which(bpv$Parameter=="tau.2"),2:4])
  zeta.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="zeta.lhs"),2:4])
  sigmaD.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="sigmaD.lhs"),2:4])
  sigmaN.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="sigmaN.lhs"),2:4])
  delta.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="delta.lhs"),2:4])
  xi.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="xi.lhs"),2:4])
  eta.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="eta.lhs"),2:4])
  mu.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="mu.lhs"),2:4])
  aP = triangleLHS_input(n, bpv[which(bpv$Parameter=="aP"),2:4])
  bP = triangleLHS_input(n, bpv[which(bpv$Parameter=="bP"),2:4])
  sigY.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="sigY.lhs"),2:4])
  sigT.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="sigT.lhs"),2:4])
  sigMc.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="sigMc.lhs"),2:4])
  sigMsc.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="sigMsc.lhs"),2:4])
  sigMnc.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="sigMnc.lhs"),2:4])
  sigL.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="sigL.lhs"),2:4])
  sigW.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="sigW.lhs"),2:4])
  betaM.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="betaM.lhs"),2:4])
  betaL.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="betaL.lhs"),2:4])
  betaY.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="betaY.lhs"),2:4])
  betaT.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="betaT.lhs"),2:4])
  betaW.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="betaW.lhs"),2:4])
  stay.1 = triangleLHS_input(n, bpv[which(bpv$Parameter=="stay.1"),2:4])
  stay.2 = triangleLHS_input(n, bpv[which(bpv$Parameter=="stay.2"),2:4])
  trans.1 = triangleLHS_input(n, bpv[which(bpv$Parameter=="trans.1"),2:4])
  trans.2 = triangleLHS_input(n, bpv[which(bpv$Parameter=="trans.2"),2:4])
  Cnp.1 = triangleLHS_input(n, bpv[which(bpv$Parameter=="Cnp.1"),2:4])
  Cnp.2 = triangleLHS_input(n, bpv[which(bpv$Parameter=="Cnp.2"),2:4])
  Cdp.1 = triangleLHS_input(n, bpv[which(bpv$Parameter=="Cdp.1"),2:4])
  Cdp.2 = triangleLHS_input(n, bpv[which(bpv$Parameter=="Cdp.2"),2:4])
  Qpn.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="Qpn.lhs"),2:4])
  Qpd.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="Qpd.lhs"),2:4])
  ICfac.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="IC.lhs"),2:4])
  
  # Update the baseline parameters not set in the excel spreadsheet
  Pc.lhs      = aP*bP      
  Psc.lhs     = aP*(1 - bP)      
  Pnc.lhs     = 1 - aP      
  
  IC.lhs  = ICfac.lhs*10^ICOrder 
  
  # Now read in the intervention parameters at baseline and for each intervention
  interventionparamvals <- read_excel("hCREiAM_IPC_Parameters.xlsx",sheet="Interventions",skip=1)
  ipv = interventionparamvals
  
  #first the baselines
  sigD.1    = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD.1"),2:4])
  sigD.2    = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD.2"),2:4])
  PS.1      = triangleLHS_input(n, ipv[which(ipv$Parameter=="PS.1"),2:4])
  PS.2      = triangleLHS_input(n, ipv[which(ipv$Parameter=="PS.2"),2:4])
  OM.1      = triangleLHS_input(n, ipv[which(ipv$Parameter=="OM.1"),2:4])
  OM.2      = triangleLHS_input(n, ipv[which(ipv$Parameter=="OM.2"),2:4])
  sigCP.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigCP.1"),2:4])
  sigCP.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigCP.2"),2:4])
  sigN1.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigN1.1"),2:4])
  sigN1.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigN1.2"),2:4])
  sigD1.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD1.1"),2:4])
  sigD1.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD1.2"),2:4])
  sigY1.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigY1.1"),2:4])
  sigY1.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigY1.2"),2:4])
  sigT1.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigT1.1"),2:4])
  sigT1.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigT1.2"),2:4])
  sigL1.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigL1.1"),2:4])
  sigL1.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigL1.2"),2:4])
  sigW1.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigW1.1"),2:4])
  sigW1.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigW1.2"),2:4])
  sigMc1.1  = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMc1.1"),2:4])
  sigMc1.2  = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMc1.2"),2:4])
  sigMsc1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMsc1.1"),2:4])
  sigMsc1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMsc1.2"),2:4])
  sigMnc1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMnc1.1"),2:4])
  sigMnc1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMnc1.2"),2:4])
  sigAe.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAe.1"),2:4])
  sigAe.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAe.2"),2:4])
  sigAm.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAm.1"),2:4])
  sigAm.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAm.2"),2:4])
  sigAp.1   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAp.1"),2:4])
  sigAp.2   = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAp.2"),2:4])
  f.1       = triangleLHS_input(n, ipv[which(ipv$Parameter=="f.1"),2:4])
  f.2       = triangleLHS_input(n, ipv[which(ipv$Parameter=="f.2"),2:4])
  Fn.1      = triangleLHS_input(n, ipv[which(ipv$Parameter=="Fn.1"),2:4])
  Fn.2      = triangleLHS_input(n, ipv[which(ipv$Parameter=="Fn.2"),2:4])
  
  #then the intervention vals
  # Chlorhexidine Bathing
  sigD.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD.1"),5:7])
  sigD.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD.2"),5:7])
  
  # Active Surveillance
  PS.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="PS.1"),5:7])
  PS.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="PS.2"),5:7])
  PS.2.ICU = triangleLHS_input(n, ipv[which(ipv$Parameter=="PS.ICU.2"),5:7])
  OM.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="OM.1"),5:7])
  OM.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="OM.2"),5:7])
  OM.2.ICU = triangleLHS_input(n, ipv[which(ipv$Parameter=="OM.ICU.2"),5:7])
  
  # Contact Precautions
  sigCP.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigCP.1"),5:7])
  sigCP.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigCP.2"),5:7])
  
  # Hand Hygiene
  sigN1.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigN1.1"),5:7])
  sigN1.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigN1.2"),5:7])
  sigD1.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD1.1"),5:7])
  sigD1.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD1.2"),5:7])
  
  # Enhanced Cleaning
  sigY1.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigY1.1"),5:7])
  sigY1.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigY1.2"),5:7])
  sigT1.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigT1.1"),5:7])
  sigT1.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigT1.2"),5:7])
  sigL1.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigL1.1"),5:7])
  sigL1.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigL1.2"),5:7])
  sigW1.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigW1.1"),5:7])
  sigW1.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigW1.2"),5:7])
  
  sigMc1.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMc1.1"),5:7])
  sigMc1.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMc1.2"),5:7])
  sigMsc1.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMsc1.1"),5:7])
  sigMsc1.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMsc1.2"),5:7])
  sigMnc1.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMnc1.1"),5:7])
  sigMnc1.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMnc1.2"),5:7])
  
  #  Judicious Use of Antibiotics
  sigAe.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAe.1"),5:7])
  sigAe.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAe.2"),5:7])
  sigAm.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAm.1"),5:7])
  sigAm.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAm.2"),5:7])
  sigAp.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAp.1"),5:7])
  sigAp.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAp.2"),5:7])
  
  #  Cohorting
  f.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="f.1"),5:7])
  f.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="f.2"),5:7])
  Fn.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="Fn.1"),5:7])
  Fn.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="Fn.2"),5:7])
  
  # Create a folders to output all model and plot files for a specific environment
  
  dir.create(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),sep="/"))
  dir.create(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),"PRCCValues",sep="/"))
  dir.create(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),"Plots",sep="/"))
  dir.create(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),"CostTables",sep="/"))
  
  ## This is where the Latin Hypercube Sample Begins
  
  for(intnum in 8:21) {
  	LHSResults = matrix(NA,nrow=n,ncol=131)
  	cnames = c("aS", "aXR", "aC", "aI", "Cnp", "Cdp", "rn", "rd", "patients", "stay", "trans", "tau")
  	anames = c("omega", "rhoO", "pi", "phi", "rhoR", "epsilon", "psi", 'theta', "zeta",
  	           "sigmaD", "sigmaN", "delta", "Qpn", "Qpd", "Qnp", "Qdp", "xi", "eta", "mu",
  	           "Pc", "Psc", "Pnc", "sigY", "sigT", "sigMc", "sigMsc", "sigMnc", "sigL", "sigW",
  	           "betaM", "betaL", "betaY", "betaT", "betaW", "IC")
  	inames = c("sigD", "PS", "OM", "sigCP", "sigN1", "sigD1", "sigY1", "sigT1", "sigL1",
  	           "sigW1", "sigMc1", "sigMsc1", "sigMnc1", "sigAe", "sigAm", "sigAp", "f", "Fn")
  	rnames = c("S", "XR", "Cu", "Cd", "I", "Y", "Tx", "M", "L", "W", "CH", "CI","IH","II","Ie","Ic","Deaths")
  	dnames = c("IDays")
  	colnames(LHSResults) <- c(paste(cnames,".1",sep=""),paste(cnames,".2",sep=""),
  	                          anames,
  	                          paste(inames,".1",sep=""),paste(inames,".2",sep=""),
  	                          paste(rnames,".1",sep=""),paste(rnames,".2",sep=""),
  	                          paste(dnames,".1",sep=""),paste(dnames,".2",sep=""))

  	## Adjusting rnames for the cohorting wards
  	if(intnum %in% 19:21){
  	  LHSResults = matrix(NA,nrow=n,ncol=135)
  	  rnames = c("S", "XR", "Cu", "Cd1", "Cd2", "I1", "I2", "Y", "Tx", "M", "L", "W", "CH", "CI","IH","II","Ie","Ic","Deaths")
  	  colnames(LHSResults) <- c(paste(cnames,".1",sep=""),paste(cnames,".2",sep=""),
  	                            anames,
  	                            paste(inames,".1",sep=""),paste(inames,".2",sep=""),
  	                            paste(rnames,".1",sep=""),paste(rnames,".2",sep=""),
  	                            paste(dnames,".1",sep=""),paste(dnames,".2",sep=""))
  	  }
  	
  	for(i in 1:n) {
  	  cat(paste(i," of ",n,"\n"))
  	  aS     <- c(aS.1[i],aS.2[i])
  	  aXR    <- c(aXR.1[i],aXR.2[i])
  	  aC     <- c(aC.1[i],aC.2[i])
  	  aI     <- c(aI.1[i],aI.2[i])
  	  
  	  Cnp    <- c(Cnp.1[i],Cnp.2[i])   #Contact rate from nurses to patients in ward i
  	  Cdp    <- c(Cdp.1[i],Cdp.2[i])   #Contact rate from doctors to patients in ward i
  	  rn     <- c(rn.1[i],rn.2[i])     # Ratio of patients to nurses in each ward
  	  rd     <- c(rd.1[i],rd.2[i])     # Ratio of patients to doctors in each ward
  	  
  	  patients <- c(patients.1[i], patients.2[i])
  	  
  	  stay   <- c(stay.1[i], stay.2[i])
  	  trans  <- c(trans.1[i], trans.2[i])
  	  tau     <- c(tau.1[i], tau.2[i])
  	  
  	  # IntialParams puts all of these values into a data frame
  	  InitialParams <- as.data.frame(cbind(aS, aXR, aC, aI, Cnp, Cdp, rn, rd,
  	                                       patients, stay, trans, tau))
  	  
  	  omega   <- omega.lhs[i]
  	  rhoO    <- rhoO.lhs[i]
  	  pi      <- pi.lhs[i]
  	  phi     <- phi.lhs[i]
  	  rhoR    <- rhoR.lhs[i]
  	  epsilon <- epsilon.lhs[i]
  	  psi     <- psi.lhs[i]
  	  theta   <- theta.lhs[i]
  
  	  zeta    <- zeta.lhs[i]
  	  sigmaD  <- sigmaD.lhs[i]
  	  sigmaN  <- sigmaN.lhs[i]
  	  delta   <- delta.lhs[i]
  	  Qpn     <- Qpn.lhs[i]
  	  Qpd     <- Qpd.lhs[i]
  	  Qnp     <- Qnp.lhs[i]
  	  Qdp     <- Qdp.lhs[i]
  	  
  	  xi      <- xi.lhs[i]
  	  eta     <- eta.lhs[i]
  	  mu      <- mu.lhs[i]
  	  Pc      <- Pc.lhs[i]
  	  Psc     <- Psc.lhs[i]
  	  Pnc     <- Pnc.lhs[i]
  	  sigY    <- sigY.lhs[i]
  	  sigT    <- sigT.lhs[i]
  	  sigMc   <- sigMc.lhs[i]
  	  sigMsc  <- sigMsc.lhs[i]
  	  sigMnc  <- sigMnc.lhs[i]
  	  sigL    <- sigL.lhs[i]
  	  sigW    <- sigW.lhs[i]
  	  
  	  betaM   <- betaM.lhs[i]
  	  betaL   <- betaL.lhs[i]
  	  betaY   <- betaY.lhs[i]
  	  betaT   <- betaT.lhs[i]
  	  betaW   <- betaW.lhs[i]
  	  
  	  IC      <- IC.lhs[i]
  	  
  	  WardIndParams <- as.data.frame(cbind(omega, rhoO, pi, phi, rhoR, epsilon, psi, theta,  
  	                                       zeta, sigmaD, sigmaN, delta, Qpn, Qpd, Qnp, Qdp, xi, 
  	                                       eta, mu, Pc, Psc, Pnc, sigY, sigT, sigMc, sigMsc, sigMnc, 
  	                                       sigL, sigW, betaM, betaL, betaY, betaT, betaW, IC))   
    
  	  # Interventions
  	  if('sigD.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])){
  			sigD <- c(sigD.1.int[i], sigD.2.int[i]) 
  		} else {
  			sigD <- c(sigD.1[i], sigD.2[i]) 
  		}
  	
  		if('sigCP.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  			sigCP   <- c(sigCP.1.int[i], sigCP.2.int[i])
  		} else {
  			sigCP   <- c(sigCP.1[i], sigCP.2[i])
  		}
  	  
  	  if('PS.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    PS   <- c(PS.1.int[i], PS.2.int[i])
  	  } else if('PS.ICU.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])){
  	    PS   <- c(PS.1[i], PS.2.ICU[i])
  	  } else {
  	    PS   <- c(PS.1[i], PS.2[i])
  	  }
  	  
  	  if('OM.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    OM   <- c(OM.1.int[i], OM.2.int[i])
  	  } else if('OM.ICU.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])){
  	    OM   <- c(OM.1[i], OM.2.ICU[i])
  	  } else {
  	    OM   <- c(OM.1[i], OM.2[i])
  	  }
  	  
  	  
  	  if('sigN1.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigN1   <- c(sigN1.1.int[i], sigN1.2.int[i])
  	  } else {
  	    sigN1   <- c(sigN1.1[i], sigN1.2[i])
  	  }
  	  
  	  
  	  if('sigD1.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigD1   <- c(sigD1.1.int[i], sigD1.2.int[i])
  	  } else {
  	    sigD1   <- c(sigD1.1[i], sigD1.2[i])
  	  }
  	  
  	  if('sigAe.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigAe   <- c(sigAe.1.int[i], sigAe.2.int[i])
  	  } else {
  	    sigAe   <- c(sigAe.1[i], sigAe.2[i])
  	  }
  	  
  	  if('sigAm.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigAm   <- c(sigAm.1.int[i], sigAm.2.int[i])
  	  } else {
  	    sigAm   <- c(sigAm.1[i], sigAm.2[i])
  	  }
  	  
  	  if('sigAp.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigAp   <- c(sigAp.1.int[i], sigAp.2.int[i])
  	  } else {
  	    sigAp   <- c(sigAp.1[i], sigAp.2[i])
  	  }
  	  
  	  if('sigY1.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigY1   <- c(sigY1.1.int[i], sigY1.2.int[i])
  	  } else {
  	    sigY1   <- c(sigY1.1[i], sigY1.2[i])
  	  }
  	  
  	  if('sigT1.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigT1   <- c(sigT1.1.int[i], sigT1.2.int[i])
  	  } else {
  	    sigT1   <- c(sigT1.1[i], sigT1.2[i])
  	  }
  
  	  if('sigL1.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigL1   <- c(sigL1.1.int[i], sigL1.2.int[i])
  	  } else {
  	    sigL1   <- c(sigL1.1[i], sigL1.2[i])
  	  }
  	  
  	  if('sigW1.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigW1   <- c(sigW1.1.int[i], sigW1.2.int[i])
  	  } else {
  	    sigW1   <- c(sigW1.1[i], sigW1.2[i])
  	  }
  	  
  	  if('sigMc1.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigMc1   <- c(sigMc1.1.int[i], sigMc1.2.int[i])
  	  } else {
  	    sigMc1   <- c(sigMc1.1[i], sigMc1.2[i])
  	  }
  	  
  	  if('sigMsc1.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigMsc1   <- c(sigMsc1.1.int[i], sigMsc1.2.int[i])
  	  } else {
  	    sigMsc1   <- c(sigMsc1.1[i], sigMsc1.2[i])
  	  }
  	  
  	  if('sigMnc1.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    sigMnc1   <- c(sigMnc1.1.int[i], sigMnc1.2.int[i])
  	  } else {
  	    sigMnc1   <- c(sigMnc1.1[i], sigMnc1.2[i])
  	  }
  	  
  	  if('f.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    f   <- c(f.1.int[i], f.2.int[i])
  	  } else {
  	    f   <- c(f.1[i], f.2[i])
  	  }
  	  
  	  if('Fn.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
  	    Fn   <- c(Fn.1.int[i], Fn.2.int[i])
  	  } else {
  	    Fn   <- c(Fn.1[i], Fn.2[i])
  	  }
  
  	  InterventionParams <- as.data.frame(cbind(sigD, PS, OM, sigCP, sigN1, sigD1, sigY1, sigT1, sigL1,
  	                                            sigW1, sigMc1, sigMsc1, sigMnc1, sigAe, sigAm, sigAp, f,
  	                                            Fn))
  	  
  	  # Params breaks up the ward-dependent parameters by ward into Params1 and
  	  # Params2, and adds the ward independent parameters at the end of the list.
  	  # The final params is a list with Params1 and 2 as
  	  # named numeric values.  That is what the names(Params..etc) is doing, naming
  	  # the values so they can be called later.
  	  Params <- list(Params1=as.numeric(cbind(InitialParams[1,], WardIndParams, InterventionParams[1,])),
  	                 Params2=as.numeric(cbind(InitialParams[2,], WardIndParams, InterventionParams[2,])))
  	  names(Params$Params1) <- c(colnames(InitialParams), colnames(WardIndParams), colnames(InterventionParams))
  	  names(Params$Params2) <- c(colnames(InitialParams), colnames(WardIndParams), colnames(InterventionParams))
  	  
  	  ## Running the ODE model at baselin to get to baseline equilibrium
  	  for(w in 1:nrow(InitialState)){
  	    if(w==1){
  	      y = as.numeric(as.vector(InitialState[w,]))
  	    } else {
  	      y = c(y,as.numeric(as.vector(InitialState[w,])))
  	    }
  	  }
  	  output1 <- lsoda(y = y, times = seq(1,70,by=1), func = CREmodel, parms = Params, DTVals = InitialState)
  	  
  	  # Run the ODE model of the given IPC for the defined number of days (if statement applies to cohorting)
  	  
  	  if(intnum %in% 19:21){
  	    EquilibriumState <- rbind.data.frame(c(output1[70, 2:5], 0, output1[70, 6], 0, output1[70, 7:18]), 
  	                                         c(output1[70, 19:22], 0, output1[70, 23], 0, output1[70, 24:35]))
  	    colnames(EquilibriumState) = c("S", "XR", "Cu", "Cd1", "Cd2", "I1", "I2", "Y", "Tx", "M", "L", "W",
  	                                   "CH", "CI", "IH", "II", "Ie", "Ic", "Deaths")
  	    
  	    for(w in 1:nrow(EquilibriumState)){
  	      if(w==1){
  	        z = as.numeric(as.vector(EquilibriumState[w,]))
  	      } else {
  	        z = c(z,as.numeric(as.vector(EquilibriumState[w,])))
  	      }
  	    }
  	    output <- lsoda(y = z, times = seq(1,days,by=1), func = CREmodel_Cohorting, parms = Params, DTVals = EquilibriumState)
  	    
  	    cat((output[days,16]*numbeds[i] + output[days,35]*numbeds[i]),"\n")
  	    
  	    LHSResults[i,] = c(unlist(InitialParams[1,]),unlist(InitialParams[2,]),unlist(WardIndParams),unlist(InterventionParams[1,]),
  	                       unlist(InterventionParams[2,]),unlist(output[length(output[,1]),2:39]),(sum(c(sum(output[,7]),sum(output[,8])))),(sum(c(sum(output[,26]),sum(output[,27])))))
  	  } else {
    	  EquilibriumState <- rbind.data.frame(output1[70, 2:18], output1[70, 19:35])
    	  colnames(EquilibriumState) = c("S", "XR", "Cu", "Cd", "I", "Y", "Tx", "M", "L", "W", "CH", "CI",
    	                                 "IH", "II", "Ie", "Ic", "Deaths")
    	  
    	  for(w in 1:nrow(EquilibriumState)){
    	    if(w==1){
    	      z = as.numeric(as.vector(EquilibriumState[w,]))
    	    } else {
    	      z = c(z,as.numeric(as.vector(EquilibriumState[w,])))
    	    }
    	  }
    	  output <- lsoda(y = z, times = seq(1,days,by=1), func = CREmodel, parms = Params, DTVals = EquilibriumState)
    	  
    	  cat((output[days,14]*numbeds[i] + output[days,31]*numbeds[i]),"\n")
    	  
    	  LHSResults[i,] = c(unlist(InitialParams[1,]),unlist(InitialParams[2,]),unlist(WardIndParams),unlist(InterventionParams[1,]),
    	                     unlist(InterventionParams[2,]),unlist(output[length(output[,1]),2:35]),(sum(output[,6])),(sum(output[,23])))
    	  
    	  }

  	  }
  	
  	 ## need to add calls for Cohorting here to add to total LHS Results
  	
  	 ## Getting baseline deaths and infections to calculate IPC infections and deaths averted
  	 if(colnames(ipv[intnum])=="BSLN"){
  	   LHSInfDeathBSLN <- cbind(LHSResults[,c("IH.1","IH.2","Deaths.1","Deaths.2")],
  	                            InfTot = rowSums(LHSResults[,c("IH.1","IH.2")])*sum(patients),
  	                            DeathTot = rowSums(LHSResults[,c("Deaths.1","Deaths.2")])*sum(patients),
  	                            InfperMil = rowSums(LHSResults[,c("IH.1","IH.2")])*(1000000/365),
  	                            DthsperMil = rowSums(LHSResults[,c("Deaths.1","Deaths.2")])*(1000000/365))
  	   }  
     
  	write.csv(LHSResults, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","LHSResults","_",format((round(100 - mean(aS.1)*100, digits=0))),"Colonized","_",paste(hpv[1,envnum]),"_",paste(colnames(ipv[intnum])), ".csv", sep=""))
    
     
    # Conducting the PRCC and saving the outputs
    if(intnum  <  19){
    	ICR = rowSums(LHSResults[,c("I.1","I.2")])
    	x.sd<-apply(LHSResults,2,sd)
    	testprcc = cbind(LHSResults[,which(x.sd[1:91]!=0&names(x.sd[1:91])!="aS.1"&names(x.sd[1:91])!="X"&names(x.sd[1:91])!="aS.2")],ICR)
    	I_PRCC <- PRCC(testprcc,sort.results=T,sort.abs = T)
    	print(I_PRCC)
    	write.csv(I_PRCC, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","PRCCValues","/","PRCC_Inf","_",paste(hpv[1,envnum]),"_",paste(colnames(ipv[intnum])), ".csv", sep=""))
    	
    	ICR = rowSums(LHSResults[,c("Cu.1","Cd.1","Cu.2","Cd.2")])
    	x.sd<-apply(LHSResults,2,sd)
    	testprcc = cbind(LHSResults[,which(x.sd[1:91]!=0&names(x.sd[1:91])!="aS.1"&names(x.sd[1:91])!="X"&names(x.sd[1:91])!="aS.2")],ICR)
    	C_PRCC <- PRCC(testprcc,sort.results=T,sort.abs = T)
    	print(C_PRCC)
    	write.csv(C_PRCC, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","PRCCValues","/","PRCC_Col","_",paste(hpv[1,envnum]),"_",paste(colnames(ipv[intnum])), ".csv", sep=""))
    	
    	ICR = rowSums(LHSResults[,c("Cu.1","Cd.1","Cu.2","Cd.2","I.1","I.2")])
    	x.sd<-apply(LHSResults,2,sd)
    	testprcc = cbind(LHSResults[,which(x.sd[1:91]!=0&names(x.sd[1:91])!="aS.1"&names(x.sd[1:91])!="X"&names(x.sd[1:91])!="aS.2")],ICR)
    	IC_PRCC <- PRCC(testprcc,sort.results=T,sort.abs = T)
    	print(IC_PRCC)
    	write.csv(IC_PRCC, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","PRCCValues","/","PRCC_InfCol","_",paste(hpv[1,envnum]),"_",paste(colnames(ipv[intnum])), ".csv", sep=""))
    } else {
      ICR = rowSums(LHSResults[,c("I1.1","I1.2","I2.1","I2.2")])
      x.sd<-apply(LHSResults,2,sd)
      testprcc = cbind(LHSResults[,which(x.sd[1:91]!=0&names(x.sd[1:91])!="aS.1"&names(x.sd[1:91])!="X"&names(x.sd[1:91])!="aS.2")],ICR)
      I_PRCC <- PRCC(testprcc,sort.results=T,sort.abs = T)
      print(I_PRCC)
      write.csv(I_PRCC, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","PRCCValues","/","PRCC_Inf","_",paste(hpv[1,envnum]),"_",paste(colnames(ipv[intnum])), ".csv", sep=""))
  
      ICR = rowSums(LHSResults[,c("Cu.1","Cd1.1","Cd2.1","Cu.2","Cd1.2","Cd2.2")])
      x.sd<-apply(LHSResults,2,sd)
      testprcc = cbind(LHSResults[,which(x.sd[1:91]!=0&names(x.sd[1:91])!="aS.1"&names(x.sd[1:91])!="X"&names(x.sd[1:91])!="aS.2")],ICR)
      C_PRCC <- PRCC(testprcc,sort.results=T,sort.abs = T)
      print(C_PRCC)
      write.csv(C_PRCC, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","PRCCValues","/","PRCC_Col","_",paste(hpv[1,envnum]),"_",paste(colnames(ipv[intnum])), ".csv", sep=""))
  
      ICR = rowSums(LHSResults[,c("Cu.1","Cd1.1","Cd2.1","Cu.2","Cd1.2","Cd2.2","I1.1","I1.2","I2.1","I2.2")])
      x.sd<-apply(LHSResults,2,sd)
      testprcc = cbind(LHSResults[,which(x.sd[1:91]!=0&names(x.sd[1:91])!="aS.1"&names(x.sd[1:91])!="X"&names(x.sd[1:91])!="aS.2")],ICR)
      IC_PRCC <- PRCC(testprcc,sort.results=T,sort.abs = T)
      print(IC_PRCC)
      write.csv(IC_PRCC, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","PRCCValues","/","PRCC_InfCol","_",paste(hpv[1,envnum]),"_",paste(colnames(ipv[intnum])), ".csv", sep=""))
    }
      
    ## Create a list of values for plots
    if(intnum == 9){
      InfAverted  <- cbind(LHSInfDeathBSLN[, "InfperMil"] - (rowSums(LHSResults[,c("IH.1","IH.2")])*(1000000/365)))
      InfAvByWrd  <- cbind((LHSInfDeathBSLN[,"IH.1"] - LHSResults[,"IH.1"])*(1000000/365),
                           (LHSInfDeathBSLN[,"IH.2"] - LHSResults[,"IH.2"])*(1000000/365))
      DthsAverted <- cbind(LHSInfDeathBSLN[,"DthsperMil"] - (rowSums(LHSResults[,c("Deaths.1","Deaths.2")])*(1000000/365)))
      
    } else if(intnum > 9){
      InfAverted  <- cbind(InfAverted,  LHSInfDeathBSLN[, "InfperMil"] - (rowSums(LHSResults[,c("IH.1","IH.2")])*(1000000/365)))
      InfAvByWrd  <- cbind(InfAvByWrd, (LHSInfDeathBSLN[,"IH.1"] - LHSResults[,"IH.1"])*(1000000/365),
                           (LHSInfDeathBSLN[,"IH.2"] - LHSResults[,"IH.2"])*(1000000/365))
      DthsAverted <- cbind(DthsAverted, LHSInfDeathBSLN[,"DthsperMil"] - (rowSums(LHSResults[,c("Deaths.1","Deaths.2")])*(1000000/365)))

    }
    
  }
  
  ## Plots for the particular HCF and environment are done through the following
  
  colnames(InfAverted)  <- colnames(ipv[9:21])
  colnames(DthsAverted) <- colnames(ipv[9:21])
  
  colnames(InfAvByWrd)   <- paste(rep(colnames(ipv[9:21]), each = 2), c(1,2), sep = ".")
  
  
  ## Exporting the raw data for the net cost and ICER Tables
  write.csv(LHSResults, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","LHSResults","_",format((round(100 - mean(aS.1)*100, digits=0))),"Colonized","_",paste(hpv[1,envnum]),"_",paste(colnames(ipv[intnum])), ".csv", sep=""))
  
  ## Plot of infections averted and deaths averted for individual IPCs
  
  colors <- brewer.pal(n = 8, "Paired")
  
  IPCnames <- c("Contact Precautions", "Chlorhexidine\n(CHG) Bathing", "Hand Hygiene", "Judicious Use",
                "Environmental Cleaning of\nNon-Water Sources", "Environmental Cleaning of\nWater Sources",
                "Cohorting")
  
  PlotDataIPCs <- data.frame(IPCnames, "Infections Averted" = unlist(colMeans(InfAverted[,c(1:2,5:6,9:11)])), "Infections Averted SE" = unlist(apply(InfAverted[,c(1:2,5:6,9:11)], 2, sd)/sqrt(n)), 
                         "Deaths Averted" = unlist(colMeans(DthsAverted[,c(1:2,5:6,9:11)])), "Deaths Averted SE" = unlist(apply(DthsAverted[,c(1:2,5:6,9:11)], 2, sd)/sqrt(n)))
  
  PlotDataIPCs$IPCnames <- factor(PlotDataIPCs$IPCnames, levels = PlotDataIPCs$IPCnames[order(PlotDataIPCs$Infections.Averted)])
  
  ggplot(PlotDataIPCs, aes(x = IPCnames, y = Infections.Averted)) +
    geom_bar(stat = "identity", color="black", position=position_dodge(), fill = colors[1]) +
    theme_classic() +
    labs(y = "Infections Averted/1,000,000 Admissions") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = 1, size = 8.5)) +
    geom_errorbar(aes(ymin = Infections.Averted - Infections.Averted.SE,
                      ymax = Infections.Averted + Infections.Averted.SE), width = 0.5)
  
  ggsave(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),"Plots","Individual IPCs Infection Averted.eps",sep="/"), width = 11, height = 8.5, units = "in")

  PlotDataIPCs$IPCnames <- factor(PlotDataIPCs$IPCnames, levels = PlotDataIPCs$IPCnames[order(PlotDataIPCs$Deaths.Averted)])
  
  ggplot(PlotDataIPCs, aes(x = IPCnames, y = Deaths.Averted)) +
    geom_bar(stat = "identity", color="black", position=position_dodge(), fill = colors[3]) +
    theme_classic() +
    labs(y = "Deaths Averted/1,000,000 Admissions") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = 1, size = 8.5)) +
    geom_errorbar(aes(ymin = Deaths.Averted - Deaths.Averted.SE,
                      ymax = Deaths.Averted + Deaths.Averted.SE), width = 0.5)
  
  ggsave(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),"Plots","Individual IPCs Deaths Averted.eps",sep="/"), width = 11, height = 8.5, units = "in")
  
  ## Plot of Bundled IPC Interventions
  
  BNDLnames <- c("Active Surveillance,\nContact Precautions\n& CHG Bathing","ICU Only Active\nSurveillance,\nContact\nPrecautions\n& CHG Bathing", 
                 "Active Surveillance\n& Judicious Use", "ICU Only Active\nSurveillance\n& Judicious Use", "Cohorting &\nActive\nSurveillance",
                 "All IPCs")
  
  BNDLPlot <- data.frame(BNDLnames, "Infections Averted" = unlist(colMeans(InfAverted[,c(3:4,7:8,12:13)])), "Infections Averted SE" = unlist(apply(InfAverted[,c(3:4,7:8,12:13)], 2, sd)/sqrt(n)), 
                         "Deaths Averted" = unlist(colMeans(DthsAverted[,c(3:4,7:8,12:13)])), "Deaths Averted SE" = unlist(apply(DthsAverted[,c(3:4,7:8,12:13)], 2, sd)/sqrt(n)))
  
  BNDLPlot$BNDLnames <- factor(BNDLPlot$BNDLnames, levels = BNDLPlot$BNDLnames[order(BNDLPlot$Infections.Averted)])
  
  ggplot(BNDLPlot, aes(x = BNDLnames, y = Infections.Averted)) +
    geom_bar(stat = "identity", color="black", position=position_dodge(), fill = colors[2]) +
    theme_classic() +
    labs(y = "Infections Averted/1,000,000 Admissions") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = 1, size = 8.5)) +
    geom_errorbar(aes(ymin = Infections.Averted - Infections.Averted.SE,
                      ymax = Infections.Averted + Infections.Averted.SE), width = 0.5)
  
  ggsave(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),"Plots","Bundeled IPCs Infection Averted.eps",sep="/"), width = 9, height = 9, units = "in")
  
  BNDLPlot$BNDLnames <- factor(BNDLPlot$BNDLnames, levels = BNDLPlot$BNDLnames[order(BNDLPlot$Deaths.Averted)])
  
  ggplot(BNDLPlot, aes(x = BNDLnames, y = Deaths.Averted)) +
    geom_bar(stat = "identity", color="black", position=position_dodge(), fill = colors[4]) +
    theme_classic() +
    labs(y = "Deaths Averted/1,000,000 Admissions") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = 1, size = 8.5)) +
    geom_errorbar(aes(ymin = Deaths.Averted - Deaths.Averted.SE,
                      ymax = Deaths.Averted + Deaths.Averted.SE), width = 0.5)
  
  ggsave(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),"Plots","Bundeled IPCs Deaths Averted.eps",sep="/"), width = 8, height = 9, units = "in")
  
  
  ## Incoporating costing into the simulated data

  cpv <- CSTparamvals[which(CSTparamvals[paste(hpv[1,envnum])]==1),1:4]
  
  ## Defining the costs of each intervention by the environmental setting
  CP.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="CP.cost"),2:4])
  CHG.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="CHG.cost"),2:4])
  ASCPCHG.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="ASCPCHG.cost"),2:4])
  ICUASCPCHG.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="ICUASCPCHG.cost"),2:4])
  HH.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="HH.cost"),2:4])

  # We did not cost for the enhanced cleaning of water sources in this model
  
  EC4.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="EC4.cost"),2:4])
  JU.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="JU.cost"),2:4])
  ASJU.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="ASJU.cost"),2:4])
  ICUASJU.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="ICUASJU.cost"),2:4])
  ASCHRT.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="ASCHRT.cost"),2:4])
  CHRT.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="CHRT.cost"),2:4])
  ALL.cost = triangleLHS_input(n, cpv[which(cpv$Parameter=="ALL.cost"),2:4])
  
  ## Defining the cost per infection averted (the GW and ICU have seperate costs in India)
  
  percost.gw  = triangleLHS_input(n, cpv[which(cpv$Parameter=="CostPerInf.GW"),2:4])
  percost.icu = triangleLHS_input(n, cpv[which(cpv$Parameter=="CostPerInf.ICU"),2:4])
  
  ## Creating a table to display the net costs of each IPC
  ## 'Enhanced Cleaning of Water Sources' is omitted
  
  costtable <- data.frame(rep(NA, n))
  for(icn in c(9:17,19:21)){
      costtemp <- data.frame(InfAvByWrd[,paste(colnames(ipv[,icn]),c(1,2),sep = ".")], 
                              percost.gw, percost.icu, 
                              DthsAverted[,colnames(ipv[,icn])],
                              mget(paste(colnames(ipv[,icn]),"cost",sep = ".")))
      # Calculating the Net Cost
      costtemp$netcost <- (costtemp[,6] - ((costtemp[,1]*costtemp[,3]) + (costtemp[,2]*costtemp[,4])))/costtemp[,5]
      colnames(costtemp)[c(5,7)] <- c(paste(colnames(ipv[,icn]),"Deaths", sep = "."), paste(colnames(ipv[,icn]),"Net", sep = "."))
      costtable <- cbind(costtable, costtemp)
  }
  costtable <- costtable[,-1]
  
  write.csv(costtable, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","CostTables","/","Net","_","Costs","_",paste(hpv[1,envnum]), ".csv", sep=""))
  
  currency <- cbind(c("HET_IN", "LET_IN", "HET_SA", "LET_SA"), c("Indian Rupees (IRS)", "Indian Rupees (IRS)", "South Africa Rand (ZAR)", "South Africa Rand (ZAR)"))
  
  NetCostNames <- c("Contact Precautions", "Chlorhexidine\n(CHG) Bathing", "Active Surveillance,\nContact Precautions\n& CHG Bathing",
                "ICU Only Active\nSurveillance,\nContact\nPrecautions\n& CHG Bathing", "Hand Hygiene", "Judicious Use",
                "Active Surveillance\n& Judicious Use", "ICU Only Active\nSurveillance\n& Judicious Use",
                "Environmental\nCleaning of\nNon-Water\nSources", "Cohorting","Cohorting &\nActive\nSurveillance","All IPCs")
  
  NetCostPlot <- data.frame(NetCostNames, "Net Cost" = unlist(colMeans(costtable[,which(colnames(costtable) %in% paste(colnames(ipv[,c(9:17,19:21)]),"Net", sep = "."))])), 
                            "Net Cost SE" = unlist(apply(costtable[,which(colnames(costtable) %in% paste(colnames(ipv[,c(9:17,19:21)]),"Net", sep = "."))], 2, sd))/sqrt(n))

  NetCostPlot$NetCostNames <- factor(NetCostPlot$NetCostNames, levels = NetCostPlot$NetCostNames[order(1:12)])

  ggplot(NetCostPlot, aes(x = NetCostNames, y = Net.Cost)) +
    geom_bar(stat = "identity", color="black", position=position_dodge(), fill = colors[5]) +
    theme_classic() +
    labs(y = paste("Net", "Cost", "in", currency[which(currency[,1] == paste0(hpv[1,envnum])),2], sep = " ")) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, vjust = 1, size = 8.5),
          axis.text.y = element_text(size = 9)) +
    geom_errorbar(aes(ymin = Net.Cost - Net.Cost.SE,
                      ymax = Net.Cost + Net.Cost.SE), width = 0.5) +
    scale_y_continuous(trans = "pseudo_log", breaks = c(-10^6, -10^5, -10^4, -10^3, -10^2, -10, 0, 10, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8),
                       labels = c(expression(-10^6), expression(-10^5), expression(-10^4), expression(-10^3), expression(-10^2), "-10", "0", "10", expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6), expression(10^7), expression(10^8)))
  
  ggsave(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),"Plots","Net Cost.eps",sep="/"), width = 14, height = 9, units = "in")
  
  ICERTable <- data.frame(
    "ASCPCHGtCHG"    = as.numeric(unlist(((costtable["ASCPCHG.cost"] - ((costtable["ASCPCHG.1"]*percost.gw) + (costtable["ASCPCHG.2"]*percost.icu)))          - (costtable["CHG.cost"] - ((costtable["CHG.1"]*percost.gw) + (costtable["CHG.2"]*percost.icu))))             / (costtable["ASCPCHG.Deaths"] - costtable["CHG.Deaths"]))),
    "ICUASCPCHGtCHG" = as.numeric(unlist(((costtable["ICUASCPCHG.cost"] - ((costtable["ICUASCPCHG.1"]*percost.gw) + (costtable["ICUASCPCHG.2"]*percost.icu))) - (costtable["CHG.cost"] - ((costtable["CHG.1"]*percost.gw) + (costtable["CHG.2"]*percost.icu))))             / (costtable["ICUASCPCHG.Deaths"] - costtable["CHG.Deaths"]))),
    "AStJU"          = as.numeric(unlist(((costtable["ASJU.cost"] - ((costtable["ASJU.1"]*percost.gw) + (costtable["ASJU.2"]*percost.icu)))                   - (costtable["JU.cost"] - ((costtable["JU.1"]*percost.gw) + (costtable["JU.2"]*percost.icu))))                / (costtable["ASJU.Deaths"] - costtable["JU.Deaths"]))),
    "ICUAStJU"       = as.numeric(unlist(((costtable["ICUASJU.cost"] - ((costtable["ICUASJU.1"]*percost.gw) + (costtable["ICUASJU.2"]*percost.icu)))          - (costtable["JU.cost"] - ((costtable["JU.1"]*percost.gw) + (costtable["JU.2"]*percost.icu))))                / (costtable["ICUASJU.Deaths"] - costtable["JU.Deaths"]))),
    "AStCHRT"        = as.numeric(unlist(((costtable["ASCHRT.cost"] - ((costtable["ASCHRT.1"]*percost.gw) + (costtable["ASCHRT.2"]*percost.icu)))             - (costtable["CHRT.cost"] - ((costtable["CHRT.1"]*percost.gw) + (costtable["CHRT.2"]*percost.icu))))          / (costtable["ASCHRT.Deaths"] - costtable["CHRT.Deaths"]))),
    "ALLtASCPCHG"    = as.numeric(unlist(((costtable["ALL.cost"] - ((costtable["ALL.1"]*percost.gw) + (costtable["ALL.2"]*percost.icu)))                      - (costtable["ASCPCHG.cost"] - ((costtable["ASCPCHG.1"]*percost.gw) + (costtable["ASCPCHG.2"]*percost.icu)))) / (costtable["ALL.Deaths"] - costtable["ASCPCHG.Deaths"])))
  )
  
  ICERTableNames <- c("Active Surveillance,\nContact Precautions\n& CHG Bathing\u2020", "ICU Only Active\nSurveillance,\nContact Precautions\n& CHG Bathing\u2020", 
                      "Active Surveillance\n& Judicious Use\u2021", "ICU Only Active\nSurveillance\n& Judicious Use\u2021", "Cohorting &\nActive Surveillance\u00A7", "All IPCs\u2016")
  
  ICERTablePlot <- data.frame(ICERTableNames, "Incremental Cost-Effectiveness Ratio" = unlist(colMeans(ICERTable)), 
                              "Incremental Cost-Effectiveness Ratio SE" = unlist(apply(ICERTable, 2, sd))/sqrt(n))
  
  ICERTablePlot$ICERTableNames <- factor(ICERTablePlot$ICERTableNames, levels = ICERTablePlot$ICERTableNames[order(1:6)])
  
  
  ggplot(ICERTablePlot, aes(x = ICERTablePlot$ICERTableNames, y = Incremental.Cost.Effectiveness.Ratio)) +
    geom_bar(stat = "identity", color="black", position=position_dodge(), fill = colors[7]) +
    theme_classic() +
    labs(y = paste("ICER", "Cost", "in", currency[which(currency[,1] == paste0(hpv[1,envnum])),2], sep = " ")) +
    geom_errorbar(aes(ymin = Incremental.Cost.Effectiveness.Ratio - Incremental.Cost.Effectiveness.Ratio.SE,
                      ymax = Incremental.Cost.Effectiveness.Ratio + Incremental.Cost.Effectiveness.Ratio.SE), width = 0.5) +
    scale_y_continuous(trans = "log10", breaks = c(10, 100, 1000, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10),
                       labels = c(expression(10), expression(10^2), expression(10^3), expression(10^4), expression(10^5), 
                                  expression(10^6), expression(10^7), expression(10^8), expression(10^9), expression(10^10))) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8.5),
          axis.text.y = element_text(size = 9))
  
  ggsave(paste(paste("hCREiAM Data",runstartday,""),colnames(hpv[envnum]),"Plots","ICER.eps",sep="/"), width = 8, height = 9, units = "in")
  
  ## Exporting the raw data for the net cost and ICER Tables
  write.csv(ICERTable, file = paste(paste("hCREiAM Data",runstartday,""),"/",colnames(hpv[envnum]),"/","CostTables","/","ICER","_",paste(hpv[1,envnum]), ".csv", sep=""))
  
  
}

