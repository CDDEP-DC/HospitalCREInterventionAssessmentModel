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
##setwd()

rm(list=ls()) # clear
set.seed(1982)
options(stringsAsFactors=F)
if (!("plyr" %in% rownames(installed.packages()))){
  install.packages("plyr")
}
library(plyr)

if (!("readxl" %in% rownames(installed.packages()))){
  install.packages("readxl")
}
library(readxl)

if (!("deSolve" %in% rownames(installed.packages()))){
  install.packages("deSolve")
}
library(deSolve)

if (!("data.table" %in% rownames(installed.packages()))){
  install.packages("data.table")
}
library(data.table)

if (!("ggplot2" %in% rownames(installed.packages()))){
  install.packages("ggplot2")
}
library(ggplot2)

if (!("triangle" %in% rownames(installed.packages()))){
  install.packages("triangle")
}
library(triangle)

if (!("lhs" %in% rownames(installed.packages()))){
  install.packages("lhs")
}
library(lhs)


## function to convert tibble read in format to comma separated list and call triangleLHS
triangleLHS_input <- function(n, vals){
	v=as.numeric(unlist(vals))
	triangleLHS(n, v[1], v[2], v[3])    
}

triangleLHS <- function(n, min, max, mode = (min + max)/2){
    p <- randomLHS(n, 1)[,1]
    qtriangle(p = p, a = min, b = max, c = mode)
}

#Source the model code
source("Model.R")

## *** Comcol file needs to be updated to use information based on the values entered in the spreadsheet *** 
source("ComCol.R")

days <- 365  #Number of days simulation will run *** Need to add as a parameter 

# NUMBER OF WARDS - The first number in each parenthesis below is the
# value for ward 1, the second number will be the value for ward 2.  All
# parentheses MUST have the same number of values!
Cu <- c(0, 0)         # Proportion of patients colonized with resistant pathogen and it's undetected
Cd <- c(0, 0)         # Proportion of patients colonized with resistant pathogen and it's detected
I  <- c(0, 0)         # Proportion of patients infected with resistant pathogen
XR <- c(0, 0)         # Proportion of patients successfully prophylaxed against resistant pathogen

Y  <- c(0, 0)         # Environmental contamination on surfaces (excluding water)
Tx <- c(0, 0)         # Environmental contamination on textiles
M  <- c(0, 0)         # Environmental contamination on medical textiles
L  <- c(0, 0)         # Environmental contamination on mobile medial textiles
W  <- c(0, 0)         # Environmental contamination on water sources

CH <- c(0, 0)         # Number of acummulated colonized patients created inside the healthcare facility
CI <- c(0, 0)         # Number of accumalated colonized patients imported into the hospital

IH <- c(0, 0)         # Number of acummulated infected patients created inside the healthcare facility
II <- c(0, 0)         # Number of accumalated infected patients imported into the hospital

Ie <- c(0, 0)         # Number of acummulated infected patients through the environment
Ic <- c(0, 0)         # Number of accumalated infected patients through contacts

Deaths <- c(0, 0)     # Total number of deaths

# Intial state puts all of these values into a data frame
InitialState <- as.data.frame(cbind(S, XR, Cu, Cd, I, Y, Tx, M, L, W, CH, CI, IH, II, Ie, Ic, Deaths))

patients <- numpatients*rowSums(InitialState[1:5])  # Total number of patients in each ward

#######################################################
# LHS 
#######################################################

# Need this for all reps by group
aS.1  <- triangleLHS(n,minc,maxc,medc) 
aS.2  <- triangleLHS(n,minc,maxc,medc) 
aC.1  <- (1-aS.1)
aC.2  <- (1-aS.2)
aI.1  <- rep(0, n)
aI.2  <- rep(0, n)
aXR.1 <- rep(0, n)
aXR.2 <- rep(0, n)

# Read in the baseline parameters -- Move parameters spreadsheet name and associated parameters to the top to initialize (perhaps add to comcol)
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
IC.lhs = triangleLHS_input(n, bpv[which(bpv$Parameter=="IC.lhs"),2:4])

# Update the baseline parameters not set in the excel spreadsheet
Pc.lhs      = aP*bP      
Psc.lhs     = aP*(1 - bP)      
Pnc.lhs     = 1 - aP      

rn.1 = triangleLHS(n, minrn.1, maxrn.1, medrn.1)
rn.2 = triangleLHS(n, minrn.2, maxrn.2, medrn.2)
rd.1 = triangleLHS(n, minrd.1, maxrd.1, medrd.1)
rd.2 = triangleLHS(n, minrd.2, maxrd.2, medrd.2)

patients.1 = rep(numpatients*rowSums(InitialState[1,1:5]),n)  # Total number of patients in each ward
patients.2 = rep(numpatients*rowSums(InitialState[2,1:5]),n)

Qnp.lhs = triangleLHS(n, minqdp, maxqdp, medqdp)
Qdp.lhs = triangleLHS(n, minqdp, maxqdp, medqdp)

IC.lhs  = IC.lhs*10^ICOrder 

# Now read in the intervention parameters at baseline and for each intervention
interventionparamvals <- read_excel("hCREiAM_IPC_Parameters.xlsx",sheet="Interventions",skip=1)
ipv = interventionparamvals

#first the baselines
sigD.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD.1"),2:4])
sigD.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD.2"),2:4])
PS.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="PS.1"),2:4])
PS.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="PS.2"),2:4])
OM.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="OM.1"),2:4])
OM.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="OM.2"),2:4])
sigCP.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigCP.1"),2:4])
sigCP.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigCP.2"),2:4])
sigN1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigN1.1"),2:4])
sigN1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigN1.2"),2:4])
sigD1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD1.1"),2:4])
sigD1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD1.2"),2:4])
sigY1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigY1.1"),2:4])
sigY1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigY1.2"),2:4])
sigT1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigT1.1"),2:4])
sigT1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigT1.2"),2:4])
sigL1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigL1.1"),2:4])
sigL1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigL1.2"),2:4])
sigW1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigW1.1"),2:4])
sigW1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigW1.2"),2:4])
sigMc1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMc1.1"),2:4])
sigMc1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMc1.2"),2:4])
sigMsc1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMsc1.1"),2:4])
sigMsc1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMsc1.2"),2:4])
sigMnc1.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMnc1.1"),2:4])
sigMnc1.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigMnc1.2"),2:4])
sigAe.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAe.1"),2:4])
sigAe.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAe.2"),2:4])
sigAm.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAm.1"),2:4])
sigAm.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAm.2"),2:4])
sigAp.1 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAp.1"),2:4])
sigAp.2 = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigAp.2"),2:4])

#then the intervention vals
sigD.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD.1"),5:7])
sigD.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="sigD.2"),5:7])

# Active Surveillance
PS.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="PS.1"),5:7])
PS.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="PS.2"),5:7])
OM.1.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="OM.1"),5:7])
OM.2.int = triangleLHS_input(n, ipv[which(ipv$Parameter=="OM.2"),5:7])

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


for(intnum in 13:14) {
	LHSResults = matrix(NA,nrow=n,ncol=127)
	cnames = c("aS", "aXR", "aC", "aI", "Cnp", "Cdp", "rn", "rd", "patients", "stay", "trans", "tau")
	anames = c("omega", "rhoO", "pi", "phi", "rhoR", "epsilon", "psi", 'theta', "zeta",
	           "sigmaD", "sigmaN", "delta", "Qpn", "Qpd", "Qnp", "Qdp", "xi", "eta", "mu",
	           "Pc", "Psc", "Pnc", "sigY", "sigT", "sigMc", "sigMsc", "sigMnc", "sigL", "sigW",
	           "betaM", "betaL", "betaY", "betaT", "betaW", "IC")
	inames = c("sigD", "PS", "OM", "sigCP", "sigN1", "sigD1", "sigY1", "sigT1", "sigL1",
	           "sigW1", "sigMc1", "sigMsc1", "sigMnc1", "sigAe", "sigAm", "sigAp")
	rnames = c("S", "XR", "Cu", "Cd", "I", "Y", "Tx", "M", "L", "W", "CH", "CI","IH","II","Ie","Ic","Deaths")
	dnames = c("IDays")
	colnames(LHSResults) <- c(paste(cnames,".1",sep=""),paste(cnames,".2",sep=""),
	                          anames,
	                          paste(inames,".1",sep=""),paste(inames,".2",sep=""),
	                          paste(rnames,".1",sep=""),paste(rnames,".2",sep=""),
	                          paste(dnames,".1",sep=""),paste(dnames,".2",sep=""))
	
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
	  } else {
	    PS   <- c(PS.1[i], PS.2[i])
	  }
	  
	  if('OM.1' %in% unlist(ipv[which(ipv[,intnum]==1),1])) {
	    OM   <- c(OM.1.int[i], OM.2.int[i])
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
	  

	  InterventionParams <- as.data.frame(cbind(sigD, PS, OM, sigCP, sigN1, sigD1, sigY1, sigT1, sigL1,
	                                            sigW1, sigMc1, sigMsc1, sigMnc1, sigAe, sigAm, sigAp))
	  
	  # Params breaks up the ward-dependent parameters by ward into Params1 and
	  # Params2, and adds the ward independent parameters at the end of the list.
	  # The final params is a list with Params1 and 2 as
	  # named numeric values.  That is what the names(Params..etc) is doing, naming
	  # the values so they can be called later.
	  Params <- list(Params1=as.numeric(cbind(InitialParams[1,], WardIndParams, InterventionParams[1,])),
	                 Params2=as.numeric(cbind(InitialParams[2,], WardIndParams, InterventionParams[2,])))
	  names(Params$Params1) <- c(colnames(InitialParams), colnames(WardIndParams), colnames(InterventionParams))
	  names(Params$Params2) <- c(colnames(InitialParams), colnames(WardIndParams), colnames(InterventionParams))
	  
	  for(w in 1:nrow(InitialState)){
	    if(w==1){
	      y = as.numeric(as.vector(InitialState[w,]))
	    } else {
	      y = c(y,as.numeric(as.vector(InitialState[w,])))
	    }
	  }
	  output1 <- lsoda(y = y, times = seq(1,70,by=1), func = CREmodel, parms = Params, DTVals = InitialState)
	  
	  EquilibriumState <- rbind.data.frame(output1[70, 2:18], output1[70, 19:35])
	  colnames(EquilibriumState) = c("S", "XR", "Cu", "Cd", "I", "Y", "Tx", "M", "L", "W", "CH", "CI", "IH", "II", "Ie", "Ic", "Deaths")
	  
	  for(w in 1:nrow(EquilibriumState)){
	    if(w==1){
	      z = as.numeric(as.vector(EquilibriumState[w,]))
	    } else {
	      z = c(z,as.numeric(as.vector(EquilibriumState[w,])))
	    }
	  }
	  output <- lsoda(y = z, times = seq(1,365,by=1), func = CREmodel, parms = Params, DTVals = EquilibriumState)
	  
	  cat((output[365,14]*392 + output[365,31]*392),"\n")
	  
	  LHSResults[i,] = c(unlist(InitialParams[1,]),unlist(InitialParams[2,]),unlist(WardIndParams),unlist(InterventionParams[1,]),unlist(InterventionParams[2,]),unlist(output[length(output[,1]),2:35]),(sum(output[,6])),(sum(output[,23])))
	  }
	
	  ## need to add calls for Cohorting here to add to total LHS Results
	  
   write.csv(LHSResults, paste("LHSResults","_",format((round(100 - mean(aS.1)*100, digits=0))),"Colonized","_",paste(enviro[1]),"_",paste(colnames(ipv[intnum])),"_",format(Sys.time(), "%Y%m%d_%H%M"), ".csv", sep=""))
}


### Need to add call to PRCC file to output which should then output to a folder
### add ggplot print calls as well to print files
