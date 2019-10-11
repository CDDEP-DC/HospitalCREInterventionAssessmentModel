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
#############################################################################################
# Filename: Model.R
#
# This file runs the main models associated with the Hospital CRE Intervention Assessment Model (hCREiAM)
# 
# The model assumes a runtime environment of lsoda the ODE solver
# Inputs: 
# Time - vector of time for the model to run
# State - Current state of the model (inital population fractions)
# Pars - Parameters needed for the model (see parameters file)
# DTVals - The structure of the two-ward model (this ensures the state variables are correctly aligned as two wards)
# 
CREmodel <- function(Time, State, Pars,DTVals) {
	
    with(as.list(c(Pars)), {
        # First set up the structure of the population
        X <- copy(DTVals)
        numrows = nrow(DTVals)
        numvals = length(DTVals[1,])
        st = 1
        # Now copy the current state of the system into the structure
        for(i in 1:numrows) {
            X[i,] = State[st:(st+numvals-1)]
            st = st+numvals
        }
        
        # Calculate the transmission parameters
        # First the transmission of the nurses and doctors
        # Equations 5 & 6 from the report with eq. 16 incorporated
        bn <- (Cnp^2*Qnp*Qpn)/((zeta + sigmaN * sigN1) + Qpn*Cnp)
        bd <- (Cdp^2*Qdp*Qpd)/((zeta + sigmaD * sigD1) + Qpd*Cdp)
        
        # Equations 5 & 6 from the report with eq. 12 & 13 incorporated
        Krn <- bn*rn*(X$Cu + sigCP*sigD*X$Cd + xi*sigCP*sigD*X$I)
        Krd <- bd*rd*(X$Cu + sigCP*sigD*X$Cd + xi*sigCP*sigD*X$I)
        
        # Now the transmission from the environment
        # These equations combine all healthcare workers
        rh  <- rd + rn
        Chp <- Cdp + Cnp
        Qhp <- (Qnp*(Cnp/rn) + Qdp*(Cdp/rd))/((Cnp/rn) + (Cdp/rd))
        
        # Equation 8 from the report
        by  <- ((X$Y/patients) / (IC + (X$Y/patients))) * Chp*rh*betaY*Qhp
        bt  <- ((X$Tx/patients)/ (IC + (X$Tx/patients)))* Chp*rh*betaT*Qhp
        bm  <- ((X$M/patients) / (IC + (X$M/patients))) * Chp*rh*betaM*Qhp
        bl  <- ((X$L/patients) / (IC + (X$L/patients))) * Chp*rh*betaL*Qhp
        bw  <- ((X$W/patients) / (IC + (X$W/patients))) * Chp*rh*betaW*Qhp

        Kre <- by + bt + bm + bl + bw 
        
        # Combine environment and HCWs
        K   <- Krn + Kre + Krd 
        
        # Re-state the LOS parameter
        gamma = 1/stay
        
        # Calculate the rates of transfer between the two wards
        Tgi = X[,1:5]*trans
        Xsum = colSums(Tgi)
        Tig = sweep(Tgi,2,Xsum)*-1
    
        # We assume a constant population, this ensures that the total population and the wards stay the same size
        intake = rowSums(X[,1:5])*gamma + X$I*delta - rowSums(Tig[,1:5]) + rowSums(Tgi[,1:5])
        
        # These are the primary ODE equations (equations in 1 from the report) for the population states
        dS  <- aS*intake  + Tig$S + omega*X$XR + rhoO*(X$Cu + X$Cd) - Tgi$S - X$S*(K + sigAe*pi + gamma)
        
        dXR <- aXR*intake + Tig$XR + sigAe*pi*(X$S + phi*(X$Cu + sigAm*X$Cd)) + X$I*sigAp*rhoR*(1 - epsilon) - Tgi$XR - X$XR*(psi*K + omega + gamma)
        
        dCu <- (1 - OM)*aC*intake  + Tig$Cu*(1 - OM) + (1 - theta)*K*(X$S + psi*X$XR) + X$I*sigAp*rhoR*epsilon - Tgi$Cu - X$Cu*(tau + sigAe*phi*pi + rhoO + PS + gamma)
        dCd <- OM*aC*intake + Tig$Cu*OM + Tig$Cd + PS*X$Cu - Tgi$Cd - X$Cd*(tau + sigAm*sigAe*phi*pi + rhoO + gamma)
        
        dI  <- aI*intake  + Tig$I + tau*(X$Cu+X$Cd) + K*theta*(X$S + psi*X$XR) - Tgi$I - X$I*(sigAp*rhoR + gamma + delta)
        
        # These are the primary ODE equations for the environmental states (equations in 9 from the report)
        dY  <- (X$Cu + sigD*X$Cd + sigD*xi*X$I)*eta*patients*Chp*(1/rh) - X$Y*(mu + sigY + sigY1)
        dTx <- (X$Cu + sigD*X$Cd + sigD*xi*X$I)*eta*patients*Chp*(1/rh) - X$Tx*(mu + sigT + sigT1)
        dM  <- (X$Cu + sigD*X$Cd + sigD*xi*X$I)*eta*patients*Chp*(1/rh)*betaM - X$M*(mu + Pc*(sigMc + sigMc1) + Psc*(sigMsc + sigMsc1) + Pnc*(sigMnc + sigMnc1))
        dL  <- (X$Cu + sigD*X$Cd + sigD*xi*X$I)*eta*patients*Chp*(1/rh)*betaL - X$L*(mu + sigL + sigL1)
        dW  <- (X$Cu + sigD*X$Cd + sigD*xi*X$I)*eta*patients*Chp*(1/rh) - X$W*(mu/450 + sigW + sigW1)

        # These maintain additional states for computing important outcome variables including deaths and infections
        dCH <- (1 - theta)*K*(X$S + psi*X$XR) + X$I*sigAp*rhoR*epsilon
        dCI <- aC*intake
        
        dIH <- tau*(X$Cu+X$Cd) + K*theta*(X$S + psi*X$XR)
        dII <- aI*intake
        
        dIe <- Kre*theta*(X$S + psi*X$XR)
        dIc <- (Krn + Krd)*theta*(X$S + psi*X$XR) 
        
        dDeaths <- X$I*delta
        
        # Now return the population structure to a single vector
        for(i in 1:numrows) {
            w = c(dS[i], dXR[i], dCu[i], dCd[i], dI[i], dY[i], dTx[i], dM[i], dL[i], dW[i], dCH[i], dCI[i], dIH[i], dII[i], dIe[i], dIc[i], dDeaths[i])
            if(i==1){
                y = w
            } else {
                y = c(y,w)
            }
        }

        # Regular return call for LSODA ODE solver
        return(list(y))
         
    })
}

# The following model is nearly the same as above, but implements patient cohorting (section 2.2.5 of the report)
# Inputs are the same
CREmodel_Cohorting <- function(Time, State, Pars,DTVals) {
  
  with(as.list(c(Pars)), {
    # First set up the structure of the population
    X <- copy(DTVals)
    numrows = nrow(DTVals)
    numvals = length(DTVals[1,])
    st = 1
    # Now copy the current state of the system into the structure
    for(i in 1:numrows) {
      X[i,] = State[st:(st+numvals-1)]
      st = st+numvals
    }
    
    # Calculate the transmission parameters
    # First the transmission of the nurses and doctors
    bn <- (Cnp^2*Qnp*Qpn)/((zeta + sigmaN * sigN1) + Qpn*Cnp)
    bd <- (Cdp^2*Qdp*Qpd)/((zeta + sigmaD * sigD1) + Qpd*Cdp)
    
    Krn <- bn*rn*(X$Cu + Fn*(sigCP*sigD*X$Cd1 + xi*sigCP*sigD*X$I1)) + bn*rn*(1-Fn)*(sigCP*sigD*X$Cd2 + xi*sigCP*sigD*X$I2)
    Krd <- bd*rd*(X$Cu + sigCP*sigD*(X$Cd1 + X$Cd2) + xi*sigCP*sigD*(X$I1 + X$I2))
    
    # Now the transmission from the environment
    rh  <- rd + rn
    Chp <- Cdp + Cnp
    Qhp <- (Qnp*(Cnp/rn) + Qdp*(Cdp/rd))/((Cnp/rn) + (Cdp/rd))
    
    by  <- ((X$Y/patients) / (IC + (X$Y/patients))) * Chp*rh*betaY*Qhp
    bt  <- ((X$Tx/patients)/ (IC + (X$Tx/patients)))* Chp*rh*betaT*Qhp
    bm  <- ((X$M/patients) / (IC + (X$M/patients))) * Chp*rh*betaM*Qhp
    bl  <- ((X$L/patients) / (IC + (X$L/patients))) * Chp*rh*betaL*Qhp
    bw  <- ((X$W/patients) / (IC + (X$W/patients))) * Chp*rh*betaW*Qhp
    
    # Combine environment and HCWs
    Kre <- by + bt + bm + bl + bw 
    
    K   <- Krn + Kre + Krd 
    
    # Re-state the LOS parameter
    gamma = 1/stay
    
    # Calculate the rates of transfer between the two wards
    Tgi = X*trans
    Xsum = colSums(Tgi)
    Tig = sweep(Tgi,2,Xsum)*-1
    
    # We assume a constant population, this ensures that the total population and the wards stay the same size
    intake = rowSums(X[,1:7])*gamma + rowSums(X[,6:7])*delta - rowSums(Tig[,1:7]) + rowSums(Tgi[,1:7])
    
    # These are the primary ODE equations for the population states
    dS  <- aS*intake  + Tig$S + omega*X$XR + rhoO*(X$Cu + (X$Cd1 + X$Cd2)) - Tgi$S - X$S*(K + sigAe*pi + gamma)
    
    dXR <- aXR*intake + Tig$XR + sigAe*pi*(X$S + phi*(X$Cu + sigAm*(X$Cd1 + X$Cd2))) + (X$I1 + X$I2)*sigAp*rhoR*(1 - epsilon) - Tgi$XR - X$XR*(psi*K + omega + gamma)
    
    dCu <- (1 - OM)*aC*intake  + Tig$Cu*(1 - OM) + (1 - theta)*K*(X$S + psi*X$XR) + (X$I1 + X$I2)*sigAp*rhoR*epsilon - Tgi$Cu - X$Cu*(tau + sigAe*phi*pi + rhoO + PS + gamma)
    
    dCd1<- OM*aC*intake*(1 - f) + Tig$Cu*OM*(1 - f) + Tig$Cd1 + PS*X$Cu*(1 - f) - Tgi$Cd1 - X$Cd1*(tau + sigAm*sigAe*phi*pi + rhoO + gamma)
    dCd2<- OM*aC*intake*f + Tig$Cu*OM*f + Tig$Cd2 + PS*X$Cu*f - Tgi$Cd2 - X$Cd2*(tau + sigAm*sigAe*phi*pi + rhoO + gamma)
    
    dI1 <- aI*intake*(1 - f)  + Tig$I1 + tau*(X$Cu+(X$Cd1 + X$Cd2))*(1 - f) + K*theta*(X$S + psi*X$XR)*(1 - f) - Tgi$I1 - X$I1*(sigAp*rhoR + gamma + delta)
    dI2 <- aI*intake*f  + Tig$I2 + tau*(X$Cu+(X$Cd1 + X$Cd2))*f + K*theta*(X$S + psi*X$XR)*f - Tgi$I2 - X$I2*(sigAp*rhoR + gamma + delta)
    
    # These are the primary ODE equations for the environmental states
    dY  <- (X$Cu + sigD*(X$Cd1 + X$Cd2) + sigD*xi*(X$I1 + X$I2))*eta*patients*Chp*(1/rh) - X$Y*(mu + sigY + sigY1)
    dTx <- (X$Cu + sigD*(X$Cd1 + X$Cd2) + sigD*xi*(X$I1 + X$I2))*eta*patients*Chp*(1/rh) - X$Tx*(mu + sigT + sigT1)
    dM  <- (X$Cu + sigD*(X$Cd1 + X$Cd2) + sigD*xi*(X$I1 + X$I2))*eta*patients*Chp*(1/rh)*betaM - X$M*(mu + Pc*(sigMc + sigMc1) + Psc*(sigMsc + sigMsc1) + Pnc*(sigMnc + sigMnc1))
    dL  <- (X$Cu + sigD*(X$Cd1 + X$Cd2) + sigD*xi*(X$I1 + X$I2))*eta*patients*Chp*(1/rh)*betaL - X$L*(mu + sigL + sigL1)
    dW  <- (X$Cu + sigD*(X$Cd1 + X$Cd2) + sigD*xi*(X$I1 + X$I2))*eta*patients*Chp*(1/rh) - X$W*(mu/450 + sigW + sigW1)
    
    # These maintain additional states for computing important outcome variables including deaths and infections
    dCH <- (1 - theta)*K*(X$S + psi*X$XR) + (X$I1 + X$I2)*sigAp*rhoR*epsilon
    dCI <- aC*intake
    
    dIH <- tau*(X$Cu+(X$Cd1 + X$Cd2)) + K*theta*(X$S + psi*X$XR)
    dII <- aI*intake
    
    dIe <- Kre*theta*(X$S + psi*X$XR)
    dIc <- (Krn + Krd)*theta*(X$S + psi*X$XR) 
    
    dDeaths <- (X$I1 + X$I2)*delta
    
    # Now return the population structure to a single vector
    for(i in 1:numrows) {
      w = c(dS[i], dXR[i], dCu[i], dCd1[i], dCd2[i], dI1[i], dI2[i], dY[i], dTx[i], dM[i], dL[i], dW[i], dCH[i], dCI[i], dIH[i], dII[i], dIe[i], dIc[i], dDeaths[i])
      if(i==1){
        y = w
      } else {
        y = c(y,w)
      }
    }
    
    # Regular return call for LSODA ODE solver
    return(list(y))
    
  })
}

