################################################################################################################################################################################
#  Description:
#' Computes GHI from GTI using the new decomposition model in the POA and the inverse transposition model. Data quality control sequence is not a part of this script.    
#
#' @title         Analytical Calculation of GHI from GTI
#' 
#                 name                        type                    description  
#' @param         Gc                          num                     Global Irradiance on a Module [W/m²]
#' @param         time                        POSIXct                 Time is in ISO-format. E.g. ISOdate(year=2013, month=01, day=01, hour = 00, min = 0, sec = 0, tz = "GMT")
#' @param         xcord                       num                     Geographic coordinate. E.g. Freiburg: 7.841147
#' @param         ycord                       num                     Geographic coordinate. E.g. Freiburg: 47.997291
#' @param         tilt                        num                     Tilt angle of module [°]
#' @param         azimuth                     num                     azimuth angle of module [°]
#' @param         albedo                      num                     Albedo value, usually 0.2, if unknown
#' @param         alpha_s                     num                     azimuth angle of sun position [Rad]!!! Position South = 180°
#' @param         gamma_s                     num                     inclination angle of sun position [Rad]!!!
#' @param         version                     char                    Defines which version (A or B) of the new model will be used.
#'
#' @return                                    num                     Global Horizontal Irradiance [W/m²]
#' 
#' @examples
#' see 'example.R'
#' @export
################################################################################################################################################################################

calc.GHI <- function(Gc,time,xcord=8,ycord=48,tilt=25,azimuth=0,albedo=0.2,gamma_s=-999,alpha_s=-999,version="A"){
  
  # Required R-packages
  library("maptools")
  library("insol")
  
  # Solar geometry (sun position) is calculated, if not provided
  doy<- daydoy(time)
  if(gamma_s[1]==-999|alpha_s[1]==-999){
    pos <- matrix(c(xcord,ycord), nrow=1)
    sunpos <-solarpos(pos,time)
    gamma_s <- (sunpos[,2])/180*pi #[Rad]
    alpha_s <- (sunpos[,1])/180*pi #[Rad]
  }
  Z <- (90-gamma_s*180/pi)/180*pi # sun's zenith in [Rad]
  
  beta <- (tilt)/180*pi         #[Rad]
  alpha <- (azimuth)/180*pi     #[Rad]
  
  # Angle of incidence (AOI) is calculated
  theta <- acos(-cos(gamma_s)*sin(beta)*cos(alpha_s-alpha)+sin(gamma_s)*cos(beta))
  
  # The horizontal extraterestrial irradiance (I0h) is calculated 
  DayAngle <- 2*pi*(doy-1)/365
  re <- 1.00011+0.034221*cos(DayAngle)+(0.00128)*sin(DayAngle)+0.000719*cos(2*DayAngle)+(7.7E-5)*sin(2*DayAngle)
  I0 <- re*1366.1
  I0h<- I0*cos(Z)
  

  ############# Decomposition model in the POA #####################################################
  
  # Values in the POA:
  I0c <- I0h*cos(theta)/sin(gamma_s)
  Kt_poa <- Gc/I0c
  
  # 3 linear segments: 
  seg1 <- which(Kt_poa<=0.3)                # segment 1
  seg2 <- which((Kt_poa>0.3)&(Kt_poa<0.78)) # segment 2
  seg3 <- which(Kt_poa>=0.78)               # segment 3
  
  Kt1 <- Kt_poa[seg1]
  Z1 <- Z[seg1]
  Kt2 <- Kt_poa[seg2]
  Z2 <- Z[seg2]
  Kt3 <- Kt_poa[seg3]
  Z3 <- Z[seg3]
  
  
  # Calculations of the model coefficients:
  if(version=="A"){
    # Version A of the new model
    
    y_a1 <- -1.79E-05*azimuth^2 + 0.0001*azimuth + 0.7635 # y and h are auxiliary variables
    h_a1 <- -0.0021*tilt + 0.9604
    a1 <- (tilt/90)*(y_a1-0.7635)+h_a1
    y_b1 <- -4.5E-05*azimuth^2 - 0.0007*azimuth - 0.5968
    h_b1 <- 5.21E-05*tilt^2 - 0.0111*tilt - 0.0191
    b1 <- (tilt/90)*(y_b1+0.5968)+h_b1
    y_c1 <- 4.27E-05*azimuth^2 + 0*azimuth + 0.3956
    h_c1 <- 0.004*tilt + 0.0367
    c1 <- (tilt/90)*(y_c1-0.3956)+h_c1
    
    y_a2 <- -2.72E-05*azimuth^2 + 0.0002*azimuth + 0.7784
    h_a2 <- -0.0069*tilt + 1.3824
    a2 <- (tilt/90)*(y_a2-0.7784)+h_a2
    y_b2 <- 1.49E-05*azimuth^2 - 0.0013*azimuth - 1.4297
    h_b2 <- -11.15E-05*tilt^2 + 0.0149*tilt - 1.8707
    b2 <- (tilt/90)*(y_b2+1.4297)+h_b2
    y_c2 <- 3.17E-05*azimuth^2 + 0.0007*azimuth + 0.7694
    h_c2 <- 6.55E-05*tilt^2 - 0.0003*tilt + 0.2692
    c2 <- (tilt/90)*(y_c2-0.7694)+h_c2
    
    y_a3 <- -3.01E-05*azimuth^2 - 0.0002*azimuth + 0.2265
    h_a3 <- 2.57E-05*tilt^2 + 0.0008*tilt - 0.049
    a3 <- (tilt/90)*(y_a3-0.2265)+h_a3
    y_b3 <- 0.68E-05*azimuth^2 + 0.0008*azimuth + 0.509
    h_b3 <- -9.19E-05*tilt^2 + 0.0075*tilt + 0.5763
    b3 <- (tilt/90)*(y_b3-0.509)+h_b3
    y_c3 <- 3.72E-05*azimuth^2 - 0.0007*azimuth - 0.4251
    h_c3 <- 8.76E-05*tilt^2 - 0.0104*tilt - 0.1947
    c3 <- (tilt/90)*(y_c3+0.4251)+h_c3
  }else{
    # Version B of the new model
    
    a1 <- -1.79E-05*azimuth^2 + 0.0001*azimuth + (-0.0021*tilt + 0.9604)
    b1 <- -4.5E-05*azimuth^2 - 0.0007*azimuth + (5.21E-05*tilt^2 - 0.0111*tilt - 0.0191)
    c1 <- 4.27E-05*azimuth^2 + 0*azimuth + (0.004*tilt + 0.0367)
    
    a2 <- -2.72E-05*azimuth^2 + 0.0002*azimuth + (-0.0069*tilt + 1.3824)
    b2 <- 1.49E-05*azimuth^2 - 0.0013*azimuth + (-11.15E-05*tilt^2 + 0.0149*tilt - 1.8707)
    c2 <- 3.17E-05*azimuth^2 + 0.0007*azimuth + (6.55E-05*tilt^2 - 0.0003*tilt + 0.2692)
    
    a3 <- -3.01E-05*azimuth^2 - 0.0002*azimuth + (2.57E-05*tilt^2 + 0.0008*tilt - 0.049)
    b3 <- 0.68E-05*azimuth^2 + 0.0008*azimuth + (-9.19E-05*tilt^2 + 0.0075*tilt + 0.5763)
    c3 <- 3.72E-05*azimuth^2 - 0.0007*azimuth + (8.76E-05*tilt^2 - 0.0104*tilt - 0.1947)
  }
  
  
  # Calculation of Kd in the POA
  Kd_calc <- c()
  Kd_calc[seg1] <- a1+b1*Kt1+c1*cos(Z1) # segment 1
  Kd_calc[seg2] <- a2+b2*Kt2+c2*cos(Z2) # segment 2
  Kd_calc[seg3] <- a3+b3*Kt3+c3*cos(Z3) # segment 3
  Kd_calc[Kd_calc<0]<-0
  Kd_calc[Kd_calc>1]<-1
  
  # Calculation of Bc, Bh and Dif_c
  Bc_calc <- (1-Kd_calc)*Gc
  amplifyDirHI <- cos(theta)/sin(gamma_s)
  Bh_calc <- Bc_calc/amplifyDirHI
  Dc_calc <- Gc-Bc_calc
  
  ############# Inverse transposition model #####################################################
  
  # Calculation of auxiliary variables:
  Surfazimuth <- (alpha*180/pi + 180)/180*pi
  COSTT <- cos(beta)*cos(Z) + sin(beta)*sin(Z)*cos(alpha_s-Surfazimuth)
  
  Factor1 <- COSTT
  Factor1[which(Factor1<0)] <- 0 
  Factor2 <- cos(Z)
  Factor2[which(Factor2<0.01745)] <- 0.01745
  Rb <- Factor1/Factor2 # direct irradiance conversion factor
  Bh_calc[which(Bh_calc<0)] <- 0 # Don't take the square root of a negative number
  AI <- Bh_calc/I0h
  SCUBE <- (sin(beta*0.5))^3
  
  Ffac <- sqrt(Bh_calc/Gc)
  Rd <- AI*Rb+(1-AI)*0.5 *(1 +cos(beta))*(1+Ffac*SCUBE) # diffuse irradiance conversion factor
  Rr <- (albedo/2)*(1-cos(beta)) # reflected irradiance conversion factor
  
  
  # Calculation of Dh:
  Dh_calc <- c()
  Dh_calc <- (Dc_calc-Bh_calc*Rr)/(Rd+Rr)
  Dh_calc[which(Dh_calc<0)]<-0
  Dh_calc[which(Dh_calc>(2*Gc/(1+cos(beta))))]<-(2*Gc/(1+cos(beta)))[which(Dh_calc>(2*Gc/(1+cos(beta))))]
  
  # Calculation of Gh:
  Gh_calc <- c()
  Gh_calc <- Dh_calc + Bh_calc
  
  return(Gh_calc)
}



