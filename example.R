# This script demonstrates how to use the function calc.GHI(), defined in 'analytical_approach.R', on a simple example

# Required R-packages
library("maptools")
library("insol")
library("ggplot2")

# Location and module's position
xco <- 7.834333  # E (for Freiburg, Germany)
yco <- 48.009026 # N (for Freiburg, Germany)
tilt <- 45
azimuth <- 0 # South

# Read the input data
fname <- file.choose() # Select the input file: 'Input_example.rds'
data <- readRDS(fname) # data is a matrix with columns: Time, GHI, GTI.
time <- data[,1]       # Timesteps
Gh <- data[,2]         # Global horizontal irradiance
Gc <- data[,3]         # Global tilted irradiance


# Solar geometry (sun position) is calculated
doy<- daydoy(time)
pos <- matrix(c(xco,yco), nrow=1)
sunpos <- solarpos(pos,time)
gamma_s <- (sunpos[,2])/180*pi  #[Rad]
alpha_s <- (sunpos[,1])/180*pi  #[Rad]
Z <- (90-gamma_s*180/pi)/180*pi # sun's zenith in [Rad]

beta <- (tilt)/180*pi           #[Rad]
alpha <- (azimuth)/180*pi       #[Rad]

# Angle of incidence (AOI) is calculated
theta <- acos(-cos(gamma_s)*sin(beta)*cos(alpha_s-alpha)+sin(gamma_s)*cos(beta))

# The horizontal extraterestrial irradiance (I0h) is calculated 
DayAngle <- 2*pi*(doy-1)/365
re <- 1.00011+0.034221*cos(DayAngle)+(0.00128)*sin(DayAngle)+0.000719*cos(2*DayAngle)+(7.7E-5)*sin(2*DayAngle)
I0 <- re*1366.1
I0h<- I0*cos(Z)


# Quality-control sequence (based on the available data):
qual1 <- which((theta*180/pi)<85)     # angle of incidence (AOI) < 85°
qual2 <- which((gamma_s*180/pi)>5)    # solar elevation angle > 5°
qual3 <- which((Gh>0)&(Gc>0))
qual4 <- which(Gh<(1.5*I0*(cos(Z))^1.2 + 100))
qual <- intersect(qual1,qual2)
qual <- intersect(qual,qual3)
qual <- intersect(qual,qual4)


# Calculation of Gh from Gc using the new approach, but only for the data that pass the quality control
Gh_sim1 <- c()
Gh_sim1 <- calc.GHI(Gc=Gc[qual],time=time[qual],xcord=xco,ycord=yco,tilt=tilt,azimuth=azimuth,
                    albedo=0.2,gamma_s=gamma_s[qual],alpha_s=alpha_s[qual],version="A")

Gh_sim <- 0*Gh          # values which are not ...
Gh_sim[qual] <- Gh_sim1 # simulated are set to zero


# Plot results
plot(time,Gc,cex=0.1,col=alpha("grey"),type="l",ylim=c(0, 900),xlab="Time [h]", 
     ylab="Gh/Gc [W/m�]")
lines(time,Gh,cex=0.1,col=alpha("black"))
lines(time,Gh_sim,cex=0.1,col=alpha("blue"))

legend("topright",legend=c("Gc", "Gh meas", "Gh sim"),
       fill=c("grey", "black", "blue"))






