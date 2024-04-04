stellar_density_profiles.py##########################
# define bar density
##########################

# Sormani et al. 2020 model 3 NSD density
#def Calc_S20_NSD_rho(x,y,z):
def Calc_S20_NSD_rho(l,b,r):
    alpha   = 0.9
    rho1    = 222.885   # 1e10Msun/kpc^3
    R1      = 0.0050617 # kpc
    n1      = 0.7194    # dimensionless
    rho2    = 169.975   # 1e10Msun/kpc^3
    R2      = 0.0246    # kpc
    n2      = 0.7933    # dimensionless
    q       = 0.37      # dimensionless
    x,y,z,vx,vy,vz = lbr2xyz(l,b,r,0,0,0)
    a       = sqrt(x**2+y**2+(z/q)**2)
    rho     = alpha*rho1*exp(-(a/R1)**n1) + alpha*rho2*exp(-(a/R2)**n2)
    return rho

# Launhardt et al. 2002 NSD density
#def Calc_L02_NSD_rho(x,y,z):
def Calc_L02_NSD_rho(l,b,r):
        nR1     = 5.0     # dimensionless
        nz1     = 1.4     # dimensionless
        rho1    = 15.228  # 10^10Msun/kpc^3
        R1      = 0.12    # kpc
        z1      = 0.045   # kpc
        nR2     = 5.0     # dimensionless
        nz2     = 1.4     # dimensionless
        rho2    = 3.888   # 10^10Msun/kpc^3
        R2      = 0.22    # kpc
        z2      = 0.045   # kpc
        log2    = log(2)
        x,y,z,vx,vy,vz = lbr2xyz(l,b,r,0,0,0)
        R       = sqrt(x**2+y**2)
        A1      = rho1*exp(-log2*(  (R/R1)**nR1 + (np.abs(z)/z1)**nz1) )
        A2      = rho2*exp(-log2*(  (R/R2)**nR2 + (np.abs(z)/z2)**nz2) )
        rho     = A1 + A2
        return rho
    
# Chatzopoulos et al. 2015 NSC density; only important for central 10pc
#def Calc_C15_NSC_rho(x,y,z):
def Calc_C15_NSC_rho(l,b,r):
    gamma = 0.71    # dimensionless
    q     = 0.73    # dimensionless
    a0    = 0.0059  # kpc
    M     = 0.0061  # 10^10Msun
    x,y,z,vx,vy,vz = lbr2xyz(l,b,r,0,0,0)
    a     = sqrt(x**2+y**2+(z/q)**2)
    if a<0.01: #If below 10pc
        rho   = (3.0-gamma)*M/(4*pi*q)*a0/((a**gamma)*((a+a0)**(4.0-gamma)));
        return rho
    else:   #only important for <10pc (will greatly overestimate the total mass, as noted in Sormani 2020)
        return 0
    
# Launhardt et al. 2002 Bar density
def Calc_L02_Bar_rho(x,y,z):
    rho0  = 0.8
    ax    = 1.1
    ay    = 0.36
    az    = 0.22
    Cperp = 1.6
    Cpar  = 3.2
    Rperp = ((np.abs(x)/ax)**Cperp + (np.abs(y)/ay)**Cperp)**(1.0/Cperp)
    Rs    = (Rperp**Cpar + (np.abs(z)/az)**Cpar)**(1.0/Cpar)
    rho   = rho0*exp(-Rs)
    return rho

# Exponential Disk Density Model
def n(l,b,r):
        n_0=2.1 #pc^-3 #2.097732758293176
        h_thin = 350. #pc
        h_r = 3500. #pc
        h_thick=1500. #pc
        x,y,z,vx,vy,vz = lbr2xyz(l,b,r,0,0,0)
        R, theta, z = cart2pol(x, y, z)
        n_pc = n_0*(np.exp(-z/h_thin) + 0.02*np.exp(-z/h_thick))*np.exp(-R/h_r)  #in pc^-3
        n_kpc = n_pc / 1e9 #in kpc^-3
        return n_kpc