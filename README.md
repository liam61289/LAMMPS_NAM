# LAMMPS_NAM
A new 'fix' for LAMMPS that includes a version of the non-adiabatic model.

INSTALLATION:

    Download lammps-28Jun14.tar.gz from https://lammps.sandia.gov/tars/ (fix NAM may not work in the more up to date versions of LAMMPS - This is currently being worked on)

    Copy contents of this repository to src folder

    Make LAMMPS as per LAMMPS instructions including any other required packages


SYNTAX:

    To use 2TMD the syntax is as follows:

        fix ID group-ID nab 2 seed C_e rho_e kappa_e gamma_p gamma_s v_0 Nx Ny Nz Bx By Bz T_infile N T_outfile

    To use NAM the syntax is as follows:

        fix ID group-ID nab 1 seed C_e rho_e kappa_e chi_p chi_s cap E1 E2 Nx Ny Nz Bx By Bz T_infile N T_outfile

    Where,
    
        seed               =  A seed for the uniform number generator
        C_e                =  The electronic heat capacity (energy/(electron*temperature) units)
        rho_e              =  The electron density (electrons/volume units)
        kappa_e            =  The electron thermal conductivity (energy/(time*distance*temperature) units)
        chi_p              =  Chi for the EPC regime (unitless)
        chi_s              =  Chi for the ES regime (unitless)
        gamma_p            =  Gamma for the EPC regime (mass/time units)
        gamma_s            =  Gamma for the ES regime (mass/time units)
        v_0                =  Velocity where ES regime is applied from (velocity units)
        cap                =  A cap on the return force as it is sometimes abnormally high (5x gamma_p is a good choice) (mass/time units)
        E1                 =  The kinetic energy at which the fermi dirac function between chi_p and chi_s starts (energy units)
        E2                 =  The kinetic energy at which the fermi dirac function between chi_p and chi_s ends (energy units)
        Nx                 
        Ny
        Nz                 =  The number of electron cells that coincide with the material (must be at least 1x1x1)
        Bx
        By
        Bz                 =  The number of electron cells either side of the material (i.e. total cells in x direction 2Bx+Nx)
        T_infile           =  A file describing the initial electron temperature of each cell in the form x y z T u (where 'u' means the temperature is free to change during the simulation and 'f' fixes the temperature)
        N                  =  How often to dump the temperature into T_outfile (0 means no dumps)
        T_outfile          =  Where to save the dump file
