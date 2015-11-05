
# This file is autogenerated by gen_code_info.py.
# DO NOT EDIT MANUALLY, CHANGES WILL BE LOST

include "myconfig.pxi"


def features():
    """Returns list of features compiled into Espresso core"""

    f = []

    IF ENGINE == 1:
        f.append("ENGINE")

    IF EVENT_DEBUG == 1:
        f.append("EVENT_DEBUG")

    IF PTENSOR_DEBUG == 1:
        f.append("PTENSOR_DEBUG")

    IF NO_INTRA_NB == 1:
        f.append("NO_INTRA_NB")

    IF CG_DNA == 1:
        f.append("CG_DNA")

    IF GRID_DEBUG == 1:
        f.append("GRID_DEBUG")

    IF CONSTRAINTS == 1:
        f.append("CONSTRAINTS")

    IF ONEPART_DEBUG == 1:
        f.append("ONEPART_DEBUG")

    IF LB_GPU == 1:
        f.append("LB_GPU")

    IF VIRTUAL_SITES == 1:
        f.append("VIRTUAL_SITES")

    IF ROTATIONAL_INERTIA == 1:
        f.append("ROTATIONAL_INERTIA")

    IF LJ_WARN_WHEN_CLOSE == 1:
        f.append("LJ_WARN_WHEN_CLOSE")

    IF BOND_ANGLE_COSSQUARE == 1:
        f.append("BOND_ANGLE_COSSQUARE")

    IF BOND_ENDANGLEDIST_HARMONIC == 1:
        f.append("BOND_ENDANGLEDIST_HARMONIC")

    IF TK == 1:
        f.append("TK")

    IF BOND_ANGLE_HARMONIC == 1:
        f.append("BOND_ANGLE_HARMONIC")

    IF GAY_BERNE == 1:
        f.append("GAY_BERNE")

    IF POLY_DEBUG == 1:
        f.append("POLY_DEBUG")

    IF BOND_ENDANGLEDIST == 1:
        f.append("BOND_ENDANGLEDIST")

    IF RANDOM_DEBUG == 1:
        f.append("RANDOM_DEBUG")

    IF GAUSSIAN == 1:
        f.append("GAUSSIAN")

    IF ESK_DEBUG == 1:
        f.append("ESK_DEBUG")

    IF COMFIXED == 1:
        f.append("COMFIXED")

    IF INTER_RF == 1:
        f.append("INTER_RF")

    IF MMM1D_GPU == 1:
        f.append("MMM1D_GPU")

    IF HAT == 1:
        f.append("HAT")

    IF PARTICLE_DEBUG == 1:
        f.append("PARTICLE_DEBUG")

    IF BUCKINGHAM == 1:
        f.append("BUCKINGHAM")

    IF FORCE_DEBUG == 1:
        f.append("FORCE_DEBUG")

    IF DPD_MASS == 1:
        f.append("DPD_MASS")

    IF MOL_CUT == 1:
        f.append("MOL_CUT")

    IF EK_DEBUG == 1:
        f.append("EK_DEBUG")

    IF VERLET_DEBUG == 1:
        f.append("VERLET_DEBUG")

    IF NPT == 1:
        f.append("NPT")

    IF EK_DOUBLE_PREC == 1:
        f.append("EK_DOUBLE_PREC")

    IF GHMC == 1:
        f.append("GHMC")

    IF DPD_MASS_LIN == 1:
        f.append("DPD_MASS_LIN")

    IF STAT_DEBUG == 1:
        f.append("STAT_DEBUG")        

    IF OIF_LOCAL_FORCES == 1:
        f.append("OIF_LOCAL_FORCES")

    IF OIF_LOCAL_FORCES == 1:
        f.append("OIF_GLOBAL_FORCES")

    IF COMFORCE == 1:
        f.append("COMFORCE")

    IF LB_BOUNDARIES == 1:
        f.append("LB_BOUNDARIES")

    IF GHOST_FLAG == 1:
        f.append("GHOST_FLAG")

    IF HERTZIAN == 1:
        f.append("HERTZIAN")

    IF IMMERSED_BOUNDARY == 1:
        f.append("IMMERSED_BOUNDARY")

    IF DPD_MASS_RED == 1:
        f.append("DPD_MASS_RED")

    IF MULTI_TIMESTEP == 1:
        f.append("MULTI_TIMESTEP")

    IF EK_ELECTROSTATIC_COUPLING == 1:
        f.append("EK_ELECTROSTATIC_COUPLING")

    IF BOND_ANGLEDIST_HARMONIC == 1:
        f.append("BOND_ANGLEDIST_HARMONIC")

    IF CATALYTIC_REACTIONS == 1:
        f.append("CATALYTIC_REACTIONS")

    IF PARTIAL_PERIODIC == 1:
        f.append("PARTIAL_PERIODIC")

    IF LENNARD_JONES == 1:
        f.append("LENNARD_JONES")

    IF ELECTROKINETICS == 1:
        f.append("ELECTROKINETICS")

    IF HALO_DEBUG == 1:
        f.append("HALO_DEBUG")

    IF BOND_ANGLE == 1:
        f.append("BOND_ANGLE")

    IF LJ_ANGLE == 1:
        f.append("LJ_ANGLE")

    IF ESR_DEBUG == 1:
        f.append("ESR_DEBUG")

    IF GHOST_FORCE_DEBUG == 1:
        f.append("GHOST_FORCE_DEBUG")

    IF OVERLAPPED == 1:
        f.append("OVERLAPPED")

    IF LB_DEBUG == 1:
        f.append("LB_DEBUG")

    IF EK_BOUNDARIES == 1:
        f.append("EK_BOUNDARIES")

    IF VIRTUAL_SITES_RELATIVE == 1:
        f.append("VIRTUAL_SITES_RELATIVE")

    IF CELL_DEBUG == 1:
        f.append("CELL_DEBUG")

    IF LATTICE_DEBUG == 1:
        f.append("LATTICE_DEBUG")

    IF FFT_DEBUG == 1:
        f.append("FFT_DEBUG")

    IF LB_BOUNDARIES_GPU == 1:
        f.append("LB_BOUNDARIES_GPU")

    IF BMHTF_NACL == 1:
        f.append("BMHTF_NACL")

    IF MOLFORCES_DEBUG == 1:
        f.append("MOLFORCES_DEBUG")

    IF EK_REACTION == 1:
        f.append("EK_REACTION")

    IF ADDITIONAL_CHECKS == 1:
        f.append("ADDITIONAL_CHECKS")

    IF MODES == 1:
        f.append("MODES")

    IF EXTERNAL_FORCES == 1:
        f.append("EXTERNAL_FORCES")

    IF TABULATED == 1:
        f.append("TABULATED")

    IF LATTICE == 1:
        f.append("LATTICE")

    IF ELECTROSTATICS == 1:
        f.append("ELECTROSTATICS")

    IF FENE_DEBUG == 1:
        f.append("FENE_DEBUG")

    IF BOND_ANGLE_OLD == 1:
        f.append("BOND_ANGLE_OLD")

    IF COMM_DEBUG == 1:
        f.append("COMM_DEBUG")

    IF SMOOTH_STEP == 1:
        f.append("SMOOTH_STEP")

    IF LANGEVIN_PER_PARTICLE == 1:
        f.append("LANGEVIN_PER_PARTICLE")

    IF METADYNAMICS == 1:
        f.append("METADYNAMICS")

    IF SOFT_SPHERE == 1:
        f.append("SOFT_SPHERE")

    IF VIRTUAL_SITES_NO_VELOCITY == 1:
        f.append("VIRTUAL_SITES_NO_VELOCITY")

    IF LJCOS == 1:
        f.append("LJCOS")

    IF LENNARD_JONES_GENERIC == 1:
        f.append("LENNARD_JONES_GENERIC")

    IF VIRTUAL_SITES_THERMOSTAT == 1:
        f.append("VIRTUAL_SITES_THERMOSTAT")

    IF MASS == 1:
        f.append("MASS")

    IF INTEG_DEBUG == 1:
        f.append("INTEG_DEBUG")

    IF DP3M == 1:
        f.append("DP3M")

    IF MORSE_DEBUG == 1:
        f.append("MORSE_DEBUG")

    IF ROTATION_PER_PARTICLE == 1:
        f.append("ROTATION_PER_PARTICLE")

    IF NEMD == 1:
        f.append("NEMD")

    IF SD_NOT_PERIODIC == 1:
        f.append("SD_NOT_PERIODIC")

    IF P3M == 1:
        f.append("P3M")

    IF VIRTUAL_SITES_DEBUG == 1:
        f.append("VIRTUAL_SITES_DEBUG")

    IF GHOSTS_HAVE_BONDS == 1:
        f.append("GHOSTS_HAVE_BONDS")

    IF BOND_ANGLE_COSINE == 1:
        f.append("BOND_ANGLE_COSINE")

    IF MOLFORCES == 1:
        f.append("MOLFORCES")

    IF GHOST_DEBUG == 1:
        f.append("GHOST_DEBUG")

    IF EWALD_GPU == 1:
        f.append("EWALD_GPU")

    IF FFTW == 1:
        f.append("FFTW")

    IF BOND_CONSTRAINT == 1:
        f.append("BOND_CONSTRAINT")

    IF TUNABLE_SLIP == 1:
        f.append("TUNABLE_SLIP")

    IF HYDROGEN_BOND == 1:
        f.append("HYDROGEN_BOND")

    IF COS2 == 1:
        f.append("COS2")

    IF COULOMB_DEBYE_HUECKEL == 1:
        f.append("COULOMB_DEBYE_HUECKEL")

    IF LJCOS2 == 1:
        f.append("LJCOS2")

    IF SD == 1:
        f.append("SD")

    IF BOND_ANGLEDIST == 1:
        f.append("BOND_ANGLEDIST")

    IF MAGGS_DEBUG == 1:
        f.append("MAGGS_DEBUG")

    IF THERMOSTAT_IGNORE_NON_VIRTUAL == 1:
        f.append("THERMOSTAT_IGNORE_NON_VIRTUAL")

    IF DPD == 1:
        f.append("DPD")

    IF LJGEN_SOFTCORE == 1:
        f.append("LJGEN_SOFTCORE")

    IF SD_USE_FLOAT == 1:
        f.append("SD_USE_FLOAT")

    IF INTER_DPD == 1:
        f.append("INTER_DPD")

    IF VIRTUAL_SITES_COM == 1:
        f.append("VIRTUAL_SITES_COM")

    IF LJ_DEBUG == 1:
        f.append("LJ_DEBUG")

    IF LB == 1:
        f.append("LB")

    IF MPI_CORE == 1:
        f.append("MPI_CORE")

    IF BOND_VIRTUAL == 1:
        f.append("BOND_VIRTUAL")

    IF OLD_DIHEDRAL == 1:
        f.append("OLD_DIHEDRAL")

    IF SD_FF_ONLY == 1:
        f.append("SD_FF_ONLY")

    IF CUDA == 1:
        f.append("CUDA")

    IF LEES_EDWARDS == 1:
        f.append("LEES_EDWARDS")

    IF SD_DEBUG == 1:
        f.append("SD_DEBUG")

    IF LB_ELECTROHYDRODYNAMICS == 1:
        f.append("LB_ELECTROHYDRODYNAMICS")

    IF P3M_DEBUG == 1:
        f.append("P3M_DEBUG")

    IF UMBRELLA == 1:
        f.append("UMBRELLA")

    IF MORSE == 1:
        f.append("MORSE")

    IF FORCE_CORE == 1:
        f.append("FORCE_CORE")

    IF DIPOLES == 1:
        f.append("DIPOLES")

    IF THERMO_DEBUG == 1:
        f.append("THERMO_DEBUG")

    IF ROTATION == 1:
        f.append("ROTATION")

    IF SHANCHEN == 1:
        f.append("SHANCHEN")

    IF TRANS_DPD == 1:
        f.append("TRANS_DPD")

    IF LE_DEBUG == 1:
        f.append("LE_DEBUG")

    IF COLLISION_DETECTION == 1:
        f.append("COLLISION_DETECTION")

    IF TWIST_STACK == 1:
        f.append("TWIST_STACK")

    IF ASYNC_BARRIER == 1:
        f.append("ASYNC_BARRIER")

    IF EXCLUSIONS == 1:
        f.append("EXCLUSIONS")

    return sorted(f)
