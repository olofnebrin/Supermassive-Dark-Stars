""" Primordial gas-phase chemistry module.

    REFERENCES:
    
    GP98: Galli & Palla (1998), 'The chemistry of the early Universe', Astronomy & Astrophysics, 335, 403-420

    Glo10: Glover+ (2010), 'Modelling CO formation in the turbulent interstellar medium', MNRAS, 404, 2-29 

    For13: Forrey (2013), 'Rate of formation of hydrogen molecules by three-body recombination during
                 primordial star formation', ApJL, 773, L25 

    Glo15: Glover (2015), 'Simulating the formation of massive seed black holes in the early Universe
                 - I. An improved chemical model', MNRAS, 451, 2082-2096

    Sm17: Smith+ (2017), 'GRACKLE: a chemistry and cooling library for astrophysics', MNRAS, 466, 2217-2234  """


function k_2body(T)
    """ 2-body reaction rates relevant for the chemistry
        of H+, H-, e-, and H2. Units: cm^3 s^-1 """

    T4 = T/1e4; TeV = T/11608.696

    # H+ + e- -> H + hν (Case B, Draine 2011, p. 138):
    
    k1 = 2.59e-13*(T4^(-0.833 - 0.034*log(T4)))

    # H0 + e- -> H- + hν (GP98, Glo15):

    k2 = 1.4e-18*(T^0.928)*exp(-T/16200.0)

    # H0 + H- -> H2 (Kre10, Sm17):

    c1 = 1.35e-9; c2 = 9.8493e-2; c3 = 3.2852e-1
    c4 = 5.5610e-1; c5 = 2.7710e-7; c6 = 2.1826
    c7 = 6.1910e-3; c8 = 1.0461; c9 = 8.9712e-11
    c10 = 3.0424; c11 = 3.2576e-14; c12 = 3.7741

    k3 = c1*(T^c2 + c3*(T^c4) + c5*(T^c6))/(1 
         + c7*(T^c8) + c9*(T^c10) + c11*(T^c12))

    # H- + H+ -> H0 + H0 (Glo10, Glo15):

    k4 = 2.4e-6*(T^(-0.5))*(1.0 + T/2e4)

    # H- + e- -> H0 + e- + e- (Glo10):

    y  = log(T/TeV)
    k5 = exp(-1.801849334e1
            + 2.36085220e0*y 
            - 2.82744300e-1*(y^2)
            + 1.62331664e-2*(y^3)
            - 3.36501203e-2*(y^4)
            + 1.17832978e-2*(y^5)
            - 1.65619470e-3*(y^6)
            + 1.06827520e-4*(y^7)
            - 2.63128581e-6*(y^8))

    return [k1 k2 k3 k4 k5]
end

function k_3body(T)
    """ 3-body reaction rates relevant for H2 formation
        as a function of temperature T. (cm^6 s^-1)

        Notes: For the rate k31 there is considerable uncertainty,
        as discussed in Grackle 2017. I have adopted the latest
        quantum mechanical calculation by Forrey 2013. For k32
        most sources seem to adopt the same rate, from Cohen &
        Westberg 1983. """

    # H0 + H0 + H0 -> H2 + H0 (For13, Sm17):

    k31 = 6e-32*(T^(-0.25)) + 2e-31*(T^(-0.5))

    # H0 + H0 + H2 -> H2 + H2 (Glo10, Sm17):

    k32 = 2.8e-31*(T^(-0.6))

    return [k31 k32]
end