This code can be used to compute the TOV equations for a spherically symmetric star having any equation of state. 
I have used geometrized units and used appropriate conversion factors to convert every quantity to powers of length (km)
In both cases, the critical density at the centre must be speified.
For polytropes, kappa and tau have to be initialised.
For other EOS, a table of rho and their corresponding p values are required as input.
In addition to mass, density and pressure profiles, I have also computed the profile of the tidal metric perturbation. The theory for the same is explained in the report.\\
P.S. An inconsistency was noted while comparing with the SpEC TOV solver's results with respect to the total mass and radius computed by the code. The inconsistency was found to arise due to the missing specific internal energy term and the correction was made for the case of polytropes. The new results are found to be consistent with SpEC. However, we are yet to add the correction for tabulated EOS.
