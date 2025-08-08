This is the OpenFoam solver for reactive transport with dissolution.

It solves for time-dependent flow (Stokes or inertial) and reactant transport. The fracture domain will be deformed based on the mineral etching amount by the acid.

An example simulation case is waterCaseOnePerfTemplate.

This version was developed with OpenFOAM-v1706. It should also work with v1712.

This solver was modified after the solver developed by Starchenko, V., Marra, C. J., and Ladd, A. J. C. (2016). Three-dimensional Simulations of Fracture Dissolution. Journal of Geophysical Research: Solid Earth, 121(9):6421-6444.

Main modifications include:

1. Keep constant flux at inlet instead of outlet.
2. Add time derivative in acid transport equation and make it time-dependent.
3. Remove Brinkman term in momentum equation.
4. Remove rescaling operation after solving flow equations.
5. For restart, the first step does not update mesh.

If you use this solver in your work, please consider citing my following publications:

1. Dong, R., Wheeler, M. F., Ma, K., and Su, H. 2020. A 3D Acid Transport Model for Acid Fracturing Treatments with Viscous Fingering. Paper presented at the SPE Annual Technical Conference and Exhibition. SPE-201465-MS. doi: https://doi.org/10.2118/201465-MS.
2. Dong, R., Wheeler, M. F., Su, H., and Ma, K. 2021. Modeling Multistage Acid Fracturing Treatments in Carbonate Reservoirs. Paper presented at the SPE Hydraulic Fracturing Technology Conference and Exhibition. SPE-204197-MS. doi: https://doi.org/10.2118/204197-MS.
3. Dong, R., Wheeler, M. F., Su, H., and Ma, K. 2021. Modeling Acid Fracturing Treatments with Multi-Stage Alternating Injection of Pad and Acid Fluids. Paper presented at the SPE Reservoir Simulation Conference. SPE-203985-MS. doi: https://doi.org/10.2118/203985-MS.
4. Dong, R., Wheeler, M. F., Su, H., and Ma, K. 2021. Modeling Acid Fracturing Treatments in Heterogeneous Carbonate Reservoirs. Paper presented at the SPE International Conference on Oilfield Chemistry. SPE-203985-MS. doi: https://doi.org/10.2118/204304-MS.
