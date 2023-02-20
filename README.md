# photons_in_water
Simulation of a photon pencil beam and interactions in a 2D water region

## Data files
**Attenuation coefficients:**
* Source: https://physics.nist.gov/PhysRefData/XrayMassCoef/ComTab/water.html
* `attenuation_coefficients.npy` is a 2D array with columns: Energy (MeV) | μ/ρ (cm²/g) | μ_en/ρ (cm²/g)

**Cross sections:**
* Source: https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
* `cross_sections.npy` is a 2D array with columns: Energy (MeV) | Rayleigh (cm²/g) | Compton (cm²/g) | photoelectric absorption (cm²/g) | pair production (cm²/g) | total (cm²/g)
