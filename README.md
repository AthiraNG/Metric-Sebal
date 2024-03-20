## Surface energy balance model (SEBAL and METRIC)
The Surface Energy Balance Algorithm over Land (SEBAL) is a widely used single-source energy balance model. This model is built on a theory that employs physical factors as well as empirical relationships. To address the limitation of the SEBAL, a METRIC-based evapotranspiration (ET) estimation model was developed. This model uses an automated internal calibration method that is almost like the SEBAL model for calculating sensible and latent energy.Both model uses climatic data to calculate Reference evapotranspiration (ET<sub>o</sub>) and it combines ground-based reference data with satellite data input to calculate the final ET. The main key differences between SEBAL and METRIC include, in the SEBAL model, the H at the cold pixel is assumed to be zero, and the LE is calculated as the difference between net radiation and G, On the other hand, in METRIC, the ET at the hot pixel is set to zero, and the cold pixel is assumed to have an ET rate of 1.05 times the reference ET rate. The hot pixel in METRIC is defined as a non-agricultural area, where the surface temperature (Ts) is not influenced by crop transpiration.

# METRIC
The METRIC model is an image processing model for remote sensing that calculates instantaneous ET values as a residual of the surface energy balance equation.
ET is calculated at each pixel at the time of satellite overpass as:

<img width="220" alt="Screenshot 2024-03-20 100335" src="https://github.com/AthiraNG/Metric-Sebal/assets/129937610/fb0099f8-2ebc-4010-9f8f-c2b8be7fe41b">

Here ETins is instantaneous ET (mmh <sup>-1</sup>), λ is the latent heat of vaporization (Jkg<sup>-1</sup>), and ρw is the density of water (kgm<sup>-3</sup>).

The reference ET fraction (ETrF) is calculated as:

<img width="140" alt="Screenshot 2024-03-20 101929" src="https://github.com/AthiraNG/Metric-Sebal/assets/129937610/0e9dad80-58db-4c51-b796-d76fe207980b">

The ETrF computed for the time of satellite overpass is assumed to be the same as the ETrF over the 24-h average. Finally, the daily ET (ET24) at each pixel is computed as:

<img width="200" alt="Screenshot 2024-03-20 102834" src="https://github.com/AthiraNG/Metric-Sebal/assets/129937610/2ac0fc54-7445-4a1d-8c57-845f952662fd">

ETr 24 is the total 24-hour daily ETr (mmday<sup>-1</sup>) calculated using the standardized FAO-56 Penman-Monteith equation.

# SEBAL

The evaporative fraction at each pixel of the image is calculated using Rn, G, and H computations at the time of satellite overpass as:

