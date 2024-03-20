# Surface Energy Balance Model 
The Surface Energy Balance Algorithm over Land (SEBAL) is a widely used single-source energy balance model. This model is built on a theory that employs physical factors as well as empirical relationships. To address the limitation of the SEBAL, a METRIC-based evapotranspiration (ET) estimation model was developed. This model uses an automated internal calibration method that is almost like the SEBAL model for calculating sensible and latent energy.Both model uses climatic data to calculate Reference evapotranspiration (ET<sub>o</sub>) and it combines ground-based reference data with satellite data input to calculate the final ET. 

## The Surface Energy Balance Algorithm over Land (SEBAL)

The evaporative fraction at each pixel of the image is calculated using Rn, G, and H computations at the time of satellite overpass as:

<img width="170" alt="Screenshot 2024-03-20 103422" src="https://github.com/AthiraNG/Metric-Sebal/assets/129937610/a684b626-dfbe-4eda-a420-782f795231d9">

The EF is the evaporative fraction, which is calculated for the time of satellite overpass is assumed to be constant over the 24-hour period of the image acquisition day. As a result, the daily ET at each image pixel is calculated as:

<img width="300" alt="Screenshot 2024-03-20 103845" src="https://github.com/AthiraNG/Metric-Sebal/assets/129937610/4f48808f-cf76-4209-b73b-68a292772bd0">

where Rn<sub>24</sub> is the daily average net radiation (W m<sup>-2</sup>), and G<sub>24</sub> is the daily average soil heat flux (W ms<sup>-2</sup>) for soil and vegetation surfaces, which are usually assumed to be zero.

## Mapping Evapotranspiration with Internalized Calibration (METRIC)

The METRIC model is an image processing model for remote sensing that calculates instantaneous ET values as a residual of the surface energy balance equation.
ET is calculated at each pixel at the time of satellite overpass as:

<img width="210" alt="Screenshot 2024-03-20 100335" src="https://github.com/AthiraNG/Metric-Sebal/assets/129937610/fb0099f8-2ebc-4010-9f8f-c2b8be7fe41b">

Here ET<sub>ins</sub> is instantaneous ET (mm h <sup>-1</sup>), λ is the latent heat of vaporization (J kg<sup>-1</sup>), and ρw is the density of water (kg m<sup>-3</sup>).

The reference ET fraction (ETr<sub>F</sub>) is calculated as:

<img width="150" alt="Screenshot 2024-03-20 101929" src="https://github.com/AthiraNG/Metric-Sebal/assets/129937610/0e9dad80-58db-4c51-b796-d76fe207980b">

The ETr<sub>F</sub> computed for the time of satellite overpass is assumed to be the same as the ETrF over the 24-h average. Finally, the daily ET (ET<sub>24</sub>) at each pixel is computed as:

<img width="210" alt="Screenshot 2024-03-20 102834" src="https://github.com/AthiraNG/Metric-Sebal/assets/129937610/2ac0fc54-7445-4a1d-8c57-845f952662fd">

ETr<sub>24</sub> is the total 24-hour daily ETr (mm day<sup>-1</sup>) calculated using the standardized FAO-56 Penman-Monteith equation.

### Output

<img width="4500" alt="Picture1]" src="https://github.com/AthiraNG/Metric-Sebal/assets/129937610/51070d79-2009-4d41-80e4-a6e7b9eb5a69">

## Resources

#### [SEBAL model in the Nansi Lake Wetland of China](https://www.sciencedirect.com/science/article/pii/S0895717710005303)
#### [METRIC Model](https://www.researchgate.net/publication/228615269_Satellite-Based_Energy_Balance_for_Mapping_Evapotranspiration_With_Internalized_Calibration_METRIC_-_Model)
#### [Comparison of Four Different Energy Balance Models for Estimating Evapotranspiration](https://www.mdpi.com/2073-4441/8/1/9)
