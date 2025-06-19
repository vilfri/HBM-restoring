# Guidelines for Temperature and Salinity Restoring in HBM

## Introduction
In the HBM ocean model, temperature and salinity restoring is a technique used to maintain realistic conditions in the simulated ocean. This guideline provides an overview of the temperature and salinity restoring process in HBM and offers instructions on how to implement it effectively.

## Temperature Restoring
Temperature restoring is the process of nudging the simulated ocean temperature towards a prescribed reference field, i.e. one decade climatological field in this case. This helps to maintain realistic temperature distributions and prevent the model from drifting away from observed conditions. Here's how to use temperature restoring in HBM:

1. Obtain CMEMS climatology file. E.g., BALTICSEA_MULTIYEAR_PHY_003_011.
2. Set up the temperature restoring parameters in the HBM configuration file. This typically involves specifying the restoring timescale and the strength of the nudging term. In the test phase, we use the summer restoring strategy. That means we only restoring the temperature during summer seasons when the Baltic inflow is not critical, i.e. June, July, August. The nudging period is one month. 
3. Configure the model to apply the temperature restoring only in specific regions of interest or throughout the entire domain, depending on your research objectives. In the test run, we restore the fields for the layer deeper than 50m for the whole Baltic Sea, including the south of Danish straits in the Inner Danish Water (IDW) domain. 
4. Initialize the model with the reference temperature field before starting the simulation.
5. During the simulation, the model will continuously adjust the simulated temperatures towards the reference field according to the specified restoring timescale and strength.

## Salinity Restoring
Salinity restoring is similar to temperature restoring but focuses on maintaining specific Baltic inflow events more realistic in the model. It nudges the simulated ocean salinity towards a prescribed reference field. The implementation would be same steps as above. Here one can choose to only restore salinity for researc purpose and specific for the Baltic Sea. However, it is always recommend to restore both temperature and salinity. 

## Tips and Best Practices
- Ensure that the reference temperature and salinity fields are of high quality and representative of the target region.
- Choose appropriate values for the restoring timescale and strength based on the timescales of the phenomena you are studying and the model's grid resolution.
- Monitor the model's behavior during the simulation to ensure that the temperature and salinity restoring is not causing unrealistic or excessive changes in the model results.
- Validate the model's output against available observations to assess the effectiveness of the temperature and salinity restoring technique.
- Experiment with different configurations and parameter values to find the optimal setup for your specific research goals.

## Conclusion
Temperature and salinity restoring in the HBM ocean model can help maintain realistic ocean conditions and prevent the model from drifting away from observed data. By following the guidelines provided in this document and considering the best practices, you can effectively utilize temperature and salinity restoring in your HBM simulations.

