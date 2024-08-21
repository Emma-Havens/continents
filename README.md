# Numerical Model of Early Continent Formation

Water content and mantle temperature are both highly influential and poorly constrained factors that affect the time of formation of early continents, their persistence at the surface, and their rate of accretion. One challenge in understanding the ultimate effect of these factors is their competing influences. Water content can both decrease viscosity and lower melting points, which has competing effects in terms of crustal preservation. High temperatures lowers rock viscosity while also decreasing water content, making its ultimate effect unclear. This repo contains a box model I built to explore how varying these factors, along with silica content, changes the timing and speed of continental formation and the survivability of generated crust.

For a highly detailed description of this model, including the geological background, mathematical foundation, numerical methods implemented, and a description of all run steps, see Numerical_Model_of_Early_Continent_Formation.pdf.

To run the simulation, run mantle_mixing_therm_res.m in MATLAB.

The MATLAB model was modified from a pre-existing model by Elvira Mulyukova. This repo also contains an earlier version of the model written in Python, which I wrote from scratch.
