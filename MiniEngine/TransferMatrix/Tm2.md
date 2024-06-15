# Layered Materials with Transfer Matrices (In Real Time)

This is a work in progress implementation of the 2-flux
transfer matrix based layered material framework presented
by Randrianandrasana et al (2021) in "Transfer Matrix Based Layered Materials Rendering".

The core of this work so far is reimplementing the 2-flux
transfer matrix model as a pixel shader for a forward rendered raster engine written with DirectX 12. 
This primarily involves adapting the model for working with preintegrated
lighting. 

From this base, I investigate the performance of the framework
in a realtime context, and what (if any) simplifications can be made
to fit with realtime constraints. 

### References