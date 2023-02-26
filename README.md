# BMCM2022

Winning submission for problem 2 for Brown's local MCM competition (Mathematical Contest in Modeling); problem 2 was concerned with constructing
an ideal concert hall in 2-dimensions, maximizing for both listener experience and physical capacity under construction cost constraints.

We did this by solving the 2d wave equation over convex domains in $\mathbb{R}^2$ with reflecting boundary conditions, viscosity, and arbitrary
point source configurations (to model being able to move speakers/ performer placement) to simulate the propagation of sound
waves through an arbitrary concert hall. After defining a metric with which would could maximize average sound *clarity* throughout the hall,
which we correlated with sound *quality*, it became a topology optimization problem.
