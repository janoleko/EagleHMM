# Notes for the report

## Data

* longitude and latitude: Converted to Step length/ speed (m/s) and Turning angle (1 Hz)
* height above mean sealevel: Calculate first difference (1 Hz)
* possible covariates:
  * temperature
  * landform (categorical variable: 6 Landforms + NA)
  * heading might be useful
  * solar time (ToD in decimal)
  * day of observation
  * elevation

Challenges of the data: 30 Second intervals of 1 Hz data then 15 Minutes (900 seconds) breaks
As basically no state switching was observed within intervals we decided to compute summary statistics for each interval.

## Preprocessing & understanding the data

To understand the behaviour of the eagle we visualized the coordinates as well as the height movement patterns for each ~ 250-300 intervals.

We came up with 3 to 4 clearly distinguishable behaviours which the HMM could capture with 3 or 4 states:

* Soaring: The coordinate plots show circles and the height plot shows upward movement
* Gliding: The coordinate plot shows fast and directed movement while the height plot shows either downwards movement or no change in altitude. This could be captured in one or two states?
* Resting: The coordinate plot shows little movement with a rather random pattern while the height plot shows barely any changes in height´ß


### Summary statistics

* Step lenght/ speed: Mean &rarr; Gamma distributed
* Turning angle: abs(mean())/pi &rarr; from von Mises distr. to beta distribution \\
&rarr; Absolute value is used because for the soaring state it does not matter if the turning direction is always left or always right.
* Height first difference: Mean &rarr; normally distributed


## Model formulation

3D-time series
