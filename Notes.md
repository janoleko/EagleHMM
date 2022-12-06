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

To understand the behaviour of the eagle we visualized the coordinates as well as the height movement patterns for each of ~ 250-300 intervals.

We came up with 3 to 4 clearly distinguishable behaviours which the HMM could maybe capture with 3 or 4 states:

* Soaring: The coordinate plots show circular patterns and the height plot shows upward movement.
* Gliding: The coordinate plot shows fast and directed movement while the height plot shows either downwards movement or no change in altitude. This could be captured in one or two states:
  * Either one state with directed movement and descending as well as some horizontal/ ascending movement,
  * or two states: One only for descending movement, the other one especially for directed movement ascending.
* Resting: The coordinate plot shows little movement with a rather random pattern while the height plot shows barely any changes in height.


### Summary statistics

* Step length/ speed (in m/s): mean() &rarr; Gamma distributed
* Turning angle: abs(mean())/pi &rarr; from von Mises distr. to beta distribution
  &rarr; Absolute value is used because for the soaring state it does not matter if the turning direction is always left or always right.
  * We deleted values where the interval only contained < 3 turning angles as there the mean is not an reliable estimate for the behaviour in the interval, which resulted in significant mass on 0.5 and 1.
* Height first difference: mean() &rarr; normally distributed


## Model formulation

3D-time series: Step length, transformed turning angle, height first difference
* assuming contemporaneous conditional independence
* The resulting component distributions are of the form: dgamma * dbeta * dnorm
* 3 State HMM: 
  * State 1: Soaring
  * State 2: Gliding/ flying
  * State 3: Resting
* 4 State HMM:
  * State 1: Soaring
  * State 2: Gliding
  * State 3: Resting
  * State 4: Flapping flight?
  
## Todo
* Likelihood umschreiben --> jeweils Schleife über Tage
* Was machen (day Variable ist kacke)?
* Heading für turning angle benutzen -> first difference (funktioniert nicht (warum?))
* scatter3D Farben checken
* landform mit aggregieren und anschauen wie die Beziehung zwischen Decoded states und landform ist: Grundsätzlich gute Idee für ersten Eindruck welche Covariaten interessant sein könnten (Multinomiale Regression)
* covariaten für transition probs
* gemeinsames HMM für mehrere Vögel?

