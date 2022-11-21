# Notes for the report

## Data

* longitude and latitude: Converted to Step length/ speed (m/s) and Turning angle (1 Hz)
* height above mean sealevel: Calculate first difference (1 Hz)
* possible covariates:
** temperature
** landform (categorical variable: 6 Landforms + NA)
** heading might be useful
** solar time (ToD in decimal)
** day of observation
** elevation

Challenges of the data: 30 Second intervals of 1 Hz data then 15 Minutes (900 seconds) breaks
As basically no state switching was observed within intervals we decided to compute summary statistics for each interval.

