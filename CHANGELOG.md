# Changelog

## 0.5.0 - 2021-07-23

### Features

-   Added current reference ramp rate limiter for all models.

### Bug fix

-   Fixed bug in HVRT logging.
-   Fixed bug that multiplied *m_limit* with *1e1* in the current control loop integral reset logic.
-   Changed controller gains for *ConstantVdc* models to make their response slightly underdamped. Previously the response was overdamped.

## 0.4.0 - 2021-06-02

### Features

-   Added minimum ride-through time before momentary cessation.
-   Added additional logging information during voltage anomalies.

### Bug fix

-   Minor bug fixes and code clean up
