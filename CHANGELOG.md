# Changelog

## 0.6.0 - 2023-09-27

### Bug fix

-   Removed support for Python 2
-   Changes to creating DER model instance using "DERModel" class.
    - Removed need to specify "modelType" when creating model instance. Model type is specified in config file now.
	- "derID" is now a required argument
-   Minor refactoring and code clean up.
-   Added Dockerfile for using high performance Julia ODE solvers library - Acknowledgement to Abel Siqueira (Netherlands eScience Center)

## 0.5.3 - 2021-09-13

### Bug fix

-   Removed pathlib module from setup.py

## 0.5.2 - 2021-09-09

### Bug fix

-   Added Numba to setup.py

## 0.5.1 - 2021-07-28

### Bug fix

-   Fixed bug due to logs folder not being found.

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
