# AD-Workflow - Aircraft Design Workflow

## Introduction
AD-Workflow is a MATLAB tool developed to guide users through the preliminary design phases of an aircraft. By defining the *Top Level Aircraft Requirements* (TLARs) and selecting reference aircraft, the tool produces a preliminary sizing of the aircraft and its wing, based on regression analysis of similar aircraft, following the methods outlined in [1].

This project was developed as an assignment for the Master's Degree course in *Aircraft Design*.

Requires MATLAB 2019b or newer.

---

## Basic Usage
To start, define the TLARs in the [`Geometry file`](tlars/Geometry.txt) and [`TLARs file`](tlars/TLARs.txt), specifying the design mission profile and some assumptions on the engine and drag characteristics of the aircraft. Lift coefficients for the various phases of the mission profile, altitudes, and non-standard conditions can be adjusted by editing the corresponding variables in the [`main script's`](RaccoltaDati_main.m) input section.

This design phase relies on statistical data from similar aircraft, as reported in [1]. Aircraft data must be provided in text files, following the format shown in the [`example A321 text file`](statistical_data/A321LR_LEAP.txt). The script derives a polar drag model based on statistical laws adapted from [1] and generates a T/W - W/S carpet plot, showing the boundaries for each requirement. You can then choose design values of T/W and W/S, and the program will recalculate the boundaries accordingly.

## Wing Design
Once the final values of T/W and W/S are set, you can define the planform parameters and profile sections. The tool lets you select the wing's planform parameters — leading edge sweep angle, taper ratio, kink position, and so on — as well as the root, kink, and tip section properties, such as profile thickness and aerodynamic characteristics, to estimate the wing's properties. These estimates are used to check whether the aircraft may encounter compressibility issues during cruise, using the theory developed in [2]. The stall path can be assessed by comparing the sectional lift coefficients to the maximum lift coefficient, obtained through an implementation of the Weissinger method taken from [3].

---

# References
[1] Roskam, Jan (1985). *Airplane Design Part I*.
[2] Ciornei, Simona (2005). *Mach number, relative thickness, sweep and lift coefficient of the wing - An empirical investigation of parameters and equations*. Tech. rep. Hamburg University of Applied Science.
[3] DeYoung, John and Harper, Charles W. (1948). "Theoretical symmetric span loading at subsonic speeds for wings having arbitrary plan form". In: *NACA TR-921*.