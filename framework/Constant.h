#pragma once

const double mass_pip = 0.13957061;
const double mass_Kp = 0.49367700;
const double mass_eta = 0.547862;
const double mass_Bs = 5.36688;
const double mass_Du = 1.86484;

// Particle radius
const double Rrad = 4.0;
const double Srad = 4.0;

// Selection criteria
const double kpi_min = mass_Kp + mass_pip;
const double kpi_max = mass_Kp + mass_eta;

const double time_min = 0.0;
const double time_max = 15.0;
