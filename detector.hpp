#ifndef _DETECTOR_HPP
#define _DETECTOR_HPP

/* Detector coordinates */
const double z1 /* m */ = 200;
const double z2 /* m */ = 201;
const double z3 /* m */ = 250;
const double z4 /* m */ = 251;

/* Detector smearing */
const double dx /* m */ =  0.01;
const double dy /* m */ =  0.01;

/* dp kick - applied at the center of the field region */
const double z_drift /* m */ = z3-z2;
const double z_kick /* m */ = (z2+z3)/2;
const double p_kick /* GeV */ = 1e-5*z_drift;

#endif // _DETECTOR
