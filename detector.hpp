#ifndef _DETECTOR_HPP
#define _DETECTOR_HPP

/* Detector coordinates */
const double z1 /* m */ = 200;
const double z2 /* m */ = 205;
const double z3 /* m */ = 250;
const double z4 /* m */ = 255;

/* Detector cuts: detector is not present in a square area around the beam direction: [-x,+x]*[-y,+y]
 * This also avoid high pz pions, which would be very poorly reconstructed */
const double xcut /* m */ = 2;
const double ycut /* m */ = 2;

/* Detector smearing */
const double dx /* m */ =  0.001;
const double dy /* m */ =  0.001;

/* dp kick - applied at the center of the field region
 * Looks much more like a delta of E field at z_kick, but makes math muuuuch simpler */
const double z_drift /* m */ = z3-z2;
const double z_kick /* m */ = (z2+z3)/2;
const double p_kick /* GeV */ = 1e-4*z_drift;

#endif // _DETECTOR
