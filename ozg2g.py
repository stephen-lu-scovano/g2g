"""
A tool to convert geographic to / from grid coordinates by the method
from Krueger &lambda;-series equations (Krueger, 1912).
"""
import math


"""Constants and Parameters

Copy & paste the values required from the table below into the
vairables, or enter your own values

Coordinate Set      Ellipsoid   Semi-major Axis     Inverse Flattening
----------------------------------------------------------------------
GDA2020/MGA2020     GRS80       6,378,137.0 m       298.257222101
GDA94/MGA94         GRS80       6,378,137.0 m       298.257222101
AGD/AMG             ANS         6,378,160.0 m       298.25
ANG                 CLARKE 1858 6,975,449.335 yd    294.26

"""
ELLIPSOID_DEFINITION = "GRS80"
SEMI_MAJOR_AXIS = 6378137.0
INVERSE_FLATTENING = 298.257222101

FLATTENING = 1/INVERSE_FLATTENING
SEMI_MINOR_AXIS = SEMI_MAJOR_AXIS * (1 - FLATTENING)
ECCENTRICITY = 2 * FLATTENING - FLATTENING * FLATTENING
E = math.sqrt(ECCENTRICITY)
SECOND_ECCENTRICITY = ECCENTRICITY / (1 - ECCENTRICITY)
E_PRIME = math.sqrt(SECOND_ECCENTRICITY)

N = (SEMI_MAJOR_AXIS - SEMI_MINOR_AXIS) / (SEMI_MAJOR_AXIS + SEMI_MINOR_AXIS)
N2 = N * N
N3 = N2 * N
N4 = N3 * N
G = SEMI_MAJOR_AXIS * (1-N) * (1-N2) * (1 + 9*N2/4 + 225*N4/64) * math.pi / 180

TM_DEFINITION = 'GDA2020/MGA2020'
FALSE_EASTING = 500000.0
FALSE_NORTHING = 10000000.0
# Central Scale factor
KAPPA_0 = 0.9996
ZONE_WIDTH = 6
# Longitude of the central meridian of zone 1 (degrees)
LNG_ZONE_1 = -177
# Longitude of western edge of zone zero (degrees)
LNG_ZONE_0 = -186
# Central meridian of zone zero (degrees)
CENTRAL_ZONE_0 = -183


def geo2grid(lat, lng):
    """Convert Geographic to Grid

    Parameters
    ----------
    lat : float
        Latitude in real number
    lng : float
        Longitud in real number
    
    Returns
    -------
    tuple
        (easting, northing, zone, gamma, kappa)
        gamma is grid convergence,
        kappa is the point scale factor.

    """
    pass

def grid2geo(easting, northing, zone):
    """[summary]

    Parameters
    ----------
    easting : float
        Easting in real number
    northing : float
        Northing in real number
    zone : integer
        Zone number

    Returns
    -------
    tuple
        (latitude, longitude, gamma, kappa)
        latitude, longitude in decimal degree
        gamma is grid convergence,
        kappa is the point scale factor.

    """

    e_prime = easting - FALSE_EASTING
    n_prime = northing - FALSE_NORTHING

    sigma = n_prime / KAPPA_0 * math.pi / (G*180)
    fp_term_1 = sigma
    fp_term_2 = (3*N/2 - 27*N3/32) * math.sin(2*sigma)
    fp_term_3 = (21*N2/16 - 55*N4/32) * math.sin(4*sigma)
    fp_term_4 = 151 * N3 / 96 * math.sin(6*sigma)
    fp_term_5 = 1096 * N4 / 512 * math.sin(8*sigma)

    fp_lat_prime = fp_term_1 + fp_term_2 + fp_term_3 + fp_term_4 + fp_term_5
    sin_j_prime = math.sin(fp_lat_prime)
    sec_j_prime = 1/math.cos(fp_lat_prime)

    rho_prime = SEMI_MAJOR_AXIS * (1-ECCENTRICITY)/math.pow((1-ECCENTRICITY*sin_j_prime*sin_j_prime), 1.5)
    nu_prime = SEMI_MAJOR_AXIS / math.sqrt(1-ECCENTRICITY*sin_j_prime*sin_j_prime)

    x_1 = e_prime / KAPPA_0 / nu_prime
    x_2 = x_1*x_1
    x_3 = x_2*x_1

    y_1 = nu_prime / rho_prime
    y_2 = y_1*y_1
    y_3 = y_2*y_1
    y_4 = y_3*y_1

    t_prime_1 = math.tan(fp_lat_prime)
    t_prime_2 = t_prime_1 * t_prime_1
    t_prime_3 = t_prime_1 * t_prime_2
    t_prime_4 = t_prime_2 * t_prime_2

    foot_point_latitude = fp_lat_prime
    lat_term_1 = -t_prime_1 / (KAPPA_0*rho_prime) * x_1 * e_prime / 2
    lat_term_2 = t_prime_1 / (KAPPA_0*rho_prime) * x_3 * e_prime / 24 * (-4*y_2+9*y_1*(1-t_prime_2)+12*t_prime_2)
    lat_term_3 = -t_prime_1 / (KAPPA_0*rho_prime) * (x_2*x_3*e_prime/720) * (8*y_4*(11-24*t_prime_2)
                    - 12*y_3*(21-71*t_prime_2) + 15*y_2*(15-98*t_prime_2+15*t_prime_4)
                    + 180*y_1*(5*t_prime_2-3*t_prime_4) + 360*t_prime_4)
    lat_term_4 = t_prime_1 / (KAPPA_0*rho_prime) * x_3*x_3*x_1*e_prime/40320 * (
                    1385 + 3633*t_prime_2 + 4095*t_prime_4 + 1575*t_prime_3*t_prime_3)
    
    # Latitude
    lat = foot_point_latitude + lat_term_1 + lat_term_2 + lat_term_3 + lat_term_4

    lng_central_deg = zone * ZONE_WIDTH + LNG_ZONE_1 - ZONE_WIDTH
    lng_central_rad = lng_central_deg / 180 * math.pi
    lng_term_1 = sec_j_prime * x_1
    lng_term_2 = -sec_j_prime * x_3/6 * (y_1+2*t_prime_2)
    lng_term_3 = sec_j_prime * x_3*x_2/120 * (-4*y_3*(1-6*t_prime_2)
                    + y_2*(9-68*t_prime_2) + 72*y_1*t_prime_2 + 24*t_prime_4)
    lng_term_4 = -sec_j_prime * x_3*x_3*x_1/5040 * (61 + 662*t_prime_2 + 1320*t_prime_4 + 720*t_prime_4*t_prime_2)

    # Longitude
    lng = lng_central_rad + lng_term_1 + lng_term_2 + lng_term_3 + lng_term_4

    gamma_1 = -t_prime_1 * x_1
    gamma_2 = t_prime_1 * x_3 / 3 * (-2*y_2+3*y_1+t_prime_2)
    gamma_3 = -t_prime_1*x_3*x_2/15 * (y_4*(11-24*t_prime_2)
                - 3*y_3*(8-23*t_prime_2) + 5*y_2*(3-14*t_prime_2)
                + 30*y_1*t_prime_2 + 3*t_prime_4)
    gamma_4 = t_prime_1*x_3*x_3*x_1 * (17+77*t_prime_2+105*t_prime_4+45*t_prime_4*t_prime_2)

    # Grid Convergence
    gamma = gamma_1 + gamma_2 + gamma_3 + gamma_4

    i1 = e_prime/KAPPA_0*e_prime/KAPPA_0 / (rho_prime * nu_prime)
    i2 = i1 * i1
    i3 = i2 * i1

    ps_1 = 1+i1/2
    ps_2 = i2/24 * (4*y_1*(1-6*t_prime_2)-3*(1-16*t_prime_2)-24*t_prime_2/y_1)
    ps_3 = i3 / 720
    ps_sum = ps_1 + ps_2 + ps_3

    # Point Scale
    kappa = ps_sum * KAPPA_0

    return (lat*180/math.pi, lng*180/math.pi, gamma, kappa)

