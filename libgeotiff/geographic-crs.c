#include <stdio.h>
#include <math.h>

#define DEG2RAD(d) ((d) * M_PI / 180.0)

/* Convert geodetic (lat, lon, alt) to ECEF (x, y, z) */
static void geodetic_to_ecef(double lat_deg, double lon_deg, double alt,
                              double *x, double *y, double *z)
{
    double lat = DEG2RAD(lat_deg);
    double lon = DEG2RAD(lon_deg);
    double sin_lat = sin(lat);
    double cos_lat = cos(lat);

    /* Radius of curvature in the prime vertical */
    double N = WGS84_A / sqrt(1.0 - WGS84_E2 * sin_lat * sin_lat);

    *x = (N + alt) * cos_lat * cos(lon);
    *y = (N + alt) * cos_lat * sin(lon);
    *z = (N * (1.0 - WGS84_E2) + alt) * sin_lat;
}

/* Convert ECEF (x, y, z) to local ENU relative to a reference
 * geodetic origin (lat0, lon0, alt0). Result is in meters. */
static void ecef_to_enu(double x, double y, double z,
                         double lat0_deg, double lon0_deg, double alt0,
                         double *east, double *north, double *up)
{
    double x0, y0, z0;
    geodetic_to_ecef(lat0_deg, lon0_deg, alt0, &x0, &y0, &z0);

    double dx = x - x0;
    double dy = y - y0;
    double dz = z - z0;

    double lat0 = DEG2RAD(lat0_deg);
    double lon0 = DEG2RAD(lon0_deg);

    double sin_lat0 = sin(lat0), cos_lat0 = cos(lat0);
    double sin_lon0 = sin(lon0), cos_lon0 = cos(lon0);

    *east  = -sin_lon0 * dx + cos_lon0 * dy;
    *north = -sin_lat0 * cos_lon0 * dx - sin_lat0 * sin_lon0 * dy + cos_lat0 * dz;
    *up    =  cos_lat0 * cos_lon0 * dx + cos_lat0 * sin_lon0 * dy + sin_lat0 * dz;
}

/* Convenience wrapper: geodetic point -> local ENU relative to origin */
static void geodetic_to_local_enu(double lat_deg, double lon_deg, double alt,
                                   double lat0_deg, double lon0_deg, double alt0,
                                   double *east, double *north, double *up)
{
    double x, y, z;
    geodetic_to_ecef(lat_deg, lon_deg, alt, &x, &y, &z);
    ecef_to_enu(x, y, z, lat0_deg, lon0_deg, alt0, east, north, up);
}

int main(void)
{
    /* Local origin (reference point) */
    double lat0 = 48.8566, lon0 = 2.3522, alt0 = 0.0; /* e.g. Paris */

    /* Point to convert */
    double lat = 48.8580, lon = 2.3540, alt = 0.0;

    double east, north, up;
    geodetic_to_local_enu(lat, lon, alt, lat0, lon0, alt0, &east, &north, &up);

    printf("East:  %.3f m\n", east);
    printf("North: %.3f m\n", north);
    printf("Up:    %.3f m\n", up);

    return 0;
}

/* Rotate local ENU delta back into ECEF delta, then add origin's ECEF */
static void enu_to_ecef(double east, double north, double up,
                         double lat0_deg, double lon0_deg, double alt0,
                         double *x, double *y, double *z)
{
    double x0, y0, z0;
    geodetic_to_ecef(lat0_deg, lon0_deg, alt0, &x0, &y0, &z0);

    double lat0 = DEG2RAD(lat0_deg);
    double lon0 = DEG2RAD(lon0_deg);
    double sin_lat0 = sin(lat0), cos_lat0 = cos(lat0);
    double sin_lon0 = sin(lon0), cos_lon0 = cos(lon0);

    /* Transpose of the ECEF->ENU rotation matrix */
    double dx = -sin_lon0 * east - sin_lat0 * cos_lon0 * north + cos_lat0 * cos_lon0 * up;
    double dy =  cos_lon0 * east - sin_lat0 * sin_lon0 * north + cos_lat0 * sin_lon0 * up;
    double dz =  cos_lat0 * north + sin_lat0 * up;

    *x = x0 + dx;
    *y = y0 + dy;
    *z = z0 + dz;
}

/* ECEF -> geodetic, Bowring's method (closed-form, no iteration needed) */
static void ecef_to_geodetic(double x, double y, double z,
                              double *lat_deg, double *lon_deg, double *alt)
{
    double a = WGS84_A;
    double e2 = WGS84_E2;
    double b = a * sqrt(1.0 - e2);
    double ep2 = (a * a - b * b) / (b * b); /* second eccentricity squared */

    double p = sqrt(x * x + y * y);
    double theta = atan2(z * a, p * b);

    double sin_t = sin(theta), cos_t = cos(theta);

    double lat = atan2(z + ep2 * b * sin_t * sin_t * sin_t,
                        p - e2 * a * cos_t * cos_t * cos_t);
    double lon = atan2(y, x);

    double sin_lat = sin(lat);
    double N = a / sqrt(1.0 - e2 * sin_lat * sin_lat);
    double h = p / cos(lat) - N;

    *lat_deg = lat * 180.0 / M_PI;
    *lon_deg = lon * 180.0 / M_PI;
    *alt = h;
}

/* Convenience wrapper: local ENU point -> geodetic (lat, lon, alt) */
static void local_enu_to_geodetic(double east, double north, double up,
                                   double lat0_deg, double lon0_deg, double alt0,
                                   double *lat_deg, double *lon_deg, double *alt)
{
    double x, y, z;
    enu_to_ecef(east, north, up, lat0_deg, lon0_deg, alt0, &x, &y, &z);
    ecef_to_geodetic(x, y, z, lat_deg, lon_deg, alt);
}

for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
        double east  = X0 + i * step;
        double north = Y0 + j * step;
        local_enu_to_geodetic(east, north, 0.0, lat0, lon0, alt0, &lat[i][j], &lon[i][j], &alt[i][j]);
    }
}

dlat/dnorth = 1 / M(lat0)              // M = meridian radius of curvature
dlon/deast  = 1 / (N(lat0) * cos(lat0)) // N = prime-vertical radius of curvature
dlat/deast  = 0   (to first order)
dlon/dnorth = 0   (to first order)

static void compute_local_jacobian(double lat0_deg,
                                    double *dlat_dnorth, double *dlon_deast)
{
    double lat0 = DEG2RAD(lat0_deg);
    double sin_lat0 = sin(lat0);
    double denom = sqrt(1.0 - WGS84_E2 * sin_lat0 * sin_lat0);

    double M = WGS84_A * (1 - WGS84_E2) / (denom * denom * denom); /* meridian radius */
    double N = WGS84_A / denom;                                    /* prime vertical radius */

    *dlat_dnorth = 1.0 / M;
    *dlon_deast  = 1.0 / (N * cos(lat0));
}

/* Fill the whole grid using one linearization at the origin */
static void fill_grid_linear(double X0, double Y0, double step, int N,
                              double lat0_deg, double lon0_deg,
                              double lat_out[][N], double lon_out[][N])
{
    double dlat_dnorth, dlon_deast;
    compute_local_jacobian(lat0_deg, &dlat_dnorth, &dlon_deast);

    for (int i = 0; i < N; i++) {
        double east = X0 + i * step;
        for (int j = 0; j < N; j++) {
            double north = Y0 + j * step;
            lat_out[i][j] = lat0_deg + (north * dlat_dnorth) * 180.0 / M_PI;
            lon_out[i][j] = lon0_deg + (east  * dlon_deast)  * 180.0 / M_PI;
        }
    }
}

double eps = 1.0; /* 1 meter probe */
double lat_e, lon_e, lat_n, lon_n, dummy_alt;
local_enu_to_geodetic(eps, 0, 0, lat0, lon0, alt0, &lat_e, &lon_e, &dummy_alt);
local_enu_to_geodetic(0, eps, 0, lat0, lon0, alt0, &lat_n, &lon_n, &dummy_alt);

double dlon_deast  = (lon_e - lon0) / eps;
double dlat_dnorth = (lat_n - lat0) / eps;
