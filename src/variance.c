#include "libairs.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/* Number of latitudes for threshold tables. */
#define NLAT_GW 19
#define NLAT_SURF 6
#define NLAT_TROP 73

/* Number of months for threshold tables. */
#define NMON 12

/* Maximum number of longitudes. */
#define NX 3600

/* Maximum number of latitudes. */
#define NY 1800

/* ------------------------------------------------------------
   Global variables...
   ------------------------------------------------------------ */

/* Latitudes for gravity wave variance thresholds. */
static double t_gw_lat[NLAT_GW]
  = { -90, -80, -70, -60, -50, -40, -30, -20, -10, 0,
  10, 20, 30, 40, 50, 60, 70, 80, 90
};

/* Gravity wave variance thresholds (ascending orbits). */
static double t_gw_asc[NMON][NLAT_GW]
  = { {0.00387, 0.00422, 0.00633, 0.0124, 0.0216, 0.0324,
       0.0553, 0.0791, 0.0501, 0.0136, 0.0134, 0.0151,
       0.0522, 0.321, 0.697, 0.776, 0.696, 0.764, 0.771},
{0.00913, 0.00942, 0.00867, 0.00897, 0.0112, 0.0168,
 0.0314, 0.0484, 0.032, 0.0128, 0.0122, 0.0134,
 0.0382, 0.124, 0.345, 0.404, 0.545, 1.16, 1.18},
{0.0845, 0.0664, 0.0384, 0.0227, 0.0147, 0.0118,
 0.0141, 0.0184, 0.0162, 0.0123, 0.0124, 0.0124,
 0.0159, 0.0509, 0.085, 0.103, 0.188, 0.367, 0.529},
{0.265, 0.297, 0.216, 0.106, 0.0666, 0.0299,
 0.0169, 0.0129, 0.0116, 0.012, 0.0135, 0.0141,
 0.0134, 0.0137, 0.017, 0.0268, 0.0259, 0.0319, 0.0323},
{0.326, 0.44, 0.628, 0.567, 0.434, 0.235,
 0.0601, 0.0214, 0.0132, 0.0113, 0.0144, 0.0185,
 0.0179, 0.0142, 0.0116, 0.00945, 0.00865, 0.00918, 0.00878},
{0.537, 0.73, 1.39, 1.75, 1.35, 0.528,
 0.188, 0.0311, 0.0133, 0.0124, 0.0205, 0.0313,
 0.0297, 0.0216, 0.0166, 0.0131, 0.00983, 0.00606, 0.0049},
{0.382, 1.15, 1.57, 2.13, 1.66, 0.851,
 0.126, 0.0204, 0.0133, 0.0135, 0.0281, 0.0385,
 0.0375, 0.0312, 0.0223, 0.0143, 0.00949, 0.0061, 0.00493},
{0.226, 0.697, 1.68, 1.56, 1.14, 0.496,
 0.0616, 0.0143, 0.0126, 0.013, 0.0216, 0.0252,
 0.0241, 0.0206, 0.0152, 0.0106, 0.00976, 0.0105, 0.00998},
{0.236, 0.489, 0.648, 0.553, 0.524, 0.21,
 0.033, 0.0129, 0.0116, 0.0129, 0.0163, 0.0165,
 0.0153, 0.014, 0.0141, 0.0185, 0.0301, 0.0591, 0.0745},
{0.046, 0.082, 0.112, 0.0806, 0.0516, 0.0469,
 0.0225, 0.0139, 0.0127, 0.0121, 0.0125, 0.0138,
 0.0176, 0.0357, 0.0563, 0.062, 0.133, 0.327, 0.3},
{0.00669, 0.00867, 0.0117, 0.0117, 0.014, 0.015,
 0.0203, 0.0213, 0.0144, 0.0116, 0.0124, 0.0179,
 0.0574, 0.185, 0.346, 0.442, 0.54, 0.669, 0.664},
{0.00355, 0.00381, 0.00658, 0.0125, 0.0217, 0.0304,
 0.0424, 0.0515, 0.0315, 0.0139, 0.0137, 0.0161,
 0.0582, 0.306, 0.999, 1.2, 1.14, 0.621, 0.448}
};

/* Gravity wave variance thresholds (descending orbits). */
static double t_gw_dsc[NMON][NLAT_GW]
  = { {0.00383, 0.00458, 0.00866, 0.019, 0.0348, 0.0598,
       0.144, 0.234, 0.135, 0.0373, 0.0325, 0.0377,
       0.0858, 0.497, 1.4, 1.32, 0.808, 0.771, 0.773},
{0.00999, 0.0123, 0.0141, 0.0148, 0.0177, 0.0286,
 0.0626, 0.102, 0.0717, 0.0302, 0.0261, 0.03,
 0.086, 0.268, 0.631, 0.716, 1.17, 1.24, 1.21},
{0.103, 0.096, 0.0715, 0.0535, 0.0343, 0.0245,
 0.025, 0.0315, 0.0303, 0.0233, 0.023, 0.0257,
 0.0353, 0.118, 0.197, 0.359, 0.541, 0.585, 0.586},
{0.272, 0.293, 0.276, 0.226, 0.146, 0.0689,
 0.0373, 0.0245, 0.0232, 0.0232, 0.0224, 0.0217,
 0.0242, 0.031, 0.0441, 0.0664, 0.0623, 0.053, 0.0361},
{0.331, 0.44, 0.641, 0.868, 0.824, 0.47,
 0.115, 0.0444, 0.0269, 0.0223, 0.0274, 0.0332,
 0.0273, 0.023, 0.0191, 0.0172, 0.0138, 0.0107, 0.00894},
{0.554, 0.716, 1.31, 2.29, 2.43, 1.05,
 0.41, 0.0651, 0.0269, 0.0257, 0.0447, 0.0622,
 0.0497, 0.0357, 0.0258, 0.0182, 0.0117, 0.00697, 0.00502},
{0.427, 0.905, 1.44, 2.78, 2.76, 1.52,
 0.278, 0.041, 0.0279, 0.0296, 0.0629, 0.0818,
 0.0758, 0.0534, 0.0356, 0.0227, 0.012, 0.00692, 0.00513},
{0.245, 0.74, 1.88, 2.32, 1.89, 0.883,
 0.122, 0.0292, 0.0264, 0.0289, 0.0516, 0.059,
 0.0495, 0.0373, 0.0268, 0.0185, 0.0163, 0.0131, 0.0103},
{0.272, 0.551, 0.812, 0.844, 0.852, 0.486,
 0.0842, 0.0269, 0.0225, 0.0239, 0.0322, 0.0324,
 0.0307, 0.0304, 0.035, 0.0484, 0.0692, 0.0956, 0.0948},
{0.0644, 0.125, 0.177, 0.135, 0.0922, 0.0899,
 0.0524, 0.0249, 0.0214, 0.0218, 0.0251, 0.0293,
 0.0403, 0.0903, 0.168, 0.246, 0.358, 0.378, 0.288},
{0.00676, 0.00923, 0.0148, 0.0195, 0.0261, 0.0286,
 0.0302, 0.0343, 0.0298, 0.024, 0.0252, 0.0403,
 0.131, 0.448, 0.681, 0.923, 0.839, 0.684, 0.629},
{0.00347, 0.00412, 0.00995, 0.0221, 0.0363, 0.0531,
 0.104, 0.168, 0.112, 0.0365, 0.0335, 0.0382,
 0.128, 0.563, 1.62, 1.87, 1.47, 0.652, 0.408}
};

/* Latitudes for zonal mean tropopause temperatures. */
static double t_trop_lat[NLAT_TROP]
  = { 90, 87.5, 85, 82.5, 80, 77.5, 75, 72.5, 70, 67.5, 65, 62.5, 60,
  57.5, 55, 52.5, 50, 47.5, 45, 42.5, 40, 37.5, 35, 32.5, 30, 27.5,
  25, 22.5, 20, 17.5, 15, 12.5, 10, 7.5, 5, 2.5, 0, -2.5, -5, -7.5,
  -10, -12.5, -15, -17.5, -20, -22.5, -25, -27.5, -30, -32.5, -35,
  -37.5, -40, -42.5, -45, -47.5, -50, -52.5, -55, -57.5, -60, -62.5,
  -65, -67.5, -70, -72.5, -75, -77.5, -80, -82.5, -85, -87.5, -90
};

/* Zonal mean tropopause temperatures. */
static double t_trop[NMON][NLAT_TROP]
  = { {211.152, 211.237, 211.434, 211.549, 211.614, 211.776, 211.974,
       212.234, 212.489, 212.808, 213.251, 213.692, 214.193, 214.591,
       214.985, 215.327, 215.658, 215.956, 216.236, 216.446, 216.738,
       216.836, 216.032, 213.607, 209.281, 205, 201.518, 198.969,
       197.123, 195.869, 195.001, 194.409, 193.985, 193.734, 193.617,
       193.573, 193.6, 193.642, 193.707, 193.856, 194.131, 194.558,
       195.121, 195.907, 196.91, 198.192, 199.744, 201.583, 203.672,
       206.012, 208.542, 211.135, 213.681, 216.085, 218.317, 220.329,
       222.071, 223.508, 224.612, 225.357, 225.761, 225.863, 225.657,
       225.287, 224.813, 224.571, 224.385, 224.3, 224.257, 224.173,
       223.786, 222.713, 222.11},
{212.593, 212.621, 212.801, 212.888, 212.912, 213.054, 213.245,
 213.512, 213.726, 213.962, 214.259, 214.508, 214.823, 215.037,
 215.297, 215.545, 215.808, 216.063, 216.323, 216.539, 216.867,
 217.051, 216.532, 214.512, 210.371, 205.658, 201.758, 198.937,
 197.047, 195.817, 194.96, 194.386, 193.993, 193.771, 193.673,
 193.635, 193.658, 193.691, 193.744, 193.872, 194.126, 194.54,
 195.085, 195.847, 196.8, 198.013, 199.489, 201.261, 203.298,
 205.596, 208.082, 210.628, 213.156, 215.563, 217.822, 219.903,
 221.745, 223.311, 224.566, 225.451, 225.947, 226.079, 225.849,
 225.406, 224.889, 224.643, 224.431, 224.246, 224.079, 223.884,
 223.42, 222.402, 221.871},
{215.529, 215.491, 215.539, 215.621, 215.691, 215.808, 215.847,
 215.881, 215.878, 215.907, 216.02, 216.113, 216.297, 216.342,
 216.38, 216.369, 216.342, 216.284, 216.185, 215.989, 215.855,
 215.626, 215.023, 213.432, 209.979, 205.886, 202.212, 199.414,
 197.488, 196.216, 195.327, 194.732, 194.347, 194.158, 194.095,
 194.079, 194.116, 194.154, 194.195, 194.302, 194.534, 194.922,
 195.461, 196.253, 197.288, 198.644, 200.309, 202.293, 204.553,
 207.033, 209.538, 211.911, 214.016, 215.862, 217.572, 219.179,
 220.655, 221.959, 223.052, 223.867, 224.344, 224.451, 224.179,
 223.706, 223.163, 222.876, 222.613, 222.385, 222.154, 221.842,
 221.304, 220.402, 220.06},
{219.921, 219.916, 219.99, 219.989, 219.916, 219.867, 219.73,
 219.522, 219.16, 218.765, 218.448, 218.144, 217.99, 217.756,
 217.553, 217.311, 217.025, 216.684, 216.241, 215.649, 215.05,
 214.302, 213.219, 211.496, 208.729, 205.649, 202.594, 200.066,
 198.144, 196.733, 195.687, 194.991, 194.586, 194.429, 194.418,
 194.443, 194.492, 194.534, 194.59, 194.718, 194.997, 195.481,
 196.165, 197.159, 198.462, 200.142, 202.154, 204.533, 207.208,
 209.848, 212.088, 213.845, 215.222, 216.348, 217.384, 218.383,
 219.313, 220.131, 220.799, 221.271, 221.479, 221.405, 221.012,
 220.4, 219.702, 219.227, 218.827, 218.434, 217.977, 217.477,
 216.783, 215.974, 215.707},
{225.363, 225.255, 225.064, 224.745, 224.351, 224, 223.551,
 222.966, 222.195, 221.435, 220.802, 220.245, 219.871, 219.424,
 218.99, 218.529, 218.013, 217.445, 216.76, 215.859, 214.723,
 213.049, 211.032, 208.767, 206.449, 204.302, 202.113, 200.187,
 198.501, 197.153, 196.117, 195.441, 195.121, 195.073, 195.146,
 195.212, 195.261, 195.288, 195.343, 195.485, 195.772, 196.284,
 197.018, 198.125, 199.624, 201.604, 204.073, 207.036, 210.193,
 212.853, 214.611, 215.635, 216.287, 216.801, 217.284, 217.716,
 218.057, 218.253, 218.282, 218.115, 217.729, 217.15, 216.376,
 215.449, 214.428, 213.574, 212.847, 212.281, 211.718, 211.211,
 210.616, 210.112, 210.056},
{228.431, 228.261, 227.966, 227.457, 226.812, 226.208, 225.518,
 224.71, 223.701, 222.762, 222.045, 221.486, 221.142, 220.761,
 220.361, 219.896, 219.34, 218.646, 217.626, 215.983, 213.624,
 210.817, 208.017, 205.73, 203.8, 202.363, 200.96, 199.778,
 198.695, 197.845, 197.166, 196.743, 196.6, 196.66, 196.809,
 196.925, 196.985, 196.996, 197.033, 197.135, 197.335, 197.754,
 198.367, 199.335, 200.693, 202.564, 205.001, 208.084, 211.473,
 214.407, 216.208, 217.018, 217.314, 217.394, 217.371, 217.234,
 216.961, 216.517, 215.878, 215.027, 213.952, 212.697, 211.274,
 209.736, 208.172, 206.872, 205.84, 205.093, 204.32, 203.816,
 203.55, 203.49, 203.606},
{229.01, 228.807, 228.45, 227.839, 227.084, 226.377, 225.589,
 224.712, 223.665, 222.724, 222.058, 221.658, 221.519, 221.376,
 221.136, 220.673, 219.926, 218.742, 216.744, 214.028, 210.994,
 208.374, 206.131, 204.563, 203.251, 202.328, 201.313, 200.411,
 199.531, 198.876, 198.356, 198.104, 198.088, 198.21, 198.385,
 198.502, 198.57, 198.601, 198.652, 198.731, 198.869, 199.207,
 199.737, 200.595, 201.802, 203.491, 205.771, 208.765, 212.241,
 215.403, 217.439, 218.251, 218.297, 217.988, 217.533, 216.941,
 216.161, 215.154, 213.887, 212.35, 210.525, 208.481, 206.287,
 204.068, 202.033, 200.405, 199.106, 198.225, 197.435, 197.02,
 197.133, 197.527, 197.808},
{226.525, 226.354, 225.996, 225.433, 224.842, 224.358, 223.818,
 223.202, 222.426, 221.723, 221.266, 220.98, 220.893, 220.707,
 220.392, 219.928, 219.182, 218.015, 216.051, 213.399, 210.617,
 208.318, 206.311, 204.838, 203.515, 202.527, 201.397, 200.423,
 199.494, 198.848, 198.385, 198.212, 198.294, 198.49, 198.707,
 198.853, 198.933, 198.967, 199.01, 199.079, 199.207, 199.537,
 200.081, 200.968, 202.215, 203.946, 206.254, 209.291, 212.876,
 216.262, 218.487, 219.387, 219.436, 219.048, 218.405, 217.527,
 216.372, 214.919, 213.152, 211.096, 208.767, 206.247, 203.609,
 201.029, 198.763, 196.961, 195.578, 194.635, 193.923, 193.54,
 193.632, 193.944, 193.912},
{223.293, 223.158, 222.945, 222.571, 222.126, 221.749, 221.362,
 220.946, 220.404, 219.946, 219.704, 219.599, 219.611, 219.429,
 219.124, 218.702, 218.063, 217.157, 215.827, 213.879, 211.352,
 208.833, 206.504, 204.728, 203.168, 201.992, 200.735, 199.74,
 198.833, 198.213, 197.801, 197.661, 197.765, 197.963, 198.182,
 198.336, 198.42, 198.456, 198.505, 198.609, 198.794, 199.19,
 199.796, 200.758, 202.089, 203.915, 206.262, 209.295, 212.807,
 216.083, 218.329, 219.47, 219.877, 219.846, 219.507, 218.85,
 217.84, 216.448, 214.652, 212.509, 210.083, 207.534, 204.982,
 202.596, 200.463, 198.769, 197.441, 196.546, 195.902, 195.472,
 195.193, 195.066, 195.006},
{219.564, 219.492, 219.415, 219.191, 218.926, 218.801, 218.691,
 218.561, 218.298, 218.06, 217.982, 217.956, 218.038, 217.954,
 217.81, 217.532, 217.08, 216.439, 215.549, 214.31, 212.725,
 210.573, 208.019, 205.585, 203.459, 201.779, 200.162, 198.879,
 197.771, 196.987, 196.459, 196.19, 196.172, 196.274, 196.435,
 196.544, 196.601, 196.644, 196.727, 196.904, 197.184, 197.696,
 198.42, 199.497, 200.934, 202.825, 205.151, 208.005, 211.279,
 214.441, 216.87, 218.493, 219.498, 220.072, 220.353, 220.336,
 219.991, 219.271, 218.142, 216.636, 214.804, 212.776, 210.636,
 208.535, 206.516, 204.825, 203.383, 202.281, 201.365, 200.561,
 199.896, 199.415, 199.382},
{215.926, 215.884, 215.897, 215.814, 215.689, 215.692, 215.707,
 215.767, 215.815, 215.92, 216.138, 216.327, 216.588, 216.668,
 216.664, 216.553, 216.373, 216.112, 215.711, 215.025, 214.106,
 212.596, 210.346, 207.503, 204.604, 202.251, 200.231, 198.607,
 197.228, 196.174, 195.382, 194.87, 194.61, 194.54, 194.579,
 194.615, 194.66, 194.709, 194.82, 195.074, 195.487, 196.103,
 196.904, 198.01, 199.43, 201.246, 203.431, 206.007, 208.905,
 211.81, 214.34, 216.36, 217.918, 219.141, 220.159, 220.965,
 221.514, 221.754, 221.637, 221.135, 220.226, 218.986, 217.475,
 215.879, 214.251, 212.918, 211.84, 211.026, 210.288, 209.553,
 208.791, 208.132, 208.053},
{212.893, 212.911, 213.03, 213.109, 213.224, 213.453, 213.653,
 213.836, 213.98, 214.166, 214.481, 214.787, 215.179, 215.435,
 215.688, 215.908, 216.084, 216.217, 216.262, 216.123, 215.819,
 214.977, 213.173, 210.214, 206.619, 203.437, 200.836, 198.843,
 197.271, 196.078, 195.164, 194.509, 194.057, 193.82, 193.742,
 193.723, 193.762, 193.813, 193.903, 194.121, 194.49, 195.016,
 195.698, 196.627, 197.82, 199.359, 201.204, 203.355, 205.78,
 208.414, 211.057, 213.521, 215.662, 217.504, 219.133, 220.544,
 221.723, 222.631, 223.274, 223.649, 223.737, 223.547, 223.053,
 222.357, 221.52, 220.948, 220.527, 220.247, 220.013, 219.726,
 219.273, 218.506, 218.144}
};

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  static wave_t *wave;

  static FILE *in, *out;

  static char pertname[LEN], set[LEN];

  static double bt[NX][NY], bt_8mu[NX][NY], bt_8mu_min[NX][NY],
    bt_8mu_max[NX][NY], dt[NX][NY], mtime[NX][NY], glat[NY], glon[NX],
    fdc[NX][NY], fwg[NX][NY], fgw[NX][NY], fcw[NX][NY],
    mean[NX][NY], min[NX][NY], max[NX][NY], var[NX][NY],
    t_dc, t_gw, dt_trop, dc_hlat = 25, dc_tlim = 250, dt230,
    nesr, gauss_fwhm, var_dh, nu, lon0, lon1, lat0, lat1,
    thresh_dc, thresh_gw, lt, help[NX * NY];

  static int asc, ix, iy, nx, ny, iarg, n[NX][NY],
    ndc[NX][NY], ngw[NX][NY], ncw[NX][NY], nwg[NX][NY],
    det_gw, det_cw, det_dc, det_wg, ilat, imon, nmin = 10,
    bg_poly_x, bg_poly_y, bg_smooth_x, bg_smooth_y,
    itrack, itrack2, ixtrack, ixtrack2, iradius = 30, output, ncid, varid,
    minid, maxid, lonid, latid, npid, dimid[10], help2[NX * NY];

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <var.tab> <pert1.nc> [<pert2.nc> ...]");

  /* Get control parameters... */
  scan_ctl(argc, argv, "SET", -1, "full", set);
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  nx = (int) scan_ctl(argc, argv, "NX", -1, "360", NULL);
  lon0 = scan_ctl(argc, argv, "LON0", -1, "-180", NULL);
  lon1 = scan_ctl(argc, argv, "LON1", -1, "180", NULL);
  ny = (int) scan_ctl(argc, argv, "NY", -1, "180", NULL);
  lat0 = scan_ctl(argc, argv, "LAT0", -1, "-90", NULL);
  lat1 = scan_ctl(argc, argv, "LAT1", -1, "90", NULL);
  bg_poly_x = (int) scan_ctl(argc, argv, "BG_POLY_X", -1, "0", NULL);
  bg_poly_y = (int) scan_ctl(argc, argv, "BG_POLY_Y", -1, "0", NULL);
  bg_smooth_x = (int) scan_ctl(argc, argv, "BG_SMOOTH_X", -1, "0", NULL);
  bg_smooth_y = (int) scan_ctl(argc, argv, "BG_SMOOTH_Y", -1, "0", NULL);
  gauss_fwhm = scan_ctl(argc, argv, "GAUSS_FWHM", -1, "0", NULL);
  var_dh = scan_ctl(argc, argv, "VAR_DH", -1, "0", NULL);
  thresh_gw = scan_ctl(argc, argv, "THRESH_GW", -1, "-999", NULL);
  thresh_dc = scan_ctl(argc, argv, "THRESH_DC", -1, "-999", NULL);
  dt_trop = scan_ctl(argc, argv, "DT_TROP", -1, "0", NULL);
  dt230 = scan_ctl(argc, argv, "DT230", -1, "0.16", NULL);
  nu = scan_ctl(argc, argv, "NU", -1, "2345.0", NULL);
  output = (int) scan_ctl(argc, argv, "OUTPUT", -1, "1", NULL);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Check grid dimensions... */
  if (nx < 1 || nx > NX)
    ERRMSG("Set 1 <= NX <= MAX!");
  if (ny < 1 || ny > NY)
    ERRMSG("Set 1 <= NY <= MAX!");

  /* Loop over perturbation files... */
  for (iarg = 3; iarg < argc; iarg++) {

    /* Read perturbation data... */
    if (!(in = fopen(argv[iarg], "r")))
      continue;
    else {
      fclose(in);
      read_pert(argv[iarg], pertname, pert);
    }

    /* Recalculate background and perturbations... */
    if (bg_poly_x > 0 || bg_poly_y > 0 ||
	bg_smooth_x > 0 || bg_smooth_y > 0 || gauss_fwhm > 0 || var_dh > 0) {

      /* Allocate... */
      ALLOC(wave, wave_t, 1);

      /* Convert to wave analysis struct... */
      pert2wave(pert, wave, 0, pert->ntrack - 1, 0, pert->nxtrack - 1);

      /* Estimate background... */
      background_poly(wave, bg_poly_x, bg_poly_y);
      background_smooth(wave, bg_smooth_x, bg_smooth_y);

      /* Gaussian filter... */
      gauss(wave, gauss_fwhm);

      /* Compute variance... */
      variance(wave, var_dh);

      /* Copy data... */
      for (ix = 0; ix < wave->nx; ix++)
	for (iy = 0; iy < wave->ny; iy++) {
	  pert->pt[iy][ix] = wave->pt[ix][iy];
	  pert->var[iy][ix] = wave->var[ix][iy];
	}

      /* Free... */
      free(wave);
    }

    /* Detection... */
    for (itrack = 0; itrack < pert->ntrack; itrack++)
      for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {

	/* Check data... */
	if (pert->time[itrack][ixtrack] < 0
	    || pert->lon[itrack][ixtrack] < -180
	    || pert->lon[itrack][ixtrack] > 180
	    || pert->lat[itrack][ixtrack] < -90
	    || pert->lat[itrack][ixtrack] > 90
	    || pert->pt[itrack][ixtrack] < -100
	    || pert->pt[itrack][ixtrack] > 100
	    || !gsl_finite(pert->bt[itrack][ixtrack])
	    || !gsl_finite(pert->pt[itrack][ixtrack])
	    || !gsl_finite(pert->var[itrack][ixtrack])
	    || !gsl_finite(pert->dc[itrack][ixtrack]))
	  continue;

	/* Get and check ascending/descending flag... */
	asc = (pert->lat[itrack > 0 ? itrack : itrack + 1][pert->nxtrack / 2]
	       > pert->lat[itrack >
			   0 ? itrack - 1 : itrack][pert->nxtrack / 2]);
	if (((set[0] == 'a' || set[0] == 'A') && !asc)
	    || ((set[0] == 'd' || set[0] == 'D') && asc))
	  continue;

	/* Check am/pm flag... */
	lt = fmod(pert->time[itrack][ixtrack], 86400.) / 3600.;
	if (((set[0] == 'm' || set[0] == 'M') && lt > 12.)
	    || ((set[0] == 'n' || set[0] == 'N') && lt < 12.))
	  continue;

	/* Get grid indices... */
	ix =
	  (int) ((pert->lon[itrack][ixtrack] - lon0) / (lon1 -
							lon0) * (double) nx);
	iy =
	  (int) ((pert->lat[itrack][ixtrack] - lat0) / (lat1 -
							lat0) * (double) ny);
	if (ix < 0 || ix >= nx || iy < 0 || iy >= ny)
	  continue;

	/* Get month index... */
	imon =
	  (int) (fmod(pert->time[0][0] / 60. / 60. / 24. / 365.25, 1.) *
		 NMON);
	if (imon < 0 || imon >= NMON)
	  continue;

	/* Get gravity wave detection threshold... */
	if (thresh_gw <= 0.0) {
	  ilat = locate_irr(t_gw_lat, NLAT_GW, pert->lat[itrack][ixtrack]);
	  if (asc)
	    t_gw = LIN(t_gw_lat[ilat], t_gw_asc[imon][ilat],
		       t_gw_lat[ilat + 1], t_gw_asc[imon][ilat + 1],
		       pert->lat[itrack][ixtrack]);
	  else
	    t_gw = LIN(t_gw_lat[ilat], t_gw_dsc[imon][ilat],
		       t_gw_lat[ilat + 1], t_gw_dsc[imon][ilat + 1],
		       pert->lat[itrack][ixtrack]);
	} else
	  t_gw = thresh_gw;

	/* Get deep convection detection threshold... */
	if (thresh_dc <= 0.0) {
	  ilat =
	    locate_irr(t_trop_lat, NLAT_TROP, pert->lat[itrack][ixtrack]);
	  t_dc =
	    LIN(t_trop_lat[ilat], t_trop[imon][ilat], t_trop_lat[ilat + 1],
		t_trop[imon][ilat + 1], pert->lat[itrack][ixtrack]) + dt_trop;
	} else
	  t_dc = thresh_dc + dt_trop;

	/* Detection of gravity waves... */
	det_gw = (pert->var[itrack][ixtrack] >= t_gw);

	/* Detection of convective waves... */
	det_cw = 0;
	if (det_gw)
	  for (itrack2 = GSL_MAX(itrack - iradius, 0);
	       itrack2 <= GSL_MIN(itrack + iradius, pert->ntrack - 1);
	       itrack2++)
	    for (ixtrack2 = GSL_MAX(ixtrack - iradius, 0);
		 ixtrack2 <= GSL_MIN(ixtrack + iradius, pert->nxtrack - 1);
		 ixtrack2++) {
	      if (det_cw)
		break;
	      det_cw = (pert->dc[itrack2][ixtrack2] <= t_dc);
	    }

	/* Detection of deep convection... */
	det_dc = (pert->dc[itrack][ixtrack] <= t_dc);

	/* Detection of wave generation... */
	det_wg = 0;
	if (det_dc)
	  for (itrack2 = GSL_MAX(itrack - iradius, 0);
	       itrack2 <= GSL_MIN(itrack + iradius, pert->ntrack - 1);
	       itrack2++)
	    for (ixtrack2 = GSL_MAX(ixtrack - iradius, 0);
		 ixtrack2 <= GSL_MIN(ixtrack + iradius, pert->nxtrack - 1);
		 ixtrack2++) {
	      if (det_wg)
		break;
	      det_wg = (pert->var[itrack2][ixtrack2] >= t_gw);
	    }

	/* Count events... */
	n[ix][iy]++;
	if (det_dc)
	  ndc[ix][iy]++;
	if (det_wg)
	  nwg[ix][iy]++;
	if (det_gw)
	  ngw[ix][iy]++;
	if (det_cw)
	  ncw[ix][iy]++;

	/* Get statistics of perturbations... */
	mean[ix][iy] += pert->pt[itrack][ixtrack];
	var[ix][iy] += gsl_pow_2(pert->pt[itrack][ixtrack]);
	max[ix][iy] = GSL_MAX(max[ix][iy], pert->pt[itrack][ixtrack]);
	min[ix][iy] = GSL_MIN(min[ix][iy], pert->pt[itrack][ixtrack]);

	/* Get statistics of brightness temperatures... */
	bt[ix][iy] += pert->bt[itrack][ixtrack];
	bt_8mu[ix][iy] += pert->dc[itrack][ixtrack];
	if (n[ix][iy] > 1) {
	  bt_8mu_min[ix][iy]
	    = GSL_MIN(bt_8mu_min[ix][iy], pert->dc[itrack][ixtrack]);
	  bt_8mu_max[ix][iy]
	    = GSL_MAX(bt_8mu_max[ix][iy], pert->dc[itrack][ixtrack]);
	} else {
	  bt_8mu_min[ix][iy] = pert->dc[itrack][ixtrack];
	  bt_8mu_max[ix][iy] = pert->dc[itrack][ixtrack];
	}

	/* Get mean time... */
	mtime[ix][iy] += pert->time[itrack][ixtrack];
      }
  }

  /* Analyze results... */
  for (ix = 0; ix < nx; ix++)
    for (iy = 0; iy < ny; iy++) {

      /* Get geolocation... */
      mtime[ix][iy] /= (double) n[ix][iy];
      glon[ix]
      = lon0 + (ix + 0.5) / (double) nx *(
  lon1 - lon0);
      glat[iy]
      = lat0 + (iy + 0.5) / (double) ny *(
  lat1 - lat0);

      /* Normalize brightness temperatures... */
      bt[ix][iy] /= (double) n[ix][iy];
      bt_8mu[ix][iy] /= (double) n[ix][iy];

      /* Get fractions... */
      fdc[ix][iy] = (double) ndc[ix][iy] / (double) n[ix][iy] * 100.;
      fwg[ix][iy] = (double) nwg[ix][iy] / (double) ndc[ix][iy] * 100.;
      fgw[ix][iy] = (double) ngw[ix][iy] / (double) n[ix][iy] * 100.;
      fcw[ix][iy] = (double) ncw[ix][iy] / (double) ngw[ix][iy] * 100.;

      /* Check number of observations... */
      if (n[ix][iy] < nmin) {
	fdc[ix][iy] = GSL_NAN;
	fwg[ix][iy] = GSL_NAN;
	fgw[ix][iy] = GSL_NAN;
	fcw[ix][iy] = GSL_NAN;
	bt_8mu[ix][iy] = GSL_NAN;
	bt_8mu_min[ix][iy] = GSL_NAN;
	bt_8mu_max[ix][iy] = GSL_NAN;
      }

      /* Check detections of deep convection at high latitudes... */
      if (fabs(glat[iy]) > dc_hlat && bt_8mu[ix][iy] <= dc_tlim) {
	fdc[ix][iy] = GSL_NAN;
	fwg[ix][iy] = GSL_NAN;
	fcw[ix][iy] = GSL_NAN;
      }

      /* Estimate noise... */
      if (dt230 > 0) {
	nesr = PLANCK(230.0 + dt230, nu) - PLANCK(230.0, nu);
	dt[ix][iy] = BRIGHT(PLANCK(bt[ix][iy], nu) + nesr, nu) - bt[ix][iy];
      }

      /* Get mean perturbation and variance... */
      mean[ix][iy] /= (double) n[ix][iy];
      var[ix][iy] =
	var[ix][iy] / (double) n[ix][iy] - gsl_pow_2(mean[ix][iy]);
    }

  /* Write ASCII file... */
  if (output == 1) {

    /* Create file... */
    printf("Write variance statistics: %s\n", argv[2]);
    if (!(out = fopen(argv[2], "w")))
      ERRMSG("Cannot create file!");

    /* Write header... */
    fprintf(out,
	    "# $1 = time [s]\n"
	    "# $2 = longitude [deg]\n"
	    "# $3 = latitude [deg]\n"
	    "# $4 = number of footprints\n"
	    "# $5 = fraction of convection events [%%]\n"
	    "# $6 = fraction of wave generating events [%%]\n"
	    "# $7 = fraction of gravity wave events [%%]\n"
	    "# $8 = fraction of convective wave events [%%]\n"
	    "# $9 = mean perturbation [K]\n"
	    "# $10 = minimum perturbation [K]\n");
    fprintf(out,
	    "# $11 = maximum perturbation [K]\n"
	    "# $12 = variance [K^2]\n"
	    "# $13 = mean surface temperature [K]\n"
	    "# $14 = minimum surface temperature [K]\n"
	    "# $15 = maximum surface temperature [K]\n"
	    "# $16 = mean background temperature [K]\n"
	    "# $17 = noise estimate [K]\n");

    /* Write results... */
    for (iy = 0; iy < ny; iy++) {
      if (iy == 0 || nx > 1)
	fprintf(out, "\n");
      for (ix = 0; ix < nx; ix++)
	fprintf(out, "%.2f %g %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
		mtime[ix][iy], glon[ix], glat[iy], n[ix][iy],
		fdc[ix][iy], fwg[ix][iy], fgw[ix][iy], fcw[ix][iy],
		mean[ix][iy], min[ix][iy], max[ix][iy], var[ix][iy],
		bt_8mu[ix][iy], bt_8mu_min[ix][iy], bt_8mu_max[ix][iy],
		bt[ix][iy], dt[ix][iy]);
    }

    /* Close file... */
    fclose(out);
  }

  /* Write netCDF file... */
  else if (output == 2) {

    /* Create netCDF file... */
    printf("Write variance statistics: %s\n", argv[2]);
    NC(nc_create(argv[2], NC_CLOBBER, &ncid));

    /* Set dimensions... */
    NC(nc_def_dim(ncid, "lat", (size_t) ny, &dimid[0]));
    NC(nc_def_dim(ncid, "lon", (size_t) nx, &dimid[1]));

    /* Add variables... */
    NC(nc_def_var(ncid, "lat", NC_DOUBLE, 1, &dimid[0], &latid));
    add_att(ncid, latid, "deg", "latitude");
    NC(nc_def_var(ncid, "lon", NC_DOUBLE, 1, &dimid[1], &lonid));
    add_att(ncid, lonid, "deg", "longitude");
    NC(nc_def_var(ncid, "var", NC_FLOAT, 2, dimid, &varid));
    add_att(ncid, varid, "K^2", "brightness temperature variance");
    NC(nc_def_var(ncid, "min", NC_FLOAT, 2, dimid, &minid));
    add_att(ncid, minid, "K", "brightness temperature minimum");
    NC(nc_def_var(ncid, "max", NC_FLOAT, 2, dimid, &maxid));
    add_att(ncid, maxid, "K", "brightness temperature maximum");
    NC(nc_def_var(ncid, "np", NC_INT, 2, dimid, &npid));
    add_att(ncid, npid, "1", "number of footprints");

    /* Leave define mode... */
    NC(nc_enddef(ncid));

    /* Write data... */
    NC(nc_put_var_double(ncid, latid, glat));
    NC(nc_put_var_double(ncid, lonid, glon));
    for (ix = 0; ix < nx; ix++)
      for (iy = 0; iy < ny; iy++)
	help[iy * nx + ix] = var[ix][iy] - POW2(dt[ix][iy]);
    NC(nc_put_var_double(ncid, varid, help));
    for (ix = 0; ix < nx; ix++)
      for (iy = 0; iy < ny; iy++)
	help[iy * nx + ix] = min[ix][iy];
    NC(nc_put_var_double(ncid, minid, help));
    for (ix = 0; ix < nx; ix++)
      for (iy = 0; iy < ny; iy++)
	help[iy * nx + ix] = max[ix][iy];
    NC(nc_put_var_double(ncid, maxid, help));
    for (ix = 0; ix < nx; ix++)
      for (iy = 0; iy < ny; iy++)
	help2[iy * nx + ix] = n[ix][iy];
    NC(nc_put_var_int(ncid, npid, help2));

    /* Close file... */
    NC(nc_close(ncid));
  }

  else
    ERRMSG("Unknown output format!");

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}
