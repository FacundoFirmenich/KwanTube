# 1JFF Medium-Horizon Validation Summary

## Configuration
- h_source: C:\Users\User\3D Objects\biofisicaquantiqaCLINE\H_1JFF.npz
- Working level: NC=5, Nk=1
- Init site: 5 (B:103)
- Horizon: 200.0 fs
- Output sampling: 5.0 fs
- Threshold: 0.001

## Runtimes
- NC_minus_1: NC=4, Nk=1, wall=15.6 s, local 30 ps diagnostic=39.1 min
- baseline: NC=5, Nk=1, wall=83.8 s, local 30 ps diagnostic=209.5 min
- NC_stress: NC=6, Nk=1, wall=325.4 s, local 30 ps diagnostic=813.4 min
- Nk_stress: NC=5, Nk=2, wall=1404.7 s, local 30 ps diagnostic=3511.7 min

## NC-stress summary
- Full-window max differences: dPop=2.92e-03, dCoh=9.60e-03, dFrob=1.72e-02
- t_dPop_max: 175.0 fs (site 7)
- t_dCoh_max: 170.0 fs
- t_dFrob_max: 180.0 fs
- Final differences: dPop=2.17e-03, dCoh=7.46e-03, dFrob=1.56e-02
- Full-window threshold verdict (dFrob < 1.0e-03): **FAIL**
- Medium-horizon convergence ratio r = 0.478

### Windowed NC-stress maxima
| Window | dPop max | dCoh max | dFrob max |
|---|---:|---:|---:|
| 0-40fs | 1.79e-06 | 5.22e-05 | 7.63e-05 |
| 0-80fs | 1.79e-05 | 1.41e-04 | 2.16e-04 |
| 0-120fs | 1.25e-04 | 3.32e-03 | 5.14e-03 |
| 0-160fs | 2.43e-03 | 9.16e-03 | 1.55e-02 |
| 0-200fs | 2.92e-03 | 9.60e-03 | 1.72e-02 |

## Nk-stress summary
- Full-window max differences: dPop=4.64e-05, dCoh=2.52e-04, dFrob=3.73e-04
- t_dPop_max: 60.0 fs (site 5)
- Final differences: dPop=1.01e-05, dFrob=1.07e-04
- Full-window threshold verdict (dFrob < 1.0e-03): **PASS**

## Manuscript-safe interpretation
- This medium-horizon validation supports the chosen working level only over the tested time interval; it does not by itself prove asymptotic convergence at all times.
