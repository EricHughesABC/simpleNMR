import pandas as pd
from pathlib import Path
import tempfile
import yaml
import os

examplepromlems_yaml_str = """2-ethyl-1-indanone:
  C13_1D:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
    Area:
      1: 2525.3
      2: 2771.44
      3: 2384.25
      4: 9412.76
      5: 8305.46
      6: 7179.36
      7: 7995.31
      8: 8156.34
      9: 9129.5
      10: 9409.77
      11: 10051.5
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
    Intensity:
      1: 381.3
      2: 737.2
      3: 542.4
      4: 2086.5
      5: 2093.5
      6: 2023.4
      7: 1945.5
      8: 2002.2
      9: 2194.1
      10: 2280.4
      11: 2246.7
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
    Width:
      1: 1.45
      2: 0.89
      3: 1.02
      4: 1.01
      5: 0.93
      6: 0.85
      7: 0.9
      8: 0.94
      9: 0.89
      10: 0.97
      11: 1.01
    ppm:
      1: 208.934
      2: 153.81
      3: 136.955
      4: 134.588
      5: 127.274
      6: 126.523
      7: 123.817
      8: 48.743
      9: 32.318
      10: 24.453
      11: 11.582
  COSY:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
    Intensity:
      1: 0.813
      2: 1.739
      3: 0.291
      4: 0.337
      5: 1.603
      6: 0.093
      7: 0.305
      8: 0.018
      9: 0.342
      10: 0.049
      11: 0.291
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
    Volume:
      1: 10.19
      2: 28.48
      3: 3.13
      4: 8.32
      5: 28.1
      6: 1.74
      7: 7.52
      8: 0.41
      9: 3.44
      10: 1.47
      11: 4.61
    Width f1:
      1: 15.89
      2: 26.61
      3: 16.96
      4: 26.62
      5: 17.82
      6: 22.7
      7: 34.44
      8: 19.26
      9: 14.61
      10: 24.05
      11: 13.55
    Width f2:
      1: 5.78
      2: 4.51
      3: 4.64
      4: 6.81
      5: 7.21
      6: 6.03
      7: 5.25
      8: 8.4
      9: 5.05
      10: 9.02
      11: 8.56
    f1 (ppm):
      1: 7.35
      2: 7.45
      3: 7.35
      4: 7.45
      5: 2.83
      6: 2.61
      7: 2.81
      8: 2.61
      9: 2.61
      10: 1.53
      11: 1.0
    f2 (ppm):
      1: 7.739
      2: 7.57
      3: 7.568
      4: 7.448
      5: 3.309
      6: 3.307
      7: 2.811
      8: 2.811
      9: 2.609
      10: 2.606
      11: 1.968
  H1_1D:
    Class:
      1: ddq
      2: td
      3: dp
      4: ddt
      5: ddq
      6: dd
      7: dddd
      8: dqd
      9: m
      10: t
    H's:
      1: 1
      2: 1
      3: 1
      4: 1
      5: 1
      6: 1
      7: 1
      8: 1
      9: 1
      10: 3
    Integral:
      1: 0.94
      2: 1.02
      3: 1.01
      4: 1.03
      5: 1.04
      6: 0.95
      7: 1.07
      8: 1.05
      9: 1.05
      10: 3.0
    J's:
      1: 0.57, 0.57, 0.56, 1.25, 7.69
      2: 1.24, 7.42, 7.38
      3: 0.95, 0.95, 0.94, 0.94, 7.68
      4: 0.86, 0.86, 7.22, 8.52
      5: 0.78, 0.78, 0.77, 7.89, 17.07
      6: 3.92, 17.09
      7: 3.92, 4.59, 7.86, 9.14
      8: 4.56, 7.50, 7.50, 7.50, 13.71
      9: 0
      10: 7.42, 7.42
    Method:
      1: Peaks
      2: Peaks
      3: Peaks
      4: Peaks
      5: Peaks
      6: Peaks
      7: Peaks
      8: Peaks
      9: Peaks
      10: Peaks
    Name:
      1: A (ddq)
      2: B (td)
      3: C (dp)
      4: D (ddt)
      5: E (ddq)
      6: F (dd)
      7: G (dddd)
      8: H (dqd)
      9: I (m)
      10: J (t)
    Range:
      1: 7.77 .. 7.72
      2: 7.60 .. 7.54
      3: 7.47 .. 7.43
      4: 7.38 .. 7.32
      5: 3.34 .. 3.28
      6: 2.84 .. 2.78
      7: 2.64 .. 2.57
      8: 2.01 .. 1.92
      9: 1.58 .. 1.49
      10: 1.03 .. 0.97
    Shift:
      1: 7.74
      2: 7.57
      3: 7.45
      4: 7.35
      5: 3.31
      6: 2.81
      7: 2.61
      8: 1.97
      9: 1.53
      10: 1.0
  H1_pureshift:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Area:
      1: 9.57
      2: 14.49
      3: 13.64
      4: 16.45
      5: 46.81
      6: 54.31
      7: 20.46
      8: 38.04
      9: 37.6
      10: 108.29
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Intensity:
      1: 0.5
      2: 0.8
      3: 0.8
      4: 0.9
      5: 2.8
      6: 2.9
      7: 1.2
      8: 2.1
      9: 1.8
      10: 5.8
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
    Width:
      1: 1.76
      2: 1.61
      3: 1.61
      4: 1.65
      5: 1.52
      6: 1.57
      7: 1.42
      8: 1.51
      9: 1.68
      10: 1.57
    ppm:
      1: 7.74
      2: 7.57
      3: 7.45
      4: 7.35
      5: 3.31
      6: 2.81
      7: 2.61
      8: 1.97
      9: 1.53
      10: 1.0
  HMBC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
      23: .nan
      24: .nan
      25: .nan
      26: .nan
      27: .nan
      28: .nan
      29: .nan
      30: .nan
      31: .nan
      32: .nan
      33: .nan
      34: .nan
      35: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
      14: None
      15: None
      16: None
      17: None
      18: None
      19: None
      20: None
      21: None
      22: None
      23: None
      24: None
      25: None
      26: None
      27: None
      28: None
      29: None
      30: None
      31: None
      32: None
      33: None
      34: None
      35: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
      23: .nan
      24: .nan
      25: .nan
      26: .nan
      27: .nan
      28: .nan
      29: .nan
      30: .nan
      31: .nan
      32: .nan
      33: .nan
      34: .nan
      35: .nan
    Intensity:
      1: 0.032
      2: 0.016
      3: 0.017
      4: 0.031
      5: 0.03
      6: 0.014
      7: 0.015
      8: 0.018
      9: 0.053
      10: 0.03
      11: 0.043
      12: 0.018
      13: 0.024
      14: 0.062
      15: 0.019
      16: 0.073
      17: 0.014
      18: 0.015
      19: 0.02
      20: 0.003
      21: 0.008
      22: 0.003
      23: 0.043
      24: 0.018
      25: 0.031
      26: 0.04
      27: 0.025
      28: 0.038
      29: 0.016
      30: 0.016
      31: 0.007
      32: 0.027
      33: 0.017
      34: 0.288
      35: 0.136
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
      14: Compound
      15: Compound
      16: Compound
      17: Compound
      18: Compound
      19: Compound
      20: Compound
      21: Compound
      22: Compound
      23: Compound
      24: Compound
      25: Compound
      26: Compound
      27: Compound
      28: Compound
      29: Compound
      30: Compound
      31: Compound
      32: Compound
      33: Compound
      34: Compound
      35: Compound
    Volume:
      1: 0.6
      2: 0.28
      3: 0.27
      4: 0.64
      5: 0.64
      6: 0.41
      7: 0.26
      8: 0.22
      9: 1.44
      10: 0.57
      11: 0.81
      12: 0.26
      13: 0.7
      14: 1.5
      15: 0.26
      16: 1.48
      17: 0.28
      18: 0.37
      19: 0.41
      20: 0.07
      21: 0.24
      22: 0.07
      23: 2.45
      24: 0.91
      25: 0.52
      26: 2.37
      27: 2.09
      28: 0.85
      29: 0.24
      30: 0.17
      31: 0.08
      32: 0.45
      33: 0.19
      34: 3.9
      35: 2.56
    Width f1:
      1: 77.73
      2: 87.93
      3: 54.24
      4: 79.55
      5: 84.84
      6: 82.91
      7: 87.04
      8: 69.68
      9: 80.34
      10: 83.6
      11: 63.39
      12: 54.81
      13: 86.42
      14: 89.84
      15: 58.75
      16: 56.13
      17: 90.27
      18: 66.97
      19: 57.02
      20: 56.84
      21: 90.05
      22: 81.86
      23: 49.79
      24: 55.19
      25: 46.73
      26: 83.32
      27: 89.58
      28: 80.52
      29: 62.6
      30: 53.74
      31: 51.26
      32: 53.43
      33: 53.47
      34: 66.13
      35: 50.52
    Width f2:
      1: 5.98
      2: 4.88
      3: 7.29
      4: 6.46
      5: 6.36
      6: 8.94
      7: 4.97
      8: 4.2
      9: 8.54
      10: 5.7
      11: 7.52
      12: 6.8
      13: 8.55
      14: 6.77
      15: 5.83
      16: 9.02
      17: 5.72
      18: 9.36
      19: 8.87
      20: 8.42
      21: 8.65
      22: 8.22
      23: 28.56
      24: 22.63
      25: 9.02
      26: 18.06
      27: 23.21
      28: 6.95
      29: 6.01
      30: 5.02
      31: 5.17
      32: 7.69
      33: 5.29
      34: 5.13
      35: 9.34
    f1 (ppm):
      1: 134.55
      2: 153.74
      3: 208.96
      4: 123.79
      5: 153.92
      6: 32.34
      7: 137.01
      8: 127.26
      9: 134.66
      10: 137.09
      11: 48.78
      12: 126.56
      13: 136.93
      14: 153.87
      15: 24.45
      16: 208.93
      17: 153.86
      18: 48.79
      19: 24.5
      20: 126.59
      21: 137.02
      22: 134.55
      23: 24.42
      24: 11.65
      25: 209.09
      26: 32.32
      27: 153.76
      28: 32.13
      29: 48.79
      30: 209.05
      31: 11.63
      32: 208.94
      33: 11.56
      34: 48.77
      35: 24.69
    f2 (ppm):
      1: 7.739
      2: 7.739
      3: 7.738
      4: 7.567
      5: 7.567
      6: 7.45
      7: 7.45
      8: 7.449
      9: 7.448
      10: 7.347
      11: 3.31
      12: 3.31
      13: 3.31
      14: 3.31
      15: 3.309
      16: 3.308
      17: 2.814
      18: 2.813
      19: 2.813
      20: 2.813
      21: 2.812
      22: 2.812
      23: 2.607
      24: 2.607
      25: 2.607
      26: 2.607
      27: 2.606
      28: 1.967
      29: 1.966
      30: 1.966
      31: 1.964
      32: 1.532
      33: 1.532
      34: 1.003
      35: 1.002
  HSQC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Intensity:
      1: 0.491
      2: 0.467
      3: 0.53
      4: 0.742
      5: -0.316
      6: -0.008
      7: 0.572
      8: -0.353
      9: -0.335
      10: 2.177
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
    Volume:
      1: 5.55
      2: 5.84
      3: 4.19
      4: 6.74
      5: -1.31
      6: -0.05
      7: 6.43
      8: -8.87
      9: -8.87
      10: 27.48
    Width f1:
      1: 53.61
      2: 62.29
      3: 43.72
      4: 40.73
      5: 46.56
      6: 49.88
      7: 43.26
      8: 64.48
      9: 63.87
      10: 61.1
    Width f2:
      1: 17.65
      2: 16.8
      3: 15.12
      4: 18.7
      5: 7.45
      6: 11.81
      7: 21.78
      8: 32.64
      9: 34.74
      10: 17.3
    f1 (ppm):
      1: 123.82
      2: 134.59
      3: 126.53
      4: 127.27
      5: 32.32
      6: 32.32
      7: 48.75
      8: 24.45
      9: 24.45
      10: 11.58
    f2 (ppm):
      1: 7.737
      2: 7.567
      3: 7.45
      4: 7.349
      5: 3.311
      6: 2.812
      7: 2.606
      8: 1.966
      9: 1.533
      10: 1.002
  NOESY:
    Annotation:
      1: .nan
      2: .nan
    Flags:
      1: None
      2: None
    Impurity/Compound:
      1: .nan
      2: .nan
    Intensity:
      1: 0.0
      2: 0.2
    Type:
      1: Compound
      2: Compound
    Volume:
      1: -1.11
      2: 1.22
    Width f1:
      1: 14.3
      2: 10.7
    Width f2:
      1: 6.43
      2: 2.73
    f1 (ppm):
      1: 2.82
      2: 1.53
    f2 (ppm):
      1: 3.31
      2: 1.97
  molecule:
    smiles:
      1: CCC2Cc1ccccc1C2=O
Cellobiose:
  C13_1D:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
    Area:
      1: 2450.5
      2: 2400.06
      3: 2434.81
      4: 2320.32
      5: 2301.38
      6: 2515.24
      7: 2527.89
      8: 1944.65
      9: 3148.26
      10: 2424.66
      11: 2406.9
      12: 2454.06
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
    Intensity:
      1: 101.6
      2: 82.0
      3: 95.3
      4: 86.8
      5: 82.8
      6: 76.5
      7: 91.2
      8: 67.7
      9: 87.0
      10: 60.3
      11: 60.3
      12: 63.3
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
    Width:
      1: 1.3
      2: 1.61
      3: 1.39
      4: 1.5
      5: 1.56
      6: 1.88
      7: 1.52
      8: 1.73
      9: 2.12
      10: 2.45
      11: 2.28
      12: 2.18
    ppm:
      1: 100.91
      2: 88.92
      3: 76.06
      4: 72.92
      5: 71.97
      6: 71.6
      7: 70.71
      8: 69.29
      9: 69.28
      10: 67.73
      11: 61.53
      12: 61.29
  COSY:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
      14: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
    Intensity:
      1: 4.52
      2: 0.24
      3: 0.01
      4: 0.13
      5: 0.05
      6: 0.07
      7: 1.96
      8: 0.1
      9: 0.17
      10: 1.22
      11: 2.87
      12: 2.42
      13: 0.2
      14: 0.24
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
      14: Compound
    Volume:
      1: 156.27
      2: 13.89
      3: 0.35
      4: 9.84
      5: 1.79
      6: 2.49
      7: 63.25
      8: 2.78
      9: 19.24
      10: 83.4
      11: 159.67
      12: 126.19
      13: 15.57
      14: 9.52
    Width f1:
      1: 13.47
      2: 19.07
      3: 11.62
      4: 31.55
      5: 15.07
      6: 16.18
      7: 11.68
      8: 12.02
      9: 27.23
      10: 13.89
      11: 12.94
      12: 12.99
      13: 22.48
      14: 10.6
    Width f2:
      1: 14.79
      2: 17.46
      3: 11.62
      4: 13.56
      5: 14.84
      6: 13.05
      7: 15.93
      8: 13.51
      9: 24.34
      10: 28.45
      11: 24.75
      12: 23.16
      13: 19.94
      14: 21.45
    f1 (ppm):
      1: 4.99
      2: 3.77
      3: 4.99
      4: 5.062
      5: 4.926
      6: 3.646
      7: 4.493
      8: 4.095
      9: 3.973
      10: 4.027
      11: 3.647
      12: 3.974
      13: 3.648
      14: 3.77
    f2 (ppm):
      1: 6.23
      2: 5.419
      3: 5.418
      4: 5.128
      5: 5.126
      6: 5.058
      7: 4.925
      8: 4.455
      9: 4.455
      10: 4.371
      11: 4.371
      12: 4.094
      13: 4.026
      14: 3.974
  H1_1D:
    Class:
      1: d
      2: dd
      3: t
      4: m
      5: dd
      6: dd
      7: d
      8: dd
      9: dd
      10: dd
      11: dd
      12: ddd
      13: dd
      14: ddd
    H's:
      1: 1
      2: 1
      3: 1
      4: 1
      5: 1
      6: 1
      7: 1
      8: 1
      9: 1
      10: 1
      11: 1
      12: 1
      13: 1
      14: 1
    Integral:
      1: 1.0
      2: 1.0
      3: 1.01
      4: 1.03
      5: 0.99
      6: 1.04
      7: 0.98
      8: 0.96
      9: 1.02
      10: 0.98
      11: 1.02
      12: 0.98
      13: 1.05
      14: 0.99
    J's:
      1: 3.73
      2: 9.26, 10.35
      3: 9.38, 9.38
      4: .nan
      5: 3.76, 10.34
      6: 7.94, 9.37
      7: 7.95
      8: 2.13, 12.25
      9: 4.34, 12.50
      10: 4.26, 12.26
      11: 2.33, 12.47
      12: 2.18, 4.29, 10.10
      13: 9.27, 10.15
      14: 2.33, 4.34, 9.98
    Method:
      1: Sum
      2: Sum
      3: Sum
      4: Sum
      5: Sum
      6: Sum
      7: Sum
      8: Sum
      9: Sum
      10: Sum
      11: Sum
      12: Sum
      13: Sum
      14: Sum
    Name:
      1: A (d)
      2: B (dd)
      3: C (t)
      4: D (m)
      5: E (dd)
      6: F (dd)
      7: G (d)
      8: H (dd)
      9: I (dd)
      10: J (dd)
      11: K (dd)
      12: L (ddd)
      13: M (dd)
      14: N (ddd)
    Range:
      1: 6.26 .. 6.20
      2: 5.45 .. 5.38
      3: 5.16 .. 5.10
      4: 5.09 .. 5.02
      5: 5.02 .. 4.96
      6: 4.95 .. 4.87
      7: 4.53 .. 4.48
      8: 4.48 .. 4.43
      9: 4.40 .. 4.34
      10: 4.13 .. 4.06
      11: 4.06 .. 4.00
      12: 4.00 .. 3.94
      13: 3.80 .. 3.74
      14: 3.68 .. 3.61
    Shift:
      1: 6.23
      2: 5.42
      3: 5.13
      4: 5.06
      5: 4.99
      6: 4.92
      7: 4.49
      8: 4.46
      9: 4.37
      10: 4.09
      11: 4.03
      12: 3.97
      13: 3.77
      14: 3.65
  H1_pureshift:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
    Area:
      1: 9.57
      2: 14.49
      3: 19.41
      4: 20.0
      5: 29.25
      6: 34.17
      7: 20.0
      8: 20.0
      9: 20.0
      10: 20.0
      11: 20.0
      12: 20.0
      13: 20.0
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
    Intensity:
      1: 0.5
      2: 0.8
      3: 1.1
      4: 1.0
      5: 1.7
      6: 2.0
      7: 1.0
      8: 1.0
      9: 1.0
      10: 1.0
      11: 1.0
      12: 1.0
      13: 1.0
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
    Width:
      1: 1.76
      2: 1.61
      3: 1.46
      4: 1.0
      5: 1.16
      6: 1.01
      7: 1.0
      8: 1.0
      9: 1.0
      10: 1.0
      11: 1.0
      12: 1.0
      13: 1.0
    ppm:
      1: 6.23
      2: 5.42
      3: 5.13
      4: 3.65
      5: 4.99
      6: 4.92
      7: 4.49
      8: 4.46
      9: 4.37
      10: 4.09
      11: 4.03
      12: 3.97
      13: 3.77
  HMBC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
      14: None
      15: None
      16: None
      17: None
      18: None
      19: None
      20: None
      21: None
      22: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
    Intensity:
      1: 0.02
      2: 0.01
      3: 0.02
      4: 0.02
      5: 0.0
      6: 0.02
      7: 0.03
      8: 0.02
      9: 0.01
      10: 0.02
      11: 0.02
      12: 0.02
      13: 0.02
      14: 0.05
      15: 0.03
      16: 0.02
      17: 0.02
      18: 0.02
      19: 0.03
      20: 0.02
      21: 0.01
      22: 0.04
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
      14: Compound
      15: Compound
      16: Compound
      17: Compound
      18: Compound
      19: Compound
      20: Compound
      21: Compound
      22: Compound
    Volume:
      1: 0.45
      2: 0.27
      3: 0.26
      4: 0.3
      5: 0.01
      6: 0.31
      7: 0.38
      8: 0.26
      9: 0.17
      10: 0.3
      11: 0.54
      12: 0.66
      13: 0.78
      14: 1.22
      15: 0.61
      16: 0.42
      17: 0.62
      18: 0.5
      19: 0.43
      20: 0.21
      21: 0.19
      22: 0.42
    Width f1:
      1: 48.03
      2: 54.23
      3: 52.07
      4: 51.13
      5: 38.74
      6: 52.86
      7: 55.57
      8: 49.3
      9: 50.56
      10: 54.53
      11: 45.15
      12: 55.54
      13: 58.26
      14: 52.98
      15: 50.94
      16: 54.22
      17: 53.53
      18: 55.71
      19: 50.93
      20: 56.09
      21: 55.84
      22: 45.07
    Width f2:
      1: 21.87
      2: 15.48
      3: 12.75
      4: 11.48
      5: 21.82
      6: 10.47
      7: 11.69
      8: 11.03
      9: 13.26
      10: 12.21
      11: 23.34
      12: 30.13
      13: 33.48
      14: 21.88
      15: 20.23
      16: 21.75
      17: 24.06
      18: 23.9
      19: 12.38
      20: 10.96
      21: 12.51
      22: 11.61
    f1 (ppm):
      1: 69.291
      2: 70.707
      3: 76.072
      4: 69.308
      5: 78.397
      6: 71.607
      7: 67.736
      8: 72.941
      9: 61.54
      10: 71.99
      11: 69.282
      12: 72.941
      13: 100.911
      14: 76.08
      15: 76.08
      16: 70.73
      17: 67.757
      18: 71.991
      19: 69.28
      20: 70.73
      21: 61.275
      22: 100.911
    f2 (ppm):
      1: 6.23
      2: 6.23
      3: 5.419
      4: 5.419
      5: 5.321
      6: 5.128
      7: 5.128
      8: 5.059
      9: 5.059
      10: 5.059
      11: 4.99
      12: 4.925
      13: 4.924
      14: 4.494
      15: 4.457
      16: 4.456
      17: 4.028
      18: 4.027
      19: 3.77
      20: 3.77
      21: 3.77
      22: 3.769
  HSQC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
      14: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
    Intensity:
      1: 0.06
      2: 0.05
      3: 0.04
      4: 0.05
      5: 0.04
      6: 0.05
      7: 0.06
      8: -0.02
      9: -0.02
      10: -0.02
      11: -0.02
      12: 0.03
      13: 0.05
      14: 0.04
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
      14: Compound
    Volume:
      1: 1.52
      2: 1.4
      3: 0.9
      4: 1.35
      5: 1.06
      6: 2.53
      7: 1.7
      8: -0.44
      9: -0.62
      10: -0.75
      11: -0.59
      12: 1.21
      13: 1.19
      14: 1.29
    Width f1:
      1: 45.72
      2: 44.09
      3: 44.59
      4: 45.63
      5: 43.06
      6: 45.35
      7: 44.81
      8: 44.27
      9: 45.32
      10: 45.48
      11: 44.41
      12: 47.29
      13: 44.81
      14: 46.27
    Width f2:
      1: 10.01
      2: 11.52
      3: 9.05
      4: 11.96
      5: 12.57
      6: 20.82
      7: 12.24
      8: 10.18
      9: 10.5
      10: 12.46
      11: 10.25
      12: 13.61
      13: 10.19
      14: 12.25
    f1 (ppm):
      1: 88.92
      2: 69.298
      3: 72.937
      4: 67.728
      5: 69.284
      6: 71.608
      7: 100.913
      8: 61.311
      9: 61.546
      10: 61.311
      11: 61.546
      12: 70.702
      13: 76.062
      14: 71.983
    f2 (ppm):
      1: 6.23
      2: 5.419
      3: 5.128
      4: 5.059
      5: 4.989
      6: 4.925
      7: 4.494
      8: 4.456
      9: 4.37
      10: 4.093
      11: 4.025
      12: 3.973
      13: 3.77
      14: 3.645
  NOESY:
    Annotation: {}
    Flags: {}
    Impurity/Compound: {}
    Intensity: {}
    Type: {}
    Volume: {}
    Width f1: {}
    Width f2: {}
    f1 (ppm): {}
    f2 (ppm): {}
  Sheet1: {}
  molecule:
    molecule:
      1: " \tC12H22O11"
    smiles:
      1: O[C@H]2[C@H](O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](OC(O)[C@@H]2O)CO
Fasiglifam:
  C13_1D:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
      23: .nan
      24: .nan
      25: .nan
      26: .nan
    Area:
      1: 1751.42
      2: 1184.97
      3: 1435.97
      4: 1689.86
      5: 1214.91
      6: 1315.41
      7: 2568.81
      8: 1226.81
      9: 1978.81
      10: 2231.45
      11: 2069.52
      12: 2088.73
      13: 2188.93
      14: 1185.6
      15: 4255.7
      16: 2682.18
      17: 2524.58
      18: 2925.4
      19: 2549.64
      20: 2825.45
      21: 2784.08
      22: 4441.58
      23: 3592.14
      24: 2555.95
      25: 2854.46
      26: 6831.01
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
      14: None
      15: None
      16: None
      17: None
      18: None
      19: None
      20: None
      21: None
      22: None
      23: None
      24: None
      25: None
      26: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
      23: .nan
      24: .nan
      25: .nan
      26: .nan
    Intensity:
      1: 668.4
      2: 580.2
      3: 503.7
      4: 459.8
      5: 412.3
      6: 478.6
      7: 821.4
      8: 447.8
      9: 410.9
      10: 504.4
      11: 479.5
      12: 325.5
      13: 617.2
      14: 487.1
      15: 1041.3
      16: 720.6
      17: 703.5
      18: 642.9
      19: 530.2
      20: 616.0
      21: 901.0
      22: 758.6
      23: 613.5
      24: 660.5
      25: 701.4
      26: 1511.8
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
      14: Compound
      15: Compound
      16: Compound
      17: Compound
      18: Compound
      19: Compound
      20: Compound
      21: Compound
      22: Compound
      23: Compound
      24: Compound
      25: Compound
      26: Compound
    Width:
      1: 1.65
      2: 1.38
      3: 1.74
      4: 2.31
      5: 1.81
      6: 1.73
      7: 1.9
      8: 1.54
      9: 2.91
      10: 2.61
      11: 2.61
      12: 3.51
      13: 2.2
      14: 1.39
      15: 2.54
      16: 2.12
      17: 2.13
      18: 2.69
      19: 2.92
      20: 2.79
      21: 1.99
      22: 3.39
      23: 3.39
      24: 2.56
      25: 2.71
      26: 2.65
    ppm:
      1: 173.58
      2: 161.15
      3: 159.58
      4: 157.43
      5: 140.72
      6: 137.83
      7: 137.09
      8: 134.41
      9: 129.23
      10: 129.06
      11: 129.02
      12: 126.31
      13: 125.02
      14: 122.42
      15: 113.78
      16: 107.45
      17: 97.4
      18: 77.58
      19: 69.8
      20: 65.88
      21: 51.05
      22: 40.71
      23: 39.61
      24: 37.55
      25: 22.55
      26: 21.21
  COSY:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Intensity:
      1: 62.52
      2: 45.66
      3: 39.67
      4: 42.56
      5: 4.79
      6: 4.93
      7: 9.55
      8: 11.34
      9: 4.98
      10: 10.95
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
    Volume:
      1: 5059.0
      2: 5259.19
      3: 4024.85
      4: 5479.18
      5: 487.57
      6: 520.17
      7: 1014.46
      8: 863.3
      9: 374.52
      10: 1112.32
    Width f1:
      1: 58.64
      2: 63.61
      3: 54.96
      4: 64.83
      5: 61.65
      6: 63.81
      7: 61.5
      8: 59.4
      9: 58.73
      10: 67.0
    Width f2:
      1: 11.85
      2: 15.55
      3: 15.85
      4: 17.05
      5: 14.19
      6: 14.21
      7: 14.84
      8: 11.01
      9: 11.0
      10: 13.02
    f1 (ppm):
      1: 7.447
      2: 7.447
      3: 6.457
      4: 7.103
      5: 3.678
      6: 4.185
      7: 3.671
      8: 2.148
      9: 2.148
      10: 2.488
    f2 (ppm):
      1: 7.385
      2: 7.062
      3: 6.48
      4: 6.478
      5: 4.685
      6: 4.682
      7: 4.192
      8: 4.097
      9: 3.284
      10: 2.695
  H1_1D:
    Class:
      1: t
      2: dt
      3: d
      4: dd
      5: dt
      6: s
      7: dd
      8: d
      9: s
      10: t
      11: dd
      12: t
      13: ddd
      14: m
      15: s
      16: dd
      17: m
      18: m
      19: s
    H's:
      1: 1
      2: 1
      3: 1
      4: 1
      5: 1
      6: 1
      7: 1
      8: 1
      9: 2
      10: 1
      11: 1
      12: 2
      13: 1
      14: 2
      15: 3
      16: 1
      17: 1
      18: 2
      19: 3
    Integral:
      1: 1.01
      2: 1.0
      3: 0.99
      4: 1.0
      5: 1.01
      6: 2.02
      7: 0.96
      8: 1.01
      9: 2.01
      10: 1.01
      11: 1.01
      12: 2.02
      13: 1.0
      14: 2.13
      15: 3.05
      16: 0.99
      17: 0.97
      18: 2.02
      19: 6.05
    J's:
      1: 7.56, 7.56
      2: 1.51, 1.51, 7.70
      3: 1.88
      4: 0.95, 8.11
      5: 1.54, 1.54, 7.51
      6: .nan
      7: 2.29, 8.12
      8: 2.24
      9: .nan
      10: 9.03, 9.03
      11: 6.78, 9.06
      12: 6.21, 6.21
      13: 6.13, 8.93, 15.13
      14: .nan
      15: .nan
      16: 5.61, 16.61
      17: .nan
      18: .nan
      19: .nan
    Method:
      1: Sum
      2: Sum
      3: Sum
      4: Sum
      5: Sum
      6: Sum
      7: Sum
      8: Sum
      9: Sum
      10: Sum
      11: Sum
      12: Sum
      13: Sum
      14: Sum
      15: Sum
      16: Sum
      17: Sum
      18: Sum
      19: Sum
    Name:
      1: F (t)
      2: G (dt)
      3: I (d)
      4: J (dd)
      5: K (dt)
      6: H (s)
      7: L (dd)
      8: M (d)
      9: E (s)
      10: D (t)
      11: C (dd)
      12: B (t)
      13: A (ddd)
      14: N (m)
      15: O (s)
      16: P (dd)
      17: Q (m)
      18: R (m)
      19: S (s)
    Range:
      1: 7.49 .. 7.41
      2: 7.41 .. 7.35
      3: 7.17 .. 7.13
      4: 7.12 .. 7.08
      5: 7.08 .. 7.03
      6: 6.77 .. 6.66
      7: 6.50 .. 6.46
      8: 6.46 .. 6.43
      9: 5.15 .. 5.05
      10: 4.74 .. 4.63
      11: 4.23 .. 4.15
      12: 4.13 .. 4.05
      13: 3.74 .. 3.62
      14: 3.31 .. 3.24
      15: 3.06 .. 3.00
      16: 2.74 .. 2.66
      17: 2.50 .. 2.44
      18: 2.20 .. 2.10
      19: 1.96 .. 1.88
    Shift:
      1: 7.447
      2: 7.385
      3: 7.143
      4: 7.103
      5: 7.061
      6: 6.715
      7: 6.478
      8: 6.456
      9: 5.099
      10: 4.681
      11: 4.189
      12: 4.094
      13: 3.675
      14: 3.28
      15: 3.036
      16: 2.695
      17: 2.476
      18: 2.148
      19: 1.922
  H1_pureshift:
    Annotation: {}
    Area: {}
    Flags: {}
    Impurity/Compound: {}
    Intensity: {}
    Type: {}
    Width: {}
    ppm: {}
  HMBC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
      23: .nan
      24: .nan
      25: .nan
      26: .nan
      27: .nan
      28: .nan
      29: .nan
      30: .nan
      31: .nan
      32: .nan
      33: .nan
      34: .nan
      35: .nan
      36: .nan
      37: .nan
      38: .nan
      39: .nan
      40: .nan
      41: .nan
      42: .nan
      43: .nan
      44: .nan
      45: .nan
      46: .nan
      47: .nan
      48: .nan
      49: .nan
      50: .nan
      51: .nan
      52: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
      14: None
      15: None
      16: None
      17: None
      18: None
      19: None
      20: None
      21: None
      22: None
      23: None
      24: None
      25: None
      26: None
      27: None
      28: None
      29: None
      30: None
      31: None
      32: None
      33: None
      34: None
      35: None
      36: None
      37: None
      38: None
      39: None
      40: None
      41: None
      42: None
      43: None
      44: None
      45: None
      46: None
      47: None
      48: None
      49: None
      50: None
      51: None
      52: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
      23: .nan
      24: .nan
      25: .nan
      26: .nan
      27: .nan
      28: .nan
      29: .nan
      30: .nan
      31: .nan
      32: .nan
      33: .nan
      34: .nan
      35: .nan
      36: .nan
      37: .nan
      38: .nan
      39: .nan
      40: .nan
      41: .nan
      42: .nan
      43: .nan
      44: .nan
      45: .nan
      46: .nan
      47: .nan
      48: .nan
      49: .nan
      50: .nan
      51: .nan
      52: .nan
    Intensity:
      1: 13.57
      2: 13.86
      3: 5.35
      4: 4.65
      5: 6.23
      6: 7.55
      7: 14.9
      8: 3.55
      9: 14.06
      10: 2.78
      11: 7.18
      12: 18.56
      13: 49.7
      14: 27.93
      15: 14.3
      16: 9.49
      17: 13.12
      18: 1.4
      19: 5.31
      20: 9.6
      21: 6.15
      22: 4.12
      23: 22.25
      24: 23.87
      25: 38.4
      26: 2.88
      27: 5.64
      28: 9.27
      29: 3.59
      30: 10.76
      31: 1.88
      32: 1.54
      33: 11.84
      34: 25.82
      35: 9.35
      36: 3.75
      37: 9.31
      38: 4.34
      39: 2.32
      40: 4.45
      41: 2.49
      42: 1.61
      43: 3.81
      44: 15.81
      45: 6.29
      46: 6.5
      47: 9.48
      48: 3.16
      49: 4.67
      50: 180.86
      51: 116.06
      52: 133.29
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
      14: Compound
      15: Compound
      16: Compound
      17: Compound
      18: Compound
      19: Compound
      20: Compound
      21: Compound
      22: Compound
      23: Compound
      24: Compound
      25: Compound
      26: Compound
      27: Compound
      28: Compound
      29: Compound
      30: Compound
      31: Compound
      32: Compound
      33: Compound
      34: Compound
      35: Compound
      36: Compound
      37: Compound
      38: Compound
      39: Compound
      40: Compound
      41: Compound
      42: Compound
      43: Compound
      44: Compound
      45: Compound
      46: Compound
      47: Compound
      48: Compound
      49: Compound
      50: Compound
      51: Compound
      52: Compound
    Volume:
      1: 1306.65
      2: 1300.87
      3: 414.31
      4: 333.54
      5: 669.03
      6: 481.13
      7: 1024.42
      8: 288.31
      9: 1340.29
      10: 262.52
      11: 830.54
      12: 1108.71
      13: 3412.41
      14: 1705.65
      15: 988.5
      16: 841.43
      17: 1101.28
      18: 72.4
      19: 626.52
      20: 612.02
      21: 443.97
      22: 242.33
      23: 1350.05
      24: 1318.27
      25: 2381.36
      26: 205.57
      27: 480.8
      28: 863.28
      29: 246.67
      30: 1206.51
      31: 403.94
      32: 130.25
      33: 1001.17
      34: 2107.69
      35: 760.0
      36: 531.05
      37: 491.13
      38: 224.04
      39: 130.92
      40: 274.26
      41: 178.09
      42: 117.68
      43: 324.43
      44: 2589.45
      45: 520.23
      46: 556.92
      47: 653.54
      48: 228.2
      49: 355.8
      50: 11657.2
      51: 7512.08
      52: 8223.64
    Width f1:
      1: 198.71
      2: 191.9
      3: 178.23
      4: 229.41
      5: 266.32
      6: 178.61
      7: 169.18
      8: 177.36
      9: 174.71
      10: 267.78
      11: 203.28
      12: 178.28
      13: 179.28
      14: 178.37
      15: 177.45
      16: 179.01
      17: 178.58
      18: 187.43
      19: 268.6
      20: 178.91
      21: 196.3
      22: 181.47
      23: 178.15
      24: 187.13
      25: 179.08
      26: 181.02
      27: 178.68
      28: 175.41
      29: 180.58
      30: 176.08
      31: 315.82
      32: 196.79
      33: 181.44
      34: 178.37
      35: 179.35
      36: 186.22
      37: 179.52
      38: 177.02
      39: 181.42
      40: 178.73
      41: 179.12
      42: 180.36
      43: 178.67
      44: 182.77
      45: 182.22
      46: 183.4
      47: 188.42
      48: 180.1
      49: 180.15
      50: 182.71
      51: 188.64
      52: 177.79
    Width f2:
      1: 17.24
      2: 17.4
      3: 15.47
      4: 11.12
      5: 14.36
      6: 12.7
      7: 14.47
      8: 16.31
      9: 19.43
      10: 12.56
      11: 20.24
      12: 11.93
      13: 13.63
      14: 12.19
      15: 13.86
      16: 17.63
      17: 16.73
      18: 9.79
      19: 15.64
      20: 12.68
      21: 13.09
      22: 11.55
      23: 12.12
      24: 10.5
      25: 12.33
      26: 14.03
      27: 16.99
      28: 18.9
      29: 13.56
      30: 22.66
      31: 24.25
      32: 15.32
      33: 16.58
      34: 16.28
      35: 16.12
      36: 27.07
      37: 10.46
      38: 10.37
      39: 11.07
      40: 12.26
      41: 14.23
      42: 14.41
      43: 16.97
      44: 31.89
      45: 16.15
      46: 16.62
      47: 13.02
      48: 14.28
      49: 15.04
      50: 12.56
      51: 12.21
      52: 12.35
    f1 (ppm):
      1: 137.862
      2: 140.721
      3: 69.864
      4: 134.473
      5: 126.389
      6: 69.864
      7: 161.137
      8: 37.592
      9: 159.585
      10: 134.392
      11: 126.389
      12: 157.462
      13: 134.473
      14: 21.257
      15: 113.853
      16: 97.48
      17: 122.428
      18: 159.667
      19: 161.137
      20: 107.524
      21: 159.667
      22: 122.51
      23: 159.584
      24: 126.382
      25: 137.813
      26: 122.458
      27: 161.188
      28: 39.619
      29: 161.188
      30: 39.619
      31: 37.572
      32: 122.458
      33: 22.544
      34: 51.046
      35: 157.434
      36: 122.458
      37: 22.544
      38: 65.919
      39: 51.046
      40: 37.572
      41: 173.575
      42: 122.458
      43: 77.646
      44: 173.575
      45: 122.458
      46: 77.569
      47: 37.572
      48: 65.919
      49: 51.046
      50: 137.118
      51: 134.434
      52: 113.822
    f2 (ppm):
      1: 7.448
      2: 7.448
      3: 7.385
      4: 7.145
      5: 7.145
      6: 7.145
      7: 7.106
      8: 7.104
      9: 7.103
      10: 7.063
      11: 7.063
      12: 6.716
      13: 6.716
      14: 6.715
      15: 6.714
      16: 6.48
      17: 6.477
      18: 6.477
      19: 6.458
      20: 6.458
      21: 6.456
      22: 6.456
      23: 5.101
      24: 5.099
      25: 5.099
      26: 4.684
      27: 4.682
      28: 4.681
      29: 4.191
      30: 4.189
      31: 4.188
      32: 4.188
      33: 4.096
      34: 4.096
      35: 4.095
      36: 3.676
      37: 3.283
      38: 3.28
      39: 3.035
      40: 2.699
      41: 2.697
      42: 2.697
      43: 2.696
      44: 2.484
      45: 2.484
      46: 2.484
      47: 2.479
      48: 2.149
      49: 2.146
      50: 1.922
      51: 1.922
      52: 1.922
  HSQC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
      14: None
      15: None
      16: None
      17: None
      18: None
      19: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
    Intensity:
      1: 11.55
      2: 12.92
      3: 24.58
      4: 15.7
      5: 14.35
      6: 60.7
      7: 11.89
      8: 18.69
      9: -36.32
      10: -5.72
      11: -7.95
      12: -26.21
      13: 9.95
      14: -20.7
      15: 45.63
      16: -1.75
      17: -4.46
      18: -14.85
      19: 101.53
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
      14: Compound
      15: Compound
      16: Compound
      17: Compound
      18: Compound
      19: Compound
    Volume:
      1: 376.14
      2: 319.88
      3: 348.02
      4: 336.25
      5: 367.5
      6: 702.07
      7: 331.09
      8: 346.12
      9: -558.4
      10: -59.22
      11: -302.33
      12: -593.06
      13: 323.59
      14: -572.26
      15: 668.6
      16: -33.61
      17: -97.92
      18: -533.77
      19: 1132.58
    Width f1:
      1: 85.02
      2: 83.05
      3: 82.84
      4: 83.63
      5: 83.94
      6: 82.59
      7: 82.52
      8: 81.72
      9: 83.3
      10: 83.97
      11: 84.83
      12: 84.31
      13: 83.12
      14: 83.85
      15: 82.75
      16: 84.09
      17: 78.2
      18: 84.93
      19: 82.99
    Width f2:
      1: 19.95
      2: 15.53
      3: 8.9
      4: 13.33
      5: 15.88
      6: 7.29
      7: 17.57
      8: 11.8
      9: 9.61
      10: 6.42
      11: 23.34
      12: 13.98
      13: 20.38
      14: 17.17
      15: 9.22
      16: 11.88
      17: 14.63
      18: 22.04
      19: 7.0
    f1 (ppm):
      1: 129.061
      2: 126.309
      3: 129.018
      4: 125.034
      5: 129.233
      6: 113.779
      7: 107.468
      8: 97.414
      9: 69.809
      10: 77.605
      11: 77.584
      12: 65.879
      13: 37.557
      14: 51.044
      15: 40.734
      16: 39.636
      17: 39.636
      18: 22.557
      19: 21.209
    f2 (ppm):
      1: 7.448
      2: 7.384
      3: 7.145
      4: 7.103
      5: 7.061
      6: 6.715
      7: 6.478
      8: 6.456
      9: 5.099
      10: 4.682
      11: 4.19
      12: 4.095
      13: 3.676
      14: 3.28
      15: 3.037
      16: 2.696
      17: 2.479
      18: 2.149
      19: 1.922
  NOESY:
    Annotation: {}
    Flags: {}
    Impurity/Compound: {}
    Intensity: {}
    Type: {}
    Volume: {}
    Width f1: {}
    Width f2: {}
    f1 (ppm): {}
    f2 (ppm): {}
  molecule:
    smiles:
      1: CC1=CC(=CC(=C1C2=CC=CC(=C2)COC3=CC4=C(C=C3)[C@@H](CO4)CC(=O)O)C)OCCCS(=O)(=O)C
REF5:
  C13_1D:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
    Area:
      1: 633.83
      2: 764.93
      3: 761.82
      4: 2095.33
      5: 2064.83
      6: 1903.08
      7: 2217.77
      8: 2179.08
      9: 816.24
      10: 2078.75
      11: 2539.95
      12: 2469.32
      13: 1910.2
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
    Intensity:
      1: 98.5
      2: 104.4
      3: 138.1
      4: 372.8
      5: 353.6
      6: 333.4
      7: 310.4
      8: 339.4
      9: 151.3
      10: 360.3
      11: 394.2
      12: 398.2
      13: 298.5
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
    Width:
      1: 1.19
      2: 1.38
      3: 1.09
      4: 1.08
      5: 1.1
      6: 1.14
      7: 1.36
      8: 1.16
      9: 0.98
      10: 1.16
      11: 1.2
      12: 1.14
      13: 1.17
    ppm:
      1: 203.09
      2: 199.28
      3: 151.22
      4: 139.05
      5: 137.72
      6: 134.91
      7: 57.79
      8: 53.72
      9: 52.46
      10: 48.91
      11: 46.36
      12: 26.58
      13: 16.59
  COSY:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
    Intensity:
      1: 0.51
      2: 0.09
      3: 0.98
      4: 0.08
      5: 0.63
      6: 0.08
      7: 0.07
      8: 0.05
      9: 0.09
      10: 0.08
      11: 2.53
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
    Volume:
      1: 112.76
      2: 16.72
      3: 261.8
      4: 12.71
      5: 147.4
      6: 18.94
      7: 16.06
      8: 12.3
      9: 22.01
      10: 20.41
      11: 701.18
    Width f1:
      1: 73.27
      2: 73.05
      3: 84.32
      4: 73.95
      5: 73.23
      6: 72.97
      7: 72.49
      8: 78.18
      9: 74.54
      10: 79.31
      11: 74.28
    Width f2:
      1: 18.21
      2: 15.57
      3: 19.17
      4: 13.33
      5: 19.25
      6: 20.68
      7: 20.64
      8: 20.81
      9: 18.87
      10: 18.89
      11: 22.64
    f1 (ppm):
      1: 1.932
      2: 3.077
      3: 5.98
      4: 3.402
      5: 2.813
      6: 3.075
      7: 1.518
      8: 1.666
      9: 1.518
      10: 1.663
      11: 1.518
    f2 (ppm):
      1: 6.482
      2: 6.077
      3: 6.073
      4: 5.98
      5: 3.403
      6: 3.403
      7: 3.403
      8: 3.403
      9: 3.074
      10: 3.074
      11: 1.663
  H1_1D:
    Class:
      1: dt
      2: dd
      3: dd
      4: dhept
      5: m
      6: d
      7: d
      8: dt
      9: dt
      10: s
    H's:
      1: 1
      2: 1
      3: 1
      4: 1
      5: 1
      6: 1
      7: 3
      8: 1
      9: 1
      10: 3
    Integral:
      1: 0.9
      2: 0.94
      3: 0.95
      4: 0.98
      5: 0.96
      6: 1.0
      7: 3.1
      8: 1.1
      9: 1.13
      10: 3.25
    J's:
      1: 0.99, 0.99, 2.06
      2: 2.92, 5.71
      3: 2.87, 5.71
      4: 1.35, 1.35, 1.41, 1.41, 1.41, 1.41, 4.32
      5: .nan
      6: 3.91
      7: 1.45
      8: 1.53, 1.53, 8.97
      9: 1.78, 1.78, 9.07
      10: .nan
    Method:
      1: Sum
      2: Sum
      3: Sum
      4: Sum
      5: Sum
      6: Sum
      7: Sum
      8: Sum
      9: Sum
      10: Sum
    Name:
      1: A (dt)
      2: B (dd)
      3: C (dd)
      4: D (dhept)
      5: E (m)
      6: F (d)
      7: G (d)
      8: H (dt)
      9: I (dt)
      10: J (s)
    Range:
      1: 6.52 .. 6.44
      2: 6.11 .. 6.04
      3: 6.03 .. 5.95
      4: 3.45 .. 3.37
      5: 3.10 .. 3.05
      6: 2.85 .. 2.78
      7: 1.98 .. 1.90
      8: 1.70 .. 1.64
      9: 1.55 .. 1.49
      10: 1.49 .. 1.41
    Shift:
      1: 6.47
      2: 6.07
      3: 5.98
      4: 3.4
      5: 3.07
      6: 2.81
      7: 1.92
      8: 1.66
      9: 1.52
      10: 1.44
  H1_pureshift:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Area:
      1: 9.57
      2: 14.49
      3: 14.49
      4: 19.41
      5: 24.33
      6: 29.25
      7: 34.17
      8: 20.0
      9: 20.0
      10: 20.0
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Intensity:
      1: 0.5
      2: 0.8
      3: 0.8
      4: 1.1
      5: 1.4
      6: 1.7
      7: 2.0
      8: 1.0
      9: 1.0
      10: 1.0
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
    Width:
      1: 1.76
      2: 1.61
      3: 1.61
      4: 1.46
      5: 1.31
      6: 1.16
      7: 1.01
      8: 1.0
      9: 1.0
      10: 1.0
    ppm:
      1: 6.47
      2: 6.07
      3: 5.98
      4: 3.4
      5: 3.07
      6: 2.81
      7: 1.92
      8: 1.66
      9: 1.52
      10: 1.44
  HMBC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
    Intensity:
      1: 0.01
      2: 0.0
      3: 0.01
      4: 0.07
      5: 0.06
      6: 0.05
      7: 0.06
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
    Volume:
      1: 0.44
      2: 0.03
      3: 0.26
      4: 1.55
      5: 1.74
      6: 1.35
      7: 1.67
    Width f1:
      1: 167.15
      2: 165.86
      3: 165.68
      4: 134.85
      5: 176.95
      6: 166.34
      7: 158.18
    Width f2:
      1: 7.54
      2: 7.31
      3: 7.34
      4: 6.7
      5: 6.76
      6: 6.73
      7: 6.71
    f1 (ppm):
      1: 151.244
      2: 203.233
      3: 139.06
      4: 53.732
      5: 52.471
      6: 203.233
      7: 57.797
    f2 (ppm):
      1: 1.925
      2: 1.925
      3: 1.925
      4: 1.444
      5: 1.444
      6: 1.444
      7: 1.443
  HSQC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
    Intensity:
      1: 0.2
      2: 0.09
      3: 0.09
      4: 0.14
      5: 0.07
      6: 0.15
      7: 0.27
      8: -0.13
      9: -0.08
      10: 0.67
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
    Volume:
      1: 4.52
      2: 4.11
      3: 4.21
      4: 4.93
      5: 1.93
      6: 5.19
      7: 9.11
      8: -7.06
      9: -1.54
      10: 20.52
    Width f1:
      1: 237.67
      2: 235.75
      3: 238.07
      4: 235.26
      5: 236.67
      6: 235.58
      7: 236.39
      8: 237.43
      9: 236.58
      10: 236.75
    Width f2:
      1: 6.36
      2: 13.16
      3: 12.82
      4: 9.75
      5: 7.87
      6: 9.44
      7: 9.4
      8: 15.41
      9: 5.57
      10: 8.41
    f1 (ppm):
      1: 139.073
      2: 137.739
      3: 134.916
      4: 48.932
      5: 53.755
      6: 57.812
      7: 16.636
      8: 46.344
      9: 46.344
      10: 26.634
    f2 (ppm):
      1: 6.475
      2: 6.071
      3: 5.98
      4: 3.401
      5: 3.074
      6: 2.813
      7: 1.926
      8: 1.663
      9: 1.518
      10: 1.445
  NOESY:
    Annotation: {}
    Flags: {}
    Impurity/Compound: {}
    Intensity: {}
    Type: {}
    Volume: {}
    Width f1: {}
    Width f2: {}
    f1 (ppm): {}
    f2 (ppm): {}
  Sheet1: {}
  molecule:
    molecule:
      1: C13H14O2
    smiles:
      1: CC1=CC(=O)C2C3CC(C=C3)C2(C)C1=O
beta-lapachone:
  C13_1D:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
    Area:
      1: 21629.27
      2: 17804.96
      3: 27571.88
      4: 65074.18
      5: 28101.66
      6: 65998.26
      7: 19825.87
      8: 72851.34
      9: 67573.82
      10: 23837.38
      11: 36215.97
      12: 63947.06
      13: 141631.72
      14: 79187.57
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
      14: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
    Intensity:
      1: 8918.8
      2: 4915.1
      3: 7575.0
      4: 26149.1
      5: 9665.4
      6: 25104.9
      7: 4701.0
      8: 26133.9
      9: 27027.0
      10: 5704.2
      11: 12545.6
      12: 24273.8
      13: 57078.5
      14: 26442.8
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
      14: Compound
    Width:
      1: 1.51
      2: 2.36
      3: 2.93
      4: 1.7
      5: 2.04
      6: 1.68
      7: 2.85
      8: 1.79
      9: 1.56
      10: 2.91
      11: 1.94
      12: 1.66
      13: 1.55
      14: 1.98
    ppm:
      1: 179.88
      2: 178.58
      3: 162.05
      4: 134.79
      5: 132.63
      6: 130.67
      7: 130.13
      8: 128.59
      9: 124.07
      10: 112.72
      11: 79.28
      12: 31.6
      13: 26.77
      14: 16.17
  COSY:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
    Intensity:
      1: 715.97
      2: 660.41
      3: 158.69
      4: 815.49
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
    Volume:
      1: 81617.38
      2: 72766.72
      3: 11817.83
      4: 45731.27
    Width f1:
      1: 72.92
      2: 70.83
      3: 51.46
      4: 72.3
    Width f2:
      1: 15.11
      2: 15.03
      3: 13.98
      4: 7.49
    f1 (ppm):
      1: 7.509
      2: 7.649
      3: 7.509
      4: 1.86
    f2 (ppm):
      1: 8.065
      2: 7.817
      3: 7.648
      4: 2.577
  H1_1D:
    Class:
      1: m
      2: m
      3: m
      4: d
      5: m
      6: m
      7: m
    H's:
      1: 1
      2: 1
      3: 1
      4: 1
      5: 2
      6: 2
      7: 3
    Integral:
      1: 0.78
      2: 0.76
      3: 0.91
      4: 0.82
      5: 2.01
      6: 2.0
      7: 5.59
    J's:
      1: .nan
      2: .nan
      3: .nan
      4: 2.53
      5: .nan
      6: .nan
      7: .nan
    Method:
      1: Sum
      2: Sum
      3: Sum
      4: Sum
      5: Sum
      6: Sum
      7: Sum
    Name:
      1: A (m)
      2: B (m)
      3: C (m)
      4: D (d)
      5: E (m)
      6: F (m)
      7: G (m)
    Range:
      1: 8.11 .. 8.02
      2: 7.85 .. 7.79
      3: 7.71 .. 7.61
      4: 7.56 .. 7.48
      5: 2.64 .. 2.51
      6: 1.91 .. 1.79
      7: 1.52 .. 1.45
    Shift:
      1: 8.073
      2: 7.817
      3: 7.665
      4: 7.522
      5: 2.584
      6: 1.841
      7: 1.466
  H1_pureshift:
    Annotation: {}
    Area: {}
    Flags: {}
    Impurity/Compound: {}
    Intensity: {}
    Type: {}
    Width: {}
    ppm: {}
  HMBC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
      8: None
      9: None
      10: None
      11: None
      12: None
      13: None
      14: None
      15: None
      16: None
      17: None
      18: None
      19: None
      20: None
      21: None
      22: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
      8: .nan
      9: .nan
      10: .nan
      11: .nan
      12: .nan
      13: .nan
      14: .nan
      15: .nan
      16: .nan
      17: .nan
      18: .nan
      19: .nan
      20: .nan
      21: .nan
      22: .nan
    Intensity:
      1: 12.49
      2: 10.04
      3: 5.38
      4: 3.73
      5: 13.34
      6: 11.95
      7: 12.68
      8: 10.46
      9: 17.03
      10: 7.82
      11: 46.95
      12: 43.86
      13: 70.77
      14: 15.41
      15: 39.77
      16: 42.77
      17: 50.4
      18: 55.91
      19: 49.5
      20: 300.96
      21: 263.07
      22: 370.21
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
      8: Compound
      9: Compound
      10: Compound
      11: Compound
      12: Compound
      13: Compound
      14: Compound
      15: Compound
      16: Compound
      17: Compound
      18: Compound
      19: Compound
      20: Compound
      21: Compound
      22: Compound
    Volume:
      1: 374.44
      2: 301.46
      3: 115.5
      4: 86.97
      5: 1005.73
      6: 291.42
      7: 219.05
      8: 264.67
      9: 480.57
      10: 175.71
      11: 2402.01
      12: 1810.43
      13: 2492.61
      14: 315.94
      15: 1170.79
      16: 1012.08
      17: 2686.94
      18: 1335.0
      19: 1338.02
      20: 7388.12
      21: 7199.02
      22: 9283.07
    Width f1:
      1: 119.15
      2: 125.12
      3: 124.94
      4: 124.5
      5: 141.97
      6: 126.22
      7: 130.19
      8: 162.36
      9: 129.31
      10: 117.89
      11: 118.21
      12: 121.83
      13: 122.13
      14: 123.2
      15: 121.75
      16: 122.1
      17: 122.81
      18: 126.04
      19: 118.75
      20: 124.3
      21: 125.52
      22: 121.52
    Width f2:
      1: 8.8
      2: 8.4
      3: 6.01
      4: 6.55
      5: 18.57
      6: 6.76
      7: 4.64
      8: 5.45
      9: 7.63
      10: 6.67
      11: 15.14
      12: 11.85
      13: 10.09
      14: 5.82
      15: 8.46
      16: 6.78
      17: 15.18
      18: 6.63
      19: 7.96
      20: 6.91
      21: 7.63
      22: 7.22
    f1 (ppm):
      1: 132.633
      2: 134.862
      3: 179.924
      4: 162.089
      5: 130.133
      6: 124.12
      7: 132.633
      8: 128.613
      9: 130.2
      10: 124.12
      11: 31.638
      12: 79.31
      13: 112.718
      14: 178.569
      15: 162.096
      16: 79.31
      17: 112.804
      18: 26.806
      19: 16.167
      20: 31.638
      21: 26.887
      22: 79.31
    f2 (ppm):
      1: 8.066
      2: 8.065
      3: 8.065
      4: 7.817
      5: 7.815
      6: 7.649
      7: 7.648
      8: 7.648
      9: 7.509
      10: 7.509
      11: 2.578
      12: 2.578
      13: 2.577
      14: 2.577
      15: 2.577
      16: 1.86
      17: 1.86
      18: 1.859
      19: 1.859
      20: 1.474
      21: 1.474
      22: 1.474
  HSQC:
    Annotation:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
    Flags:
      1: None
      2: None
      3: None
      4: None
      5: None
      6: None
      7: None
    Impurity/Compound:
      1: .nan
      2: .nan
      3: .nan
      4: .nan
      5: .nan
      6: .nan
      7: .nan
    Intensity:
      1: 14.15
      2: 13.71
      3: 7.58
      4: 10.34
      5: -21.98
      6: -26.27
      7: 86.78
    Type:
      1: Compound
      2: Compound
      3: Compound
      4: Compound
      5: Compound
      6: Compound
      7: Compound
    Volume:
      1: 389.35
      2: 371.78
      3: 242.72
      4: 338.72
      5: 639.02
      6: 780.72
      7: 1923.26
    Width f1:
      1: 82.49
      2: 81.9
      3: 81.17
      4: 83.31
      5: 81.98
      6: 81.62
      7: 84.0
    Width f2:
      1: 17.37
      2: 17.24
      3: 20.55
      4: 20.47
      5: 18.46
      6: 18.96
      7: 13.74
    f1 (ppm):
      1: 128.595
      2: 124.082
      3: 134.78
      4: 130.675
      5: 16.191
      6: 31.601
      7: 26.78
    f2 (ppm):
      1: 8.064
      2: 7.817
      3: 7.648
      4: 7.51
      5: 2.577
      6: 1.86
      7: 1.473
  NOESY:
    Annotation: {}
    Flags: {}
    Impurity/Compound: {}
    Intensity: {}
    Type: {}
    Volume: {}
    Width f1: {}
    Width f2: {}
    f1 (ppm): {}
    f2 (ppm): {}
  molecule:
    smiles:
      1: CC1(C)CCC2=C(O1)C1=CC=CC=C1C(=O)C2=O
"""

exampleproblems_df = {}
# with open("test_df_dict.yaml") as f:
#     exampleproblems_dict = yaml.safe_load(f)
#     for k, v in exampleproblems_dict.items():
#         exampleproblems_df[k] = {k2: pd.DataFrame(v2) for k2, v2 in v.items()}

exampleproblems_dict = yaml.safe_load(examplepromlems_yaml_str)
for k, v in exampleproblems_dict.items():
    exampleproblems_df[k] = {k2: pd.DataFrame(v2) for k2, v2 in v.items()}
exampleproblems_names = list(exampleproblems_df.keys())

print(exampleproblems_names)

exampleproblem_name = exampleproblems_names[1]


def create_excel_directory(exampleproblem_name:str)->tempfile.TemporaryDirectory:
    return tempfile.TemporaryDirectory(prefix=f'{exampleproblem_name}_')

def create_excel_file_path(tmpdir:tempfile.TemporaryDirectory)->Path:
    tdirpathname = Path(tmpdir.name)
    return tdirpathname / f'{tdirpathname.name}.xlsx'

def create_tmp_excel_file(exampleproblem_name:str, df_dict:dict):
    tmpdir = create_excel_directory(exampleproblem_name)
    excel_pathname = create_excel_file_path(tmpdir)

    with pd.ExcelWriter(excel_pathname) as writer:
      for k, v in df_dict.items():
          v.to_excel(writer, sheet_name=k)
    return excel_pathname, tmpdir



if __name__ == "__main__":

    excel_path, tmpdir  = create_tmp_excel_file(exampleproblem_name, exampleproblems_df[exampleproblem_name])

    print(excel_path, "\n", tmpdir.name)

    tmpdir.cleanup()
    