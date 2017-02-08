#Introduction
The is page where I document output of Assessment Metric script for testing purpose when modifying the script to select the optimal path if >1 max Chr Id occur.
#Feb 8 2017
The below result was from 3909 test run with testfile1m:
```
TRINITY_DN44051_c0_g1_i1        Chr02 Chr06 Chr08
TRINITY_DN44051_c0_g1_i1        fragment1       Chr06   52      254     108960163       108959961       1.16e-91        340 -  1       203
TRINITY_DN44051_c0_g1_i1        fragment1       Chr06   52      254     108976720       108976518       1.16e-91        340 -  1       203
TRINITY_DN44051_c0_g1_i1        fragment1       Chr02   52      254     58737657        58737455        1.16e-91        340 -  1       203
TRINITY_DN44051_c0_g1_i1        fragment2       Chr08   1       51      77623059        77623009        2.73e-17        93.3-  1       51
TRINITY_DN44051_c0_g1_i1        1       Problem: different chromosome ID
TRINITY_DN44051_c0_g1_i1        2       Problem: different chromosome ID
TRINITY_DN44051_c0_g1_i1        3       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 scaffold_114 scaffold_51 scaffold_875
TRINITY_DN49775_c0_g1_i1        fragment1       Chr01   43      660     126383056       126382440       0.0     1097    -1  6  18
TRINITY_DN49775_c0_g1_i1        fragment2       Chr01   1       48      37255673        37255720        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr01   1       48      60710484        60710437        3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr03   1       48      122145168       122145215       3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr03   1       48      132294580       132294627       3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr03   1       48      132717911       132717864       3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr09   1       48      71193185        71193138        3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr07   1       48      34530173        34530220        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr07   1       48      43556321        43556368        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr07   1       48      48891084        48891131        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr05   1       48      35314431        35314384        3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr05   1       48      52376881        52376834        3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr05   1       48      70713901        70713948        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr01   1       48      113097210       113097257       3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr05   1       48      135475715       135475668       3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr04   1       48      581103  581150  3.25e-15        87.8    1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr04   1       48      25742320        25742367        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr04   1       48      52262706        52262753        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr04   1       48      101973268       101973221       3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr04   1       48      123273951       123273998       3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr02   1       48      6402395 6402348 3.25e-15        87.8    -1      48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr02   1       48      20156261        20156308        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr02   1       48      43436032        43435985        3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       scaffold_875    1       48      6380    6427    3.25e-15        87.8    1   4  8
TRINITY_DN49775_c0_g1_i1        fragment2       Chr01   1       48      174662817       174662864       3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       scaffold_51     1       48      913747  913700  3.25e-15        87.8    -1  4  8
TRINITY_DN49775_c0_g1_i1        fragment2       scaffold_114    1       48      273940  273987  3.25e-15        87.8    1   4  8
TRINITY_DN49775_c0_g1_i1        fragment2       Chr08   1       48      63347399        63347352        3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr08   1       48      105463331       105463378       3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr06   1       48      65631879        65631926        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr03   1       48      39025478        39025525        3.25e-15        87.81  48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr03   1       48      76319023        76318976        3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        fragment2       Chr03   1       48      110622514       110622467       3.25e-15        87.8-  1       48
TRINITY_DN49775_c0_g1_i1        1       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        2       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        3       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        4       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        5       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        6       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        7       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        8       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        9       Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        10      Problem: different chromosome ID
TRINITY_DN49775_c0_g1_i1        11      Problem: different chromosome ID
TRINITY_DN53613_c0_g2_i1        Chr05 Chr08
TRINITY_DN53613_c0_g2_i1        fragment1       Chr08   1       254     59463767        59463514        8.28e-126       453 -  1       254
TRINITY_DN53613_c0_g2_i1        fragment3       Chr05   249     285     91317657        91317621        6.40e-07        59.0-  1       37
TRINITY_DN53613_c0_g2_i1        1       Problem: different chromosome ID
TRINITY_DN53613_c0_g2_i1        2       Problem: different chromosome ID
TRINITY_DN53617_c0_g1_i1        Chr02 Chr08
TRINITY_DN53617_c0_g1_i1        fragment1       Chr08   1       216     60080255        60080470        1.03e-106       390 1  216
TRINITY_DN53617_c0_g1_i1        fragment3       Chr02   216     347     60769063        60769194        4.68e-54        215 1  132
TRINITY_DN53617_c0_g1_i1        1       Problem: different chromosome ID
TRINITY_DN53617_c0_g1_i1        2       Problem: different chromosome ID

```
