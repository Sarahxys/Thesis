# Introduction
The is page where I document output of Assessment Metric script for testing purpose when modifying the script to select the optimal path if >1 max Chr Id occur.

# Percent of transcripts with >1MaxChrID
* Tropicalis: 6%
* Laevis: 10%

# Feb 8 2017
* The below result was from 3909 test run with testfile1m. 
* Problem observed: the second transcript should not be marked as having more than two maxChr as only Chr01 occurs in both fragments.
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
* Problem identified in script: the script checked $max after ChrID, fragID, and tightID were stored in temptight_hash while it should check $max and update $maxChrID everytime after adding a new tight/fragment. 
* Solution: The problem were solved by putting the check for $max/$maxChrID in the loop and below are the print out. Problem solved.
```
the max number is 1 and the max ChrID is Chr06
Chr02   fragment1
Chr06   fragment1
Chr08   fragment2
TRINITY_DN44051_c0_g1_i1        Chr02 Chr06 Chr08
TRINITY_DN44051_c0_g1_i1        fragment1       Chr06   52      254     108960163       108959961       1.16e-91        340
TRINITY_DN44051_c0_g1_i1        fragment1       Chr06   52      254     108976720       108976518       1.16e-91        340
TRINITY_DN44051_c0_g1_i1        fragment1       Chr02   52      254     58737657        58737455        1.16e-91        340
TRINITY_DN44051_c0_g1_i1        fragment2       Chr08   1       51      77623059        77623009        2.73e-17        93.3
TRINITY_DN44051_c0_g1_i1        1       Problem: different chromosome ID
TRINITY_DN44051_c0_g1_i1        2       Problem: different chromosome ID
TRINITY_DN44051_c0_g1_i1        3       Problem: different chromosome ID
the max number is 1 and the max ChrID is Chr08
Chr05   fragment3
Chr08   fragment1
TRINITY_DN53613_c0_g2_i1        Chr05 Chr08
TRINITY_DN53613_c0_g2_i1        fragment1       Chr08   1       254     59463767        59463514        8.28e-126       453
TRINITY_DN53613_c0_g2_i1        fragment3       Chr05   249     285     91317657        91317621        6.40e-07        59.0
TRINITY_DN53613_c0_g2_i1        1       Problem: different chromosome ID
TRINITY_DN53613_c0_g2_i1        2       Problem: different chromosome ID
the max number is 1 and the max ChrID is Chr08
Chr02   fragment3
Chr08   fragment1
TRINITY_DN53617_c0_g1_i1        Chr02 Chr08
TRINITY_DN53617_c0_g1_i1        fragment1       Chr08   1       216     60080255        60080470        1.03e-106       390
TRINITY_DN53617_c0_g1_i1        fragment3       Chr02   216     347     60769063        60769194        4.68e-54        215
TRINITY_DN53617_c0_g1_i1        1       Problem: different chromosome ID
TRINITY_DN53617_c0_g1_i1        2       Problem: different chromosome ID

```
# Feb 09 2017
* problem 1
-> fragment5, 8, and 9 shouldn't co-exist/selected
-> below modification was made to correct this error
```
elsif((($tqstart = $tempqstart)&&($tqend < $tempqend))||
		   		(($tqstart > $tempqstart)&&($tqend = $tempqend))||
		   		(($tqstart > $tempqstart)&&($tqend < $tempqend))){
#fragment comparison case #4: temp cover the region of the existing fragment, store the id of the fragment in case it is long fragment that overlap with 3+ existing fragment -> compare e-value
```
```
the max number is 1 and the max ChrID is Chr09
Chr05   fragment9
Chr06   fragment5
Chr07   fragment8
Chr09   fragment1
TRINITY_DN34746_c0_g1_i1        Chr05 Chr06 Chr07 Chr09
TRINITY_DN34746_c0_g1_i1        fragment1       Chr09   96      218     26869622        26869743        1.18e-52        210     1       123
TRINITY_DN34746_c0_g1_i1        fragment5       Chr06   1       71      12496497        12496427        1.65e-25        120     -1      71
TRINITY_DN34746_c0_g1_i1        fragment8       Chr07   1       76      79934849        79934774        1.65e-25        120     -1      76
TRINITY_DN34746_c0_g1_i1        fragment9       Chr05   1       85      96479241        96479156        1.65e-25        120     -1      86
TRINITY_DN34746_c0_g1_i1        1       Problem: different chromosome ID
TRINITY_DN34746_c0_g1_i1        2       Problem: different chromosome ID
TRINITY_DN34746_c0_g1_i1        3       Problem: different chromosome ID
TRINITY_DN34746_c0_g1_i1        4       Problem: different chromosome ID

```
* Problem 2
-> cosider as error?
```
the max number is 1 and the max ChrID is scaffold_42
scaffold_113    fragment3
scaffold_314    fragment4
scaffold_42     fragment1
TRINITY_DN34783_c8_g3_i5        scaffold_113 scaffold_314 scaffold_42
TRINITY_DN34783_c8_g3_i5        fragment1       scaffold_42     1       190     728923  729113  1.44e-87        327     1       191
TRINITY_DN34783_c8_g3_i5        fragment3       scaffold_113    187     229     418849  418891  5.05e-11        73.4    1       43
TRINITY_DN34783_c8_g3_i5        fragment4       scaffold_314    185     229     82186   82142   5.05e-11        73.4    -1      45
TRINITY_DN34783_c8_g3_i5        1       Problem: different orientation
TRINITY_DN34783_c8_g3_i5        2       Problem: different orientation
TRINITY_DN34783_c8_g3_i5        3       Problem: different orientation
```

* Problem 3
-> correction made: added the below when checking chromosome error.
```
	if ($firstfragname =~ "scaffold"){
		$error_hash{$qseqID} = "Contain scaffold";
	}
```
```
the max number is 1 and the max ChrID is scaffold_1795
Chr05   fragment1
scaffold_1795   fragment1
TRINITY_DN34791_c2_g1_i1        Chr05 scaffold_1795
TRINITY_DN34791_c2_g1_i1        fragment1       scaffold_1795   8       540     4398    3891    0.0     830     -1      533
TRINITY_DN34791_c2_g1_i1        fragment1       Chr05   8       540     144380469       144380972       0       803     1       533
TRINITY_DN34791_c2_g1_i1        1       no error
TRINITY_DN34791_c2_g1_i1        2       no error
path selected:1 and error message is no error

```
* Problem 4
-> correct for it? stated it as one of the limitation?
```
the max number is 1 and the max ChrID is Chr08
Chr03   fragment3
Chr05   fragment3
Chr07   fragment3
Chr08   fragment1
TRINITY_DN34747_c0_g1_i1        Chr03 Chr05 Chr07 Chr08
TRINITY_DN34747_c0_g1_i1        fragment1       Chr08   51      599     38585133        38585681        0.0     991     1       549
TRINITY_DN34747_c0_g1_i1        fragment3       Chr07   1       50      7755254 7755303 2.41e-16        91.5    1       50
TRINITY_DN34747_c0_g1_i1        fragment3       Chr05   1       50      61719614        61719663        2.41e-16        91.5    1       50
TRINITY_DN34747_c0_g1_i1        fragment3       Chr03   1       50      17834069        17834020        2.41e-16        91.5    -1      50
TRINITY_DN34747_c0_g1_i1        1       Problem: different chromosome ID
TRINITY_DN34747_c0_g1_i1        2       Problem: different chromosome ID
TRINITY_DN34747_c0_g1_i1        3       Problem: different chromosome ID
TRINITY_DN34747_c0_g1_i1        4       Problem: different chromosome ID

```

