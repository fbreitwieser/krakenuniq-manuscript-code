# HyperLogLog test

The HyperLogLog implementation was tested both on random numbers and actual database k-mers. While filling the HyperLogLog counters, true counts and intermediate estimates are outputed at factors of 10 and tenths of factors of 10 - i. e. 1, 10, 100, 110, 120, ..., 1000, 1100, 1200, ... . The binaries to generate and count the numbers are part of KrakenHLL, but are not built by default. To build them, call `make allall` in the `src/` directory.

## Results on random data

We generated 100 times 10 million random numbers using a 64-bit Mersenne Twister 19937 as random number generator, seeded with C++ std::random_device. As we have 100 estimates at each true count, we can plot the variance of the relative error of the estimates.

Command line to generate files:
```{sh}
for P in 10 12 14 16 18; do ./count_unique -p $P -r 10000000 -s -t -y -f -e -x 100 > random/hll-p$P.txt; done
```

Head of file:
```
observed	estimate_heule	estimate_flajolet	estimate_ertl	rel_error_heule	rel_error_flajolet	rel_error_ertl	who_won
1	1	1	1	0	0	0	equal
10	10	10	10	0	0	0	equal
100	100	104	100	0	0.04	0	equal
110	110	115	110	0	0.0454545	0	equal
120	120	124	120	0	0.0333333	0	equal
130	130	134	130	0	0.0307692	0	equal
140	140	146	140	0	0.0428571	0	equal
150	150	156	150	0	0.04	0	equal
160	160	167	160	0	0.04375	0	equal
```

## Results on real k-mers

We sampled about one million k-mers from the RefSeq database described in the manuscript and estimated their counts with HyperLogLog counters with precision between 10 and 18. Each k-mer was selected with a probability of 1.01 million / # of k-mers in the database.

Command line to generate the file:
```{sh}
./test_hll_on_db ../dbs/refseq-oct2017-k31/database.kdb 1 1010000
```

Head of file:
```
precision	true_count	estimate	ertl_estimate
10	1	1	1
11	1	1	1
12	1	1	1
13	1	1	1
14	1	1	1
15	1	1	1
16	1	1	1
17	1	1	1
18	1	1	1
```

## Plotting of results

Plots were generated using R and ggplot2, see the R script for details.
