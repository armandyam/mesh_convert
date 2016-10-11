#!/bin/bash
./omega_h --axes axes_1.vtu --adapt_log adapting_1 ~/src/ajay/inputs/uniform_simple/netgen.vol ~/src/ajay/inputs/uniform_simple/metric_file.mtr adapted_1.vol &> log_1
tail log_1
./omega_h --limit 1.1 --axes axes_2.vtu --adapt_log adapting_2 ~/src/ajay/inputs/problematic/netgen.vol ~/src/ajay/inputs/problematic/adj_hess_metric.mtr adapted_2.vol &> log_2
tail log_2
./omega_h --smooth 4 --axes axes_3.vtu --adapt_log adapting_3 ~/src/ajay/inputs/netgen_2888_ep0p005_p2.vol ~/src/ajay/inputs/metric_2888_ep0p005_p2.mtr adapted_2.vol &> log_3
tail log_3
