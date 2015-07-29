# GASPI-Applications

SLURM
=====

interactive execution:
```shell
srun -n 2 /path/to/GASPI-Applications/startup-scripts/gaspi_slurm.sh /path/to/gpi2-app
```

GASPI-UTILS
============
Its a collection of helper functions which I often used to implement GASPI applications.
GASPI-UTILS consists of the following functions and mechanisms:
* segment creation with an internal segment id counter
* deletion of all created segments
* check queue size
* function that returns a queue which has a number of free entries or blocks if there is no suitable queue
* blocking waitsome function, which calls a waitsome and reset
* flushing of queues (in a range)
* Marcos to check errors of GASPI function calls
