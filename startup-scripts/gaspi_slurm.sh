#!/bin/sh

function mysort { for i in ${node_array[@]}; do echo "$i"; done | sort -n; }

program=$*

export GASPI_SOCKET="$SLURM_LOCALID"
export GASPI_MFILE="$HOME/machinefile_$GASPI_SOCKET"

NODES=`srun hostname`

node_array=($NODES)

sorted_array=( $(mysort) )

echo "${sorted_array[@]}" | tr ' ' '\n' > $GASPI_MFILE

export GASPI_MASTER=${sorted_array[0]}

if [ 0 -eq "$SLURM_PROCID" ];then
  export GASPI_TYPE="GASPI_MASTER"
  echo "Master is $GASPI_MASTER"
  echo "<<start application>>"
else
  export GASPI_TYPE="GASPI_WORKER"
fi

$program

rm $GASPI_MFILE
