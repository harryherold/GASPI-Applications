#!/bin/sh

function mysort { for i in ${node_array[@]}; do echo "$i"; done | sort -n; }

program=$*

export GASPI_SOCKET="$SLURM_LOCALID"
export GASPI_MFILE="$HOME/machinefile_$SLURM_PROCID"

NODES=`srun hostname`

node_array=($NODES)

sorted_array=( $(mysort) )
export GASPI_MASTER=${sorted_array[0]}

if [ 0 -eq "$GASPI_SOCKET" ];then
  echo "${sorted_array[@]}" | tr ' ' '\n' > $GASPI_MFILE
fi

if [ 0 -eq "$SLURM_PROCID" ];then
  export GASPI_TYPE="GASPI_MASTER"
  echo "Master is $GASPI_MASTER"
  echo "<<start application>>"
else
  export GASPI_TYPE="GASPI_WORKER"
  echo "<<start worker $SLURM_PROCID>>"
fi

$program

if [ 0 -eq "$GASPI_SOCKET" ];then
  rm $GASPI_MFILE
fi
