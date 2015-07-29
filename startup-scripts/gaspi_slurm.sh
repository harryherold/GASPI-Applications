#!/bin/sh

function mysort { for i in ${node_array[@]}; do echo "$i"; done | sort -n; }

program=$*

# requires to variable NODES in batch script

TMP_PATTERN=".gpi2.XXXXXXXX"
export GASPI_SOCKET="$SLURM_LOCALID"
export GASPI_MFILE=$(mktemp --tmpdir=/tmp ${TMP_PATTERN})

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
