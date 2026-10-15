#!/bin/csh -f
# Usage: run_job.csh <jobname>
cd /root/gamess
source install.info
set JOB = $1
set VERNO = 00
set GMSPATH = /root/gamess
set JOB_PATH = /root/gamess
set SCR = /root/gamess/scratch
set USERSCR = /root/gamess/restart
rm -f $SCR/$JOB.F* $USERSCR/$JOB.*
mkdir -p $SCR $USERSCR
cp -f $JOB.inp $SCR/$JOB.F05
source gms-files.csh
$GMSPATH/gamess.$VERNO.x "$JOB -scr $SCR < /dev/null" > /root/gamess/$JOB.out
echo "RUN_EXIT=$status"
