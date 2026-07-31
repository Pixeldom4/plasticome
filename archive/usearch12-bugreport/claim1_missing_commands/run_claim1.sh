#!/usr/bin/env bash
# Claim 1: usearch12.0-beta does not implement the exhaustive all-vs-all / e-value commands
# that usearch11.0.667 provides. Full raw output of each invocation is saved under logs/.
set -u
U11=../../../bin/usearch11; U12=../../../bin/usearch12; FA=../inputs/mre_n84_n35.fasta
mkdir -p logs
verdict () {  # $1=logfile -> one-line verdict
  if grep -qi 'unknown command-line option' "$1"; then echo "REJECTED: $(grep -i 'unknown command-line option' "$1" | head -1)"
  else echo "RECOGNIZED (parses command; runs or requests expected args)"; fi
}
printf "%-16s | %-52s | %s\n" "command" "usearch v11.0.667" "usearch v12.0-beta"
printf "%-16s-+-%-52s-+-%s\n" "----------------" "----------------------------------------------------" "----------------------------------------------------"
for cmd in allpairs_local allpairs_global ublast calc_distmx usearch_local usearch_global; do
  $U11 -$cmd $FA -output /dev/null -userout /dev/null >logs/u11_$cmd.log 2>&1 || true
  $U12 -$cmd $FA -output /dev/null -userout /dev/null >logs/u12_$cmd.log 2>&1 || true
  printf "%-16s | %-52s | %s\n" "-$cmd" "$(verdict logs/u11_$cmd.log)" "$(verdict logs/u12_$cmd.log)"
done
