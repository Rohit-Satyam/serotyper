f1=$1
f2=$2

{
  # header: keep all cols from file1 + all cols from file2 except the key column
  paste <(head -n1 "$f1") <(head -n1 "$f2" | cut -f2-)

  # full outer join on col1 (requires sorting), fill missing with empty string
  join -t $'\t' -a 1 -a 2 -e '' -o auto -1 1 -2 1 \
    <(tail -n +2 "$f1" | sort -t $'\t' -k1,1) \
    <(tail -n +2 "$f2" | sort -t $'\t' -k1,1)
} > Serotyper_report.tsv
