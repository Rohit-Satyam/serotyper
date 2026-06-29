seqkit fx2tab --threads $2 -l $1 | awk '
BEGIN { sum_log = 0; count = 0 }
{
    sum_log += log($4)
    count++
}
END { printf "%d\n", exp(sum_log / count) + 0.5 }
'
