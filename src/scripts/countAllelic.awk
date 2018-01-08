#samtools view H3K36me3_fib01_A08393_hap1.bam | awk -f ~/Allele-specific/tools/countAllelic.awk

BEGIN {
  min = 0
  max = 0
}

{
  idx=$5
  if ($3 ~ /C57BL6J/)
    idx =-idx

  cnt[idx]++;

  if (min >= idx)
    min = idx
  if (max <= idx)
    max = idx
}

END {
  for (x = min; x <= max; x++) {
    if (x in cnt)
      print x, cnt[x]
    else
      print x, 0
  }
}
