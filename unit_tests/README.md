# Deduper Unit Test Files

```input.sam``` contains sorted input lines that contain the following information:

Header scheme: {KEEP or TOSS}:{VALID or INVALID umi}:STRAND:UMI

| Alignment Number | Duplicate? | Keep or Toss | Reason |
|---|---|---|---|
| 1 | No | Keep | unique |
| 2 | Yes | Toss | identical to 1 |
| 3 | Yes | Toss | identical to 1 after adjusting for soft-clipping |
| 4 | No | Keep | different position to 1 |
| 5 | No | Keep | different position to 1 after adjusting for soft-clipping|
| 6 | No | Keep | different umi to 1 |
| 7 | No | Keep | different strand to 1 |
| 8 | Yes | Toss | identical to 1, but further down to make sure order doesn't matter |
| 9 | No | Keep | different chromosome than 1 |
| 10 | No | Keep | unique |
| 11 | No | Toss | invalid UMI |
| 12 | Yes | Toss | same as 11 |
| 13 | No | Keep | same as 11 but on other strand |
| 14 | Yes | Toss | same as 13 after minus position adjustment |
| 15 | Yes | Toss | same as 13 after minus position adjustment |
| 16 | Yes | Toss | same as 13 after minus position adjustment |
| 17 | Yes | Toss | same as 9, further down to test ordering |
| 18 | No | Toss | invalid UMI |


When run though ```Hays_deduper.py``` it I should get ```output.sam``` as the resulting file and the following text output:
```
Number of alignments kept: 8
Number of duplicates removed: 8
Number of unknown umis discarded: 2
```
