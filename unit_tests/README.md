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

```test.sam``` is an additional test .sam file that was provided.


# Challenge information

## UMI error correction

```input_umi_correct.sam``` contains 1 additional line compared to the original ```input.sam``` to test the ```Hays_deduper_umicorrect.py``` script. Here's what's different.

| Alignment Number | Duplicate? | Keep or Toss | Reason |
|---|---|---|---|
| 11 | No | Toss | invalid UMI, can't be corrrected |
| 18 | No | Toss | invalid UMI that is corrected and now its a duplicate |
| 19 | No | Keep | invalid UMI that is corrected and not a duplicate |

When run through ```Hays_deduper_umicorrect.py``` with the -e flag is set to True (UMI correction), you should expect the file output to look like ```output_umi_correct.sam``` and get the resulting text output:
```
Number of alignments kept: 9
Number of duplicates removed: 9
Number of unknown umis discarded: 1
Number of unknown umis corrected: 2
```


## Choice of duplicate written out

```input_choice.sam``` contains a couple lines to test the ```Hays_deduper_choice.py``` script.

| Alignment Number | Duplicate? | Keep or Toss | Reason |
|---|---|---|---|
| 1 | Yes | Toss | same as 3 but lower average quality |
| 2 | Yes | Toss | same as 3 but shorter length |
| 3 | No | Keep | highest quality and longest length |

When run though ```Hays_deduper_choice.py``` it I should get ```output_choice.sam``` as the resulting file and the following text output:
```
Number of alignments kept: 1
Number of duplicates removed: 2
Number of unknown umis discarded: 0
```
