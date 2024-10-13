# BI624 Deduper - Lab Notebook
Jules Hays
## Part 1 Due: 10/17/24
## Part 2 Due: 

Python version: 3.12

Environments used: package (package version x.x)

---

### 10/12/24
### Initial Setup

I logged into Talapas and copied the whole genome sequencing data into a directory located at the followng path: /projects/bgmp/jkhay/bioinfo/Bi624/Deduper-julhays

The GitHub repo is located online here: https://github.com/julhays/Deduper-julhays


### Part 1 - Pseudocode
Goal: Write up a strategy for writing a Reference Based PCR Removal Tool. Given a sorted sam file of uniquely mapped reads, remove all PCR duplicates (retain only a single copy of each read). Develop a strategy that avoids loading everything into memory.

* Define the problem

PCR is a necessary step in the DNA or RNA-seq library preparation workflow to get a strong enough signal for detection during sequencing. Unfortunately, this PCR amplification can result in having 2 or more copies of a read that came from the same molecule of DNA/RNA isolated from the sample. This can be detrimental in downstream processes such as differential expression analysis or splicing analysis because these analyses can not distinguish between a PCR duplicate or a biologically relevant duplicate where the original sample had multiple transcripts of a gene. Therefore, it is important to remove PCR duplicates before performing these analyses. Deduplication often occurs with the SAM file after aligning the reads to a genome or transcriptome, which is referred to as reference based dedepulication,because then you only have to compare reads aligned in the same position which deceases the search space. In a SAM file, duplicates will have the same chromosome, the same position on the chromosome (when adjusting for soft-clipping), and be on the same strand. You can distinguish between PCR duplicates and biological duplicates as well with something called an UMI, which is a unique molecular index that is randomly added to the library. If reads have the same alignment position and UMI, it is more likely to be a PCR duplicate and all but 1 should be removed. 

The goal of this script is to create a reference based PCR duplicate removal tool for single-end data given a list of 96 known UMIs. This tool inputs a sorted SAM file that has sequencing reads uniquely mapped to a reference genome/transcriptome, and outputs a SAM file that retains only 1 copy of PCR deuplicated reads.


* Write examples:  

Below are the cases that could be encountered within the SAM files:
1. everything same (chromosome, position, strand, UMI) - duplicate, keep 1
2. everything is the same but the soft clipped adjustment makes position the same - duplicates, keep 1
3. everything is the same but different positions - not duplicates, keep both
4. everything is the same but different UMIs - not duplicates, keep both
5. everything is the same but on different strands - not duplicates, keep both
6. everything is the same but on different chromosomes - not duplicates, keep both

These cases and the outcome are accounted for in an input SAM file called ```input.sam``` and an output SAM file called ```ouput.sam``` located in this folder.

* Develop your algorithm using pseudocode

```
1. store input SAM header into a text file called header.txt

2. sort input by UMI then chromosome with bash and grep out headers

3. read in known UMIs into a set, known_UMIs

4. read in header file and write to output file, dedup.sam

5. initialize 2 empty sets, plus_positions and minus_positions, to store position information for plus and minus strand

6. initialize a temp_UMI and temp_chromo variable to track what UMI/chromo combo you are currently working on

7. read in a line of the input file with a while loop:

#check if UMI in list of known UMIs, if its not I'm not really sure what we are supposed to do, I'll have to ask (for now I'm just outputting to a separate file)
if UMI not in known_UMIs:
    continue with the next line
    (or for the bonus, I could calculate the hamming distance of the UMI to all the known UMIs, select the known UMI with the closest hamming distance and replace it with that. But this would need to be done before the lines are read into this)

#check to see if we've moved onto a new UMI/chromo pair
if read in UMI != temp_UMI or if read in chromo != temp_chromo
    set temp_UMI var equal to UMI of line
    set plus_positions equal to empty set
    set minus_positions equal to empty set

#check to see if pos needs to be adjusted for soft clipping
if 'S' on the start of the CIGAR string:
    position = adjust_position_for_soft_clip(line, "minus")
else:
    position = extract position from line

#check what strand
if bit 16 is true in bitwise flag:
    #if position already in set, then this is a duplicate and will be ignored.nIf not it will be outputted
    if position not in minus_positions set:
        add position to minus_positions set
        write line to output file (dedup.sam)

#same process as above but for the plus strand
else:
    if position not in plus_positions set:
        add position to plus_positions set
        write line to output file (dedup.sam)

8. repeat step 7 for all lines in the input file
```

* Determine high level functions  

```
def adjust_position_for_soft_clip(line: str, strand: str) -> int:
    ```Takes a line of a SAM file that has been soft clipped on     the front and adjusts for the soft clipping by subtracting the  number by the 'S' from the current position. Outputs the newly adjusted position.```
    return adjusted_position
Input: NS500451:154:HWKTMBGXX:1:11101:22945:1315:AGTGCTGT	0	2	100 36	4S71M	*	0	0	GGTGTCATAAACAACAGGCTCTTTGCTGTCGGGTTTCTTTGGCGGAGCCTTGAACCAGGTGAGAGTTGGGA	66EEEEEEEE6EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEEE	MD:Z:71	NH:i:1	HI:i:1	NM:i:0	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU
Expected output: 96
```


---


In order to sort the file by UMI then chromosome with bash we will use the following command:
```
cat pseudocode/input.sam | grep -v "@" | sort -t ':'  -k 8 | sort -k 3 -s
```