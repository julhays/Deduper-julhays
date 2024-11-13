# BI624 Deduper - Lab Notebook
Jules Hays
## Part 1 Due: 10/17/24
## Part 2 Due: 

Python version: 3.12

Packages used: samtools 1.20

File directory:
* ```unit_tests/``` -> input and output unit test encompassing all scenarios
* ```Hays_deduper.py``` -> deduplicating script
* ```STL96.txt``` -> list of known UMIs
* ```lab_notebook.md``` -> notebook for this assignment
* ```part1_pseudocode.md``` -> my intial psuedocode approach for tackling this problem
* ```run_dedup.sh``` -> sbatch file to sort input sam file and run ```Hays_deduper.py```

This challenge branch also includes:
* ```Hays_deduper_umicorrect.py``` -> a new script that incorperates an option for UMI error correction of known UMIs
* ```Hays_deduper_choice.py``` -> a new script version that chooses which duplicate to write out based on what alignment has the highest quality and if the quality is same then the longest length
*```Hays_deduper_randomers.py``` -> a new script that incorperates an option to use 8 nt randomers instead of known UMIs
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

---

### 10/18/24
### Part 2 - Pseudocode peer review

Hi Varsheni,

I like how you outlined all the cases you might encounter in your unit test folder. That is very helpful to see all the possibilities your code will account for. Your flow chart is really pretty as well.

Your algorithm makes a lot of sense to me and I can follow the logic clearly. Sorting by UMI and chromosome is a great idea! I like your idea of using a set to store the unique read information. It is also good that you clear it each time you encounter a new chromosome and UMI, which will save you memory. How do you plan to check if the S is on the start or the end of the CIGAR string?

I believe your algorithm does everything it's supposed to do outlined in part 1. Good job!

Your functions all seem reasonable as well. I like how you included multiple examples of inputs and outputs to thoroughly show what your function is going to do. Including a function to extract chromosome info, as well as the UMI and strand is useful so I may implement that in my code.

Overall, I think you have a great strategy. Good luck coding this up!

-Jules


Hi Zach,

I like how you outlined all the duplicate scenarios you might encounter in input_sam_info.md Links to an external site.. That is very helpful to see all the possibilities your code will account for.

Your problem statement could have a bit more detail, such as why we want to remove duplicates or how we tell what is a duplicate. I lost points on a previous assignment for not having a detailed about problem statement so I don't want the same to happen to you!

Your algorithm makes a lot of sense to me and I can follow the logic clearly. I like your idea of using a set to store the unique read information. It is also good that you clear it each time you encounter a new chromosome, which will save you memory. One suggestion is you could incorporate the the cigar_pos into step 8 so that you can adjust the position if needed before adding it to the tuple rather than having to go back into the tuple and readjust it (tuples are immutable so it would cause you to have to create a whole new tuple if you did this, which would take time). On that same note, maybe store all the pieces of information from steps 6, 7, 8, and 9 as a local variable and then add them all to a tuple at the same time to avoid having to remake the tuple each time you add one thing.

I believe your algorithm does everything it's supposed to do outlined in part 1. Good job!

Your functions all seem reasonable as well. I like how you included multiple examples of input and output to give a more thorough overview of what your function will do. Maybe you could also include a function to extract the UMI information from the first column?

Overall, I think you have a great strategy. Good luck coding this up!

-Jules



Hi Mahmoud,

Your problem statement is very detailed so I can clearly understand what you are trying to do. One suggestion here is when you say "PCR duplicates should map to the same 5' start of read" that is not always true because the POS field gives the position of the leftmost end of the read (which is the 5' for the + strand and the 3' end for the - strand).

Your algorithm makes a lot of sense to me and I can follow the logic clearly. I like your idea of using a set to store the unique read information. It is also good that you clear it each time you encounter a new chromosome, which will save you memory. The only piece that is confusing to me is the "adjust the start position based on the strand direction". I could be wrong, but I believe the adjusted start position is just the POS number - the number of clipped nts (the number in front of the S in the CIGAR string). Does the strandedness change how you deal with this? One other suggestion is since you are going to have a new set for each chromosome, you probably don't need to store the chromosome information in each individual tuple. You could just have it as a global variable instead.

I believe your algorithm does everything it's supposed to do outlined in part 1. Good job!

Your functions all seem reasonable as well. The find_UMI one is a good idea because I did not consider how I will extract the UMI from the NAME field. Will the find_pos function make the adjustment to the POS based on soft clipping?

Overall, I think you have a great strategy. Good luck coding this up!

-Jules


---
### 10/24/24
### The feedback I recieved

Amelia: Overall, great logic that will capture every case. I like that you clear memory at both a new UMI and a new chromosome - this is code that won't use a ton of memory!

Fantastic overview of the project goal. You obviously have a great grasp of the context of this project and why this script should be used. I also like the different test cases you outlined. You kept it simple to follow what we're looking for - great job.

I'm not sure why you would want to store the header in a separate file - I think you would want to just write all the header lines directly to your single de-duplicated output file without that middle step.

I'm not sure you can sort by UMI since the UMI is all the way at the end of the QNAME column. You might only be able to sort by chromosome and go from there. I do really like the idea of working solely on one UMI and then one chromosome - but it might not be feasible with the sort command. If I'm wrong or you've tested it, then this is a great approach!!

I like outputting invalid UMIs to a separate file, very helpful for the end user. Hamming Distance comment is a compelling error-correction method! That would be cool to implement.

There's some weirdness with correcting a "-" strand that has soft-clipping that I don't see outlined here. That is a complicated concept so if you have further questions about that, let me know!


Claire: Hi Jules,

I really appreciate how thorough your pseudocode is and it is very easy to read. I can easily follow the logic! I especially appreciate how in depth your problem statement was and how you thoroughly accounted for all the different scenarios. I was wondering if you had considered what to do if the S on the CIGAR string is not at the beginning and is actually at the end. You could maybe account for this by maybe taking in CIGAR string as a string and upon hitting an S, stopping and int() -ing everything prior to the S. If it is soft clipped on the left, then maybe you could use an if-else statement to process it and if the int() fails then you could just move on. Just a suggestion! I like your function for calculating the position based on soft clipping. Could you also include a function for strandedness? This could be helpful.

In general, you did a really good job and your code is very thorough and very easy to follow.

Good luck on writing your code!

Claire


The logic here is great, clear and easy to follow, I really like the idea of sorting the initial SAM file by UMI and then chromosome I think that's a smart way of parsing through the file and making the python simpler I wouldn't change anything about your information flow or the way that you're storing and writing out the info.

The one note i have it makes sure to do the extra soft clipping steps for the - strand, to add the M match number, right side soft clipping, and deletions (and then -1 to adjust for 1 based counting) instead of just sticking to the same + strand left side S check. I think that will be a super easy adjustment though.

Also! I'm pretty sure that if the UMI is not in the set of known UMIs we are just disregarding the entire read - this could be different for a challenge question though.

I think that the one function you have here is totally reasonable, everything else makes a lot of sense to be wrapped in the actual loop so that looks great - also love the input output example I will be stealing that formatting!


### Part 3 - Lets start coding

I made a top level file called ```Hays_deduper.py```.

I started writing some of the code body and functions.

---
### 10/28/24
### Part 3 cont

Still coding up ```Hays_deduper.py```.

Questions:
* should unknown UMIs be discarded?
get rid of them but count them
* how to calculate memory usage?
max rss 
usr bin time output of other command in talapas lecture to look it up with job id



I made a script, test.py, to test some regex stuff with my cigar string.
```
chmod 755 test.py
```

---
### 11/2/24
### Part 3 cont

I'm ready to start testing and debugging my script.

```
chmod 755 Hays_deduper.py
```

I fixed lots of errors, then had to adjust some functions so they worked right.

This still isn't working properly, but I suspect my input file is messed up. Knowing what I know now, I am going to redo ```input.sam```.


Description of unit test:

'''input.sam''' contains input lines that contain the following information:

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


When run though ```Hays_deduper.py``` it I should get '''output.sam''' as the resulting file and the following text output:
```
Number of alignments kept: 8
Number of duplicates removed: 8
Number of unknown umis discarded: 2
```

Ok, I am getting these numbers now. Let me check that my outputs are the same:

```
diff dedup_out.sam part1-pseudocode/output.sam
```
This returned nothing, so the files match.

Now, I need to set up the argparse. I will also write an sbatch script called ```run_dedup.sh``` to add all the arguments and run the file.


For the google form, I need to run ```/projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam```. This file does not appear to be sorted, so I added the following sort command to my sbatch script:
```
samtools sort -O sam /projects/bgmp/shared/deduper/C1_SE_uniqAlign.sam -o C1_SE_uniqAlign_sorted.sam
```

I will run this file through my script and I need to report the number of header lines, number of unique reads, number of wrong UMIs, and number of removed duplicates. I also need to report the number of reads per chromosome in the form <chrom_name><tab><count>. Finally, I need to report infomation about how much memory my script used to deduplicate and how long it ran for.

I ran the file:
```
sbatch run_dedup.sh 
Submitted batch job 23107673
```

It errored on the first run because I made the chromosome information an integer but mitochondrial MT "chromosome" can't be an int so I will remove that int casting because it is not necessary.

Running again:
```
sbatch run_dedup.sh 
Submitted batch job 23107681

Number of alignments kept: 13719048
Number of duplicates removed: 4467362
Number of unknown umis discarded: 0
```

It successfully ran and the output file is called ```dedup_out.sam```

Here are all the bash commands I used to obtain / where I found that information:
```
$ srun -A bgmp -p bgmp -N 1 -c 1 -t 1:00:00 --pty bash

# Header lines
$ cat dedup_out.sam | grep "^@" | wc -l
65

# Unique reads - output information
Number of alignments kept: 13719048

# Wrong UMIs - output information
Number of unknown umis discarded: 0

#Duplicates removed - output information
Number of duplicates removed: 4467362

# Number of reads per chromosome
$ cat dedup_out.sam | grep -v "^@" | cut -f 3 | uniq -c | awk -v OFS='\t' '{print $2, $1}' | sort -n

GL456210.1      5
GL456211.1      6
GL456212.1      4
GL456221.1      4
GL456233.2      656
GL456239.1      1
GL456354.1      1
GL456367.1      3
GL456368.1      3
GL456370.1      21
GL456379.1      2
GL456382.1      1
GL456383.1      1
GL456389.1      1
GL456390.1      1
GL456396.1      17
JH584295.1      111
JH584299.1      3
JH584304.1      294
MT      202002
MU069434.1      3
MU069435.1      5450
X       317853
Y       2247
1       697508
2       2787018
3       547615
4       589839
5       562160
6       510818
7       1113183
8       576463
9       627488
10      564903
11      1220389
12      359951
13      467659
14      387239
15      437465
16      360923
17      517566
18      290506
19      571665

#Memory (in GB) - /usr/bin/time info
Maximum resident set size (kbytes): 557364
557364 kbytes = 0.56GB

#Run time (H:mm:ss) - /usr/bin/time info
Elapsed (wall clock) time (h:mm:ss or m:ss): 1:03.11
```

Final report:
```
Header lines: 65
Unique reads: 13719048
Wrong UMIs: 0
Duplicates removed: 4467362
Number of reads per chromosome:

GL456210.1      5
GL456211.1      6
GL456212.1      4
GL456221.1      4
GL456233.2      656
GL456239.1      1
GL456354.1      1
GL456367.1      3
GL456368.1      3
GL456370.1      21
GL456379.1      2
GL456382.1      1
GL456383.1      1
GL456389.1      1
GL456390.1      1
GL456396.1      17
JH584295.1      111
JH584299.1      3
JH584304.1      294
MT      202002
MU069434.1      3
MU069435.1      5450
X       317853
Y       2247
1       697508
2       2787018
3       547615
4       589839
5       562160
6       510818
7       1113183
8       576463
9       627488
10      564903
11      1220389
12      359951
13      467659
14      387239
15      437465
16      360923
17      517566
18      290506
19      571665

Memory (in GB): 0.56
Run time (H:mm:ss): 0:01:03
```

---
### 11/3/24
### Submitting

I submitted the Qualatrics form on the canvas assignment.

I made a little modification to my code that allows for any CIGAR strand to be inputted (including P and H) even though we don't do anything with those characters at least it won't error out. I also added a line that catches any additional characters that aren't valid in a CIGAR string or that I didn't account for and prints a helpful message to let the user know.

I will make sure everything is pushed to GitHub for the base assignment.

I might try to tackle some of the bonuses.

---
### 11/12/24
### Challenge problems - UMI error correction

I made a new branch on my GitHub called bonus. I am going to make a copies of my script to tackle the bonus problems. 

I will start with error correction of UMIs, with a script called ```Hays_deduper_umicorrect.py```.. My strategy is to include a function that if an UMI is not in the list, there is a function that calculates the UMI in the list with the closest hamming distance to the UMI with errors, and replaced the error UMI with that one for the duplicate analysis. If there are more than 1 UMIs that have the closest hamming distance then the UMI will be considered unknown still.

Ok, I got the function done and made a new unit test called ```input_umi_correct.sam``` and ```output_umi_correct.sam```to test that it works. I included a flag in my argpasrse that takes a boolean to determine if that option will be incorperated or not. I will run it on the actual file with the following command:
```
/usr/bin/time -v ./Hays_deduper_umicorrect.py -u STL96.txt -f <input_file> -o dedup_out.sam -e True
```

I got the expected output for my unit test, which is:
```
Number of alignments kept: 9
Number of duplicates removed: 9
Number of unknown umis discarded: 1
Number of unknown umis corrected: 2
```

Time to run the given test file. I would expect the numbers to be the same as the old script because there were no discarded UMIs in that file. So I should still expect to see 0 discarded lines and 0 corrected UMIs. Here are the results:

```
sbatch run_dedup.sh 
Submitted batch job 23252206

Number of alignments kept: 13719048
Number of duplicates removed: 4467362
Number of unknown umis discarded: 0
Number of unknown umis corrected: 0
```
This is what I expected to see, so I believe my addition to my script works!


### Challenge problems - Choice of duplicate written out

Next, I want to try to implement the choice of duplicate written to the file in a script called ```Hays_deduper_choice.py```. My strategy for this is to output the alignment that has a higher quality, and if the quality is the same I will break the tie with the alignment that is a longer length. If there is still a tie at that point I will arbitrarily pick one. To implement that, I think I need to change my sets to dictionaries where the key is the same as the set, and the value is a list of the quality, length, and sam file line of the alignment. This will be too hard to implement as an option with an argparse so the script is just gonna treat this as default behavior.

Ok, I got the script done and made a new unit test called ```input_choice.sam``` and ```output_choice.sam```to test that it works. I will run it on the actual file with the following command:
```
/usr/bin/time -v ./Hays_deduper_choice.py -u STL96.txt -f <input_file> -o dedup_out.sam
```

I got the expected output for my unit test, which is:
```
Number of alignments kept: 8
Number of duplicates removed: 8
Number of unknown umis discarded: 3
```
And the output file matches.

Time to run the given test file. I would expect the numbers to be the same as the old script because this should not change the number of alignments that are removed. So I should still expect to see 0 discarded lines and 0 corrected UMIs. Here are the results:

```
sbatch run_dedup.sh 
Submitted batch job 23252344

Number of alignments kept: 13719048
Number of duplicates removed: 4467362
Number of unknown umis discarded: 0
```
This is what I expected to see. Additionally, I ran a diff command on the original output vs my new output to see if any lines were changed and there were many lines that were different which tells me a different duplicate is being kept. So I believe my addition to my script works!



### Challenge problems - Option for randomers vs known UMIs

Next, I will tackle randomers vs known UMIs. I am going to make another script called ```Hays_deduper_randomers.py``` that has an option to use 8 nt randomers if a list of known UMIs isn't given. My strategy is to have a check that the UMI is the right length (8 nt) and only has A, G, T, C to be considered a valid UMI if a list of known UMIs is not provided.

Ok, I got the script done and tested it on my ```input_umi_correct.sam```, which should produce the same output as  ```output_umi_correct.sam```to test that it works. I will run it on the actual file with the following command:
```
/usr/bin/time -v ./Hays_deduper_randomers.py -f <input_file> -o dedup_out.sam
```

I got the expected output for my unit test, which is:
```
Number of alignments kept: 9
Number of duplicates removed: 8
Number of unknown umis discarded: 2
```
Time to run the given test file. I would expect the numbers to be the same as the old script because there were no discarded UMIs in that file. So I should still expect to see 0 discarded lines and 0 corrected UMIs. Here are the results:

```
sbatch run_dedup.sh 
Submitted batch job 23252400

Number of alignments kept: 13719048
Number of duplicates removed: 4467362
Number of unknown umis discarded: 0
```
This is what I expected to see, so I believe my addition to my script works!


Should I tackle paired end next???????