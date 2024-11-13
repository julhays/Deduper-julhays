#!/usr/bin/env python

# Author: Jules Hays

# Deduper

import argparse
import re

#set up argpase
def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Remove PCR duplicates from a SAM file of interest")
    parser.add_argument("-f", "--file", help="Designates absolute file path to sorted input sam file", type=str, required=True)
    parser.add_argument("-o", "--outfile", help="designates absolute file path to output the deduplicated sam file", type=str, required = True)
    parser.add_argument("-u", "--umi", help="designates file path containing the list of UMIs", type=str, default='STL96.txt')
    parser.add_argument("-e", "--error_correct_umi", help="designates if incorrect UMIs should be error corrected rather than thrown out", type=bool, default= False)
    return parser.parse_args()

#call get_args to create args object
args = get_args()

#set global variables and assign them to the user inputted values at the function call
sam_file: str = args.file
output: str = args.outfile
umi_file_path: str = args.umi
umi_correct: bool = args.error_correct_umi


#hard coded variables for testing
#sam_file = '/projects/bgmp/jkhay/bioinfo/Bi624/Deduper-julhays/part1-pseudocode/input.sam'
#output = 'dedup_out.sam'
#umi_file_path = 'STL96.txt'


#define necessary functions

#read in umis to make a set of known umis
def make_umi_set(umi_file: str) -> set:
    '''Inputs text file that lists known UMIs with a return in between. Outputs a set that 
    contains each UMI'''
    umi_set = set()   
    with open(umi_file, "rt") as fh:
        for line in fh:
            line = line.strip('\n')
            umi_set.add(line)
    return umi_set


#correct error umis
def correct_umi(wrong_umi: str, umi_set: set) -> str or -1:
    '''Inputs UMI that doesn't match UMIs in the UMI set. Calculates the hamming distance 
    between each known UMI and the errored UMi returns the UMI with the smallest hamming 
    distance. If the smallest hamming distance is for more than 1 UMI, then the corrected 
    UMI cannot be determined and the function will return -1'''

    #set a variables to store information about the lowest hamming distance as you loop through
    lowest_distance = len(wrong_umi) + 1  #set something higher than the highest ppossible hamming distance
    corrected_umi = ""   #stores the lowest hamming distance UMI so far
    duplicate_low_distance = False  #keeps track of if more than 1 UMI had the lowest difference

    #loop through the umi set
    for known_umi in umi_set:

        #calculate the hamming distance between the UMI and the known UMI
        distance = sum(u1 != u2 for u1, u2 in zip(wrong_umi, known_umi))

        #checks if distance is lower than current lowest distance
        if distance < lowest_distance:

            #sets the variables to contain info about the new closest UMI
            lowest_distance = distance
            corrected_umi = known_umi
            duplicate_low_distance = False  #resets this flag when you find a new lowest distance

        #checks if distance is the same as the current distance, because then you wouldn't be able to tell which UMI is right
        elif distance == lowest_distance:
            duplicate_low_distance = True

        #move onto next UMI in the set if the distance isn't <= current lowest distance
        else:
            continue
    
    #if the UMI had multiple matches, return a -1 so my main code known to discard it
    if duplicate_low_distance:
        return -1
    #returns the UMI with the lowest hamming distance
    else:
        return corrected_umi


#calculates the 5' start position
def calc_pos(cigar: str, pos: int, strand: str) -> int:
    '''Takes the cigar string, leftmost start, and strand information from an entry of a SAM file. 
    Adjusts for left side soft clipping for the plus strand, and matches/mismatches, right side soft 
    clipping, deletions, and skipped regions for the minus strand to calculate the 5' start positions'''

    #initalialize dictionary to store all cigar letter values, key = letter, value = sum of numbers associated with letter
    cigar_dict = {'S_left': 0, 'M': 0, 'D': 0, 'I': 0, 'N': 0, 'S': 0, 'H': 0, 'P':0}
    first_entry = True  #hold the value of if its on the first part of the cigar strand
    
    #create a list of tuples with each cigar letter and number pair
    for num, letter in re.findall(r'(\d+)([A-Za-z])?', cigar):
        
        #check if its the leftside S which will appear first
        if letter == 'S' and first_entry:
            cigar_dict['S_left'] = int(num) #add S left value to dictionary
            first_entry = False #permanantly make first entry false after first pass
        
        #case for all other CIGAR string letters
        elif letter in cigar_dict:
            cigar_dict[letter] += int(num) #add to sum of that letter entry
            first_entry = False #permanantly make first entry false after first pass
        
        #else in case a different or weird character is encountered
        else:
            print("Invalid CIGAR string character encountered")

    #calculate adjusted position for plus strand: pos - S_left
    if strand == "plus":
        position = pos - cigar_dict['S_left']
        
    #calculate adjusted position for minus strand: pos + M + S_right + D + N - 1
    elif strand == "minus":
        position = pos + cigar_dict['M'] + cigar_dict['S'] + cigar_dict['D'] + cigar_dict['N'] - 1

    #this shouldn't ever occur but just a check
    else:
        print("not a valid strand")

    return position


#START CODE BODY

#read in known UMIs into a set
umi_set = make_umi_set(umi_file_path)

#set variables to count unknown umi's, outputted lines, and removed duplicates
unknown_umi_count = 0
outputted_line_count = 0
removed_dup_count = 0
corrected_umis = 0

#initialize sets to store outputted SAM line info
plus_info = set() #info will be stored in a tuple of (UMI, posiion) for each entry
minus_info = set()

#sets mutable, can't access by position
#if thing in set: do this -> O(1) lookup


#initialize variable to store current chromosome
curr_chromo = 0

#open sam file to read
with open(sam_file, "rt") as fh:
    #open output write file
    with open(output, "wt") as out:
        #read in each line
        for line in fh:

            #check if its a header line, if so write to output
            if line.startswith("@"):
                out.write(line)

            #all other lines
            else:
                line = line.strip('\n') #strip new line
                info = line.split('\t') #split each column into list entry
                
                #first check to see if a new chromosome is encountered
                chromo = info[2]
                if chromo != curr_chromo: #if is, clear sets and reset curr chromo
                    plus_info = set()
                    minus_info = set()
                    curr_chromo = chromo

                #isolate needed information from entry
                QNAME = info[0]
                FLAG = int(info[1])
                POS = int(info[3])
                CIGAR = info[5]

                #get needed umi information from QNAME, UMI is the last 8 characters
                umi = QNAME[-8:]

                #check if UMI in list of known UMIs, if its not discard
                if umi not in umi_set:
                    #check if error correct UMI option is on
                    if umi_correct:
                        #calculate correct UMI
                        new_umi = correct_umi(umi, umi_set)

                        #check if a valid umi was found
                        if new_umi == -1:
                            unknown_umi_count += 1
                            continue
                        else:
                            umi = new_umi
                            corrected_umis += 1

                    #if UMI correct option not on, just discard the UMI if not in set
                    else:
                        unknown_umi_count += 1
                        continue

                #determine what strand the query is on
                if ((FLAG & 16) == 16):  #minus strand
                    #calculate 5' start position
                    position = calc_pos(CIGAR, POS, "minus")

                    #create tuple of duplication information
                    dup_info = (umi, position)

                    #check to see if tuple already in minus set
                    if dup_info in minus_info:
                        #if it is, move on because a dup has already been written out
                        removed_dup_count += 1
                        continue
                    else: #write to output because a dup hasn't been encountered yet
                        minus_info.add(dup_info)
                        output = f'{line}\n'
                        out.write(output)
                        outputted_line_count += 1

                else: #plus strand
                    #calculate 5' start position
                    position = calc_pos(CIGAR, POS, "plus")

                    #create tuple of duplication information
                    dup_info = (umi, position)

                    #check to see if tuple already in plus set
                    if dup_info in plus_info:
                        #if it is, move on because a dup has already been written out
                        removed_dup_count += 1
                        continue
                    else: #write to output because a dup hasn't been encountered yet
                        plus_info.add(dup_info)
                        output = f'{line}\n'
                        out.write(output)
                        outputted_line_count += 1


print(f'Number of alignments kept: {outputted_line_count}')
print(f'Number of duplicates removed: {removed_dup_count}')
print(f'Number of unknown umis discarded: {unknown_umi_count}')

if umi_correct:
    print(f'Number of unknown umis corrected: {corrected_umis}')
