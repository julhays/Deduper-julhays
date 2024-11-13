# Deduper Challenge

I addressed the following challenges for this assignment:
- Error correction of UMIs
- Choice of duplicate written to file
- Known UMIs vs randomers

File directory:
* ```unit_tests/``` -> input and output unit test encompassing all scenarios
* ```Hays_deduper.py``` -> deduplicating script
* ```STL96.txt``` -> list of known UMIs
* ```lab_notebook.md``` -> notebook for this assignment
* ```part1_pseudocode.md``` -> my intial psuedocode approach for tackling this problem
* ```run_dedup.sh``` -> sbatch file to sort input sam file and run ```Hays_deduper.py```

This challenge branch also includes:
* ```Hays_deduper_umicorrect.py``` -> a new script that incorperates an option for UMI error correction of known UMIs and also an option to use 8 nt randomers instead of UMIs
* ```Hays_deduper_choice.py``` -> a new script version that chooses which duplicate to write out based on what alignment has the highest quality and if the quality is same then the longest length
* ```Hays_deduper_randomers.py``` -> a new script that incorperates an option to use 8 nt randomers instead of known UMIs

