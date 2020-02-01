# This is a tutorial for identifying sequences from a genome and adding them to a multiple sequence alingment
You can use this method to add sequences from any genome to add to your multiple sequence alignment (MSA) generated using RAD-Seq, target enrichment and whole genome data. In this tutorial, I will provide the basic steps that can be modified as per your requirements. 

We will use the 4,388 genes dataset from Liu et al. 2017 (https://www.pnas.org/content/114/35/E7282#sec-10) and identify these loci in *Cynopterus brachyotis* genome (Accession ID:SSHV00000000.1). Once we have identified the target regions, we will extract these sequences from *Cynopterus brachyotis* genome and add to our existing MSA for downstream analysis. Do note that this code only works on linux platform. 

## Step 1: BLAST
We will download the *Cynopterus brachyotis* genome from NCBI and make a local database using BLAST. To install BLAST please refer to (https://www.ncbi.nlm.nih.gov/books/NBK279671/))

```
makeblastdb -in cynopterus.fasta -out cynopterus -dbtype nucl -hash_index 
```

Now that the local database is set, we will BLAST the 4,388 loci from Liu et al. 2017 to this database. We will use only the human sequences as query for the BLAST. These human sequences can be obtained from biomart using the list of cDNA names available in Liu et al. 2017.

```
blastn -db cynopterus -query loci_list.fa -out cynopterus_blast_results.out -outfmt 6 -max_target_seqs 1
```
We can open and view the BLAST results in an excel sheet or any other text editor. We can filter the BLAST results based on the E-value criterion and retain loci with single hits only for downstream processing. Once we have identified these loci, we should locate the start and stop positions for further processing and also separate the BLAST results based on sense and antisense strands. We will use these information to generate a bed file to extract DNA sequences from the genome.

## Step 2: Extracting sequences from the genome
We will use bedtools (https://bedtools.readthedocs.io/en/latest/) for this step. We will use the start and stop position information from the previous step and generate a bed file (see the example below)

```
SSHV01000001.1	5411593	5414289
SSHV01000001.1	5399783	5401607
SSHV01000011.1	5382392	5383038
SSHV01000006.1	5378319	5378879
```

We will use the getfasta command within bedtools to extract the sequences that were retained in Step 1 and save these in a file in the fasta format
```
bedtools getfasta -fi cynopterus.fasta -bed cynopterus.bed -fo cynopterus_out.txt
```

## Step 3: Adding desired headers to the fasta file generated in Step 2
We will now change the header in the fasta file following the code provided in the biostars webpage (https://www.biostars.org/p/103089/). To implement this step we will first make a text file with desired header names (see example below). Do note that is code works with fasta file in which the sequence spans a single line only.

```
awk 'NR%2==0' cynopterus_out.txt | paste -d'\n' header_file.txt - > cynopterus_rename.txt
```

Example of the header file with desired names

```
>cynopterus_0001
>cynopterus_0002
>cynopterus_0003
>cynopterus_0004
>cynopterus_0005
>cynopterus_0006
>cynopterus_0007
>cynopterus_0008
>cynopterus_0009
>cynopterus_0010
>cynopterus_0011

```

## Step 4: Split file
We will now split the renamed fasta file into multiple files wherein each file consists of a single sequence. The name of each file will be the same as the fasta header. The code for this was obtained from the following website (https://stackoverflow.com/questions/11818495/split-a-fasta-file-and-rename-on-the-basis-of-first-line).


```
awk '/^>/ {OUT=substr($0,2) ".fa"}; OUT {print >OUT}' cynopterus_rename.txt
```

## Step 5: Rename fasta files based on the human cDNA sequence name
We will now rename the files using move command. Before we run this step, it is better to make a copy of all the files generated in the previous step in a separate folder and work on the files from this new folder. 
Within this new folder we will proceed to rename each file based on its human query sequence. We will use the blast results to get the list of the cDNA names. Ensure that the cDNA names in the list is in the same order as the files in the new folder. Linux will process files in an ascending oder.
In Liu et al. 2017, MSA file of each target locus is named after its human cDNA sequence. Our code will ensure that the file name of the *Cynopterus brachyotis* sequence from each locus in essentially the same as the filname of that locus in Liu et al. 2017. 

```
for file in *; do read line;  mv -v "${file}" "${line}";  done < list.txt
```

Example list file comprising cDNA names

```
ENSG00000105829
ENSG00000164144
ENSG00000125885
ENSG00000089195
ENSG00000089199
ENSG00000125772
ENSG00000183840
ENSG00000143458
ENSG00000143363
ENSG00000143409
```

## Step 6: Add the *Cynopterus brachyotis* sequence to cDNA files
We will use a loop to concatenate the sequences obtained from Step 5 to the 4,388 loci of Liu et al. 2017. For this, keep the renamed fasta files from Step 5 in a separate folder called 'blast_results'. We this excercise, we assume that the cDNA sequences from Liu et al. 2017 are in a folder called 'input_seq' and the results will be saved in a folder called 'output_results'.

Do note that the list file is the same as the one used for the previous step.

```
while read name; 
do cat blast_results/$name input_seq/$name > output_results/output_$name.fa
done < list.txt 
```
Once the loop is complete manually check few alignments to ensure that the correct sequences were added.

## Step 7: Realign the edited files from previous step
We will now realign the sequences within each alignment file using the MAFFT program (https://mafft.cbrc.jp/alignment/software/)

```
for fasta_file in $(ls *.fa)
do
mafft --reorder --adjustdirection --auto $fasta_file > $fasta_file.out
done
```

This will successfully add sequences from any candidate genome (in the present case, it is *Cynopterus brachyotis*) to already existing multiple sequence alignments.
