#!/bin/bash

get_input() {
    # Function to parse arguments
    # Specifying usage message
    usage="Usage: pipeline.bash -i <input directory> -o <output directory> -[v]
              Listeria gene prediction pipeline. The options available are:
                        -i : Directory for genome sequences [required]
                        -o : Output directory [required]
                        -h : Print usage instructions"

  # Getopts block, will take in the arguments as inputs and assign them to variables
        while getopts "i:o:h" option; do
                case $option in
                        i) input_directory=$OPTARG;;
                        o) output_directory=$OPTARG;;
                        h) echo "$usage"
                              exit 0;;
                       \?) echo "Invalid option."
                          "$usage"
                              exit 1;;
                esac
        done

  #Check for presence of required arguments
  if [ ! "$input_directory" ] || [ ! "$output_directory" ]
  then
    echo "ERROR: Required arguments missing!"
    echo "$usage"
    exit 1
  fi

  #Check if input directory is a directory

  if [ ! -d $input_directory ]
  then echo "ERROR: Not a valid directory"
  echo "$usage"
  exit 1
  fi

  #Check if output file is already present, give option to rewrite.
	if [ -d $output_directory ]
        then
		echo "Output directory already exists, would you like to overwrite? Reply with y/n"
		read answer
		case $answer in
			y) echo "Overwriting folder $output_directory in subsequent steps";;
			n) echo "Folder overwrite denied, exiting pipeline"
				exit 1;;
			\?) echo "Incorrect option specified, exiting pipeline"
				exit 1;;
		esac
	fi

}

make_temp(){

  #Make temp Directory
  mkdir -p temp

  #Make temp input Directory
  mkdir temp/inputs

  #Copy input directory contents into temp_directory
  cp $input_directory/* temp/inputs/

  #Make output directory
  mkdir -p $output_directory

  #Parse input directory to get list of genomes
  ls temp/inputs/ | xargs -L1 bash -c 'a=${0%.*};ext=${0#*.};echo $a >> temp/genomes_list.txt;echo $ext >>temp/genomes_list.txt'

}

run_ab_initio(){

  #Make directory within temp for prodigal genes
  mkdir temp/prodigal_results

  cat temp/genomes_list.txt | xargs -L2 bash -c 'prodigal -i temp/inputs/"$0"."$1" -f gff -o temp/prodigal_results/"$0".gff'

  #Make directory within temp for genemark genes
  mkdir temp/genemark_results

  #Run genemark on genes
  cat temp/genomes_list.txt | xargs -L2 bash -c 'lib/gms2_linux_64/gms2.pl --seq temp/inputs/"$0"."$1" --genome-type bacteria --output temp/genemark_results/"$0".gff --format gff '

  #Make directory for glimmer results
  mkdir temp/glimmer_results

  #Run glimmer
  cat temp/genomes_list.txt | xargs -L2 bash -c 'lib/glimmer3.02/scripts/g3-from-scratch.csh temp/inputs/"$0"."$1" temp/glimmer_results/"$0"'

  #Convert glimmer output to gff
  cat temp/genomes_list.txt | xargs -L2 bash -c 'lib/glimmer_glimmer_to_gff.py temp/glimmer_results/"$0".predict > temp/glimmer_results/"$0".gff'

}

run_rna (){

  #Make directory within temp for merged results
  mkdir temp/final_rna_results

  #Make directory within temp for trascan genemark_results
  mkdir temp/aragorn

  #Run trnascan on genes
  cat temp/genomes_list.txt | xargs -L2 bash -c 'aragorn -fo -o temp/aragorn/"$0"_aragorn.fna temp/inputs/"$0"."$1"'

  #Make directory within temp for rnammer
  mkdir temp/rnammer

  #Make
  cat temp/genomes_list.txt | xargs -L2 bash -c 'lib/RNAmmer/rnammer -S bac -m lsu,ssu,tsu -multi -gff temp/rnammer/"$0"_rnammer.gff -f temp/rnammer/"$0"_rnammer.fna < temp/inputs/"$0"."$1"'

  cat temp/genomes_list.txt | xargs -L2 bash -c 'cat temp/rnammer/"$0"_rnammer.fna temp/aragorn/"$0"_aragorn.fna > temp/final_rna_results/"$0"_RNA.fna'

  #Move RNA results to final directory
  mv temp/final_rna_results/* $output_directory/
}

merge_predictions(){
  #Make directory within temp for merged results
  mkdir temp/merged_results

  #Find the genes that overlap between the prodigal and genemark predictions (with a minimum fraction of overlap of 0.95)
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/prodigal_results/"$0".gff -b temp/genemark_results/"$0".gff -r -f 0.95 > temp/merged_results/"$0"_prodigal_genemark.gff'

  #Find the genes that overlap between the glimmer and genemark predictions (with a minimum fraction of overlap of 0.95)
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/glimmer_results/"$0".gff -b temp/genemark_results/"$0".gff -r -f 0.95 > temp/merged_results/"$0"_glimmer_genemark.gff'

  #Find the genes that overlap between the prodigal and glimmer predictions (with a minimum fraction of overlap of 0.95)
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/prodigal_results/"$0".gff -b temp/glimmer_results/"$0".gff -r -f 0.95 > temp/merged_results/"$0"_prodigal_glimmer.gff'

  #Find the genes that are predicted by all three tools with a minimum fraction overlap of 0.95 (progigal, genemark, and glimmer)
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/merged_results/"$0"_prodigal_genemark.gff -b temp/glimmer_results/"$0".gff -r -f 0.95 > temp/merged_results/"$0"_prodigal_genemark_glimmer.gff'

  #Find the genes predicted by prodigal and genemark that are not in glimmer
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/merged_results/"$0"_prodigal_genemark.gff -b temp/merged_results/"$0"_prodigal_genemark_glimmer.gff -v  > temp/merged_results/"$0"_prodigal_genemark_only.gff'

  #Find the genes predicted by glimmer and genemark that are not in glimmer
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/merged_results/"$0"_glimmer_genemark.gff -b temp/merged_results/"$0"_prodigal_genemark_glimmer.gff -v  > temp/merged_results/"$0"_glimmer_genemark_only.gff'

  #Find the genes predicted by prodigal and glimmer that are not in genemark
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools intersect -a temp/merged_results/"$0"_prodigal_glimmer.gff -b temp/merged_results/"$0"_prodigal_genemark_glimmer.gff -v  > temp/merged_results/"$0"_prodigal_glimmer_only.gff'

  #Merge all 4 DNA prediction files (prodigal_genemark_only, prodigal_glimmer_only, genemark_glimmer_only, and prodigal_genemark_glimmer)
  cat temp/genomes_list.txt | xargs -L2 bash -c 'cat temp/merged_results/"$0"_prodigal_genemark_glimmer.gff temp/merged_results/"$0"_prodigal_genemark_only.gff temp/merged_results/"$0"_glimmer_genemark_only.gff temp/merged_results/"$0"_prodigal_glimmer_only.gff | sort -n -k 4 > temp/merged_results/"$0"_cds.gff'

  #Get fasta file for the coding regions
  cat temp/genomes_list.txt | xargs -L2 bash -c 'bedtools getfasta -s -fi temp/inputs/"$0"."$1" -bed temp/merged_results/"$0"_cds.gff -fo temp/merged_results/"$0"_cds.fna'
 
  #Get amino acid seqs (first reading frame) from cds file
  cat temp/genomes_list.txt | xargs -L2 bash -c 'transeq temp/merged_results/"$0"_cds.fna -outseq temp/merged_results/"$0"_cds.faa -trim -clean'


  #Move final files to output directory
  mv temp/merged_results/*_cds.faa $output_directory/
  mv temp/merged_results/*_cds.fna $output_directory/
  mv temp/merged_results/*_cds.gff $output_directory/

}

main(){
  get_input "$@"
  make_temp
  run_ab_initio
  run_rna
  merge_predictions
  rm -r temp
}

main "$@"

