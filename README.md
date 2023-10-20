# dockerized_long_read_genome_assembly
a complete workflow that can be dockerized for the long read assembly, it allows for the genome assembly update as well as it allows for the assembly from the start. if you have the illumina reads it allows for the genome mapping also. long read dockerized pipeline it will easier to execute and also better for the docker application. it will even make your genome browser tracks. 

You can select the case to make your genome browsers tracks and direct plugins into the jbrowse for the visualization or snp calling tracts for the visualization. Finished code: 2023-10-20, next update will be with the genome visualization tracts. it does everything for the pacbio and oxford nanopore sequencing. More support for oxford nanopore in the next release.

```
# block of code is here for viewing
read -r -p "please provide the species name:": species
read -r -p "are you assembling the genome first time:": first
if [[ $first ]]
then 
    read -r -p "please provide the path for the genome assembly reads:" reads
    read -r -p "please provide the genome file for the genome completeness estimation:": genomestimation
    read -r -p "do you have a reference genome or a contaminant file for the read removal:" read_removal
    read -r -p "how many reads you want to keep while mapping:" mapping
fi
read -r -p "are you updating an existing assembly:": update
if [[ $update ]]
then 
    echo "please provide the path for the previous assembly:": previousassembly
    echo "please provide the path for the directory containing the new reads:": newreadsdirectory
fi 
read -r -p "please select the choice of the assembler:": assembler
if [[ $assembler == "canu" ]] 
then 
    read -r -p "please provide the threads you want to use:": threads
    read -r -p "please provide the genome size for the calibration:": calibrate
    read -r -p "please provide the sensitity that you want to use for the genome assembly:": sensitivity
    read -r -p "please provide the minimum length for the read selection:": read_selection
    read -r -p "please provide the memory bindings for the meryl:": merylmemory
elif [[ $assembler == "mecatref" ]]
then 
   read -r -p "please provide the name of the reference genome fasta file:": refgenome
elif [[ $assembler == "flye" ]]
then 
    read -r -p "please provide the minimum overlap for the long reads:": flyeoverlap
fi
read -r -p "do you want to make the genome alignments for the comparison:": alignment_tracks
if [[ $alignment_tracks ]]
then 
    read -r -p "please provide the reference genome for the comparison:": comparison
    read -r -p "please provide the genome gff files for the comparison:": comparisongff
    read -r -p "please provide the threshold sensitivity, if not provided then it will set the
                        default threshold sensitivity of 70%:": threshold
    read -r -p "please provide the coverage for the alignments, if no coverage provided then
                        a default coverage of 60% will be applied:": coverage 
    read -r -p "which format of alignments you want for the modplot sam alignments or the 
                        general format which is a tabular format:": fileformat
fi
read -r -p "do you want the annotations to be done:": annotations
read -r -p "do you want to create the genome tracks:": tracks


if [[ $first == "yes" ]] && 
                    [[ $assembler == "canu" ]] && 
                                [[ $alignment_tracks ]]
then 
    mkdir $(pwd)/reads_assembly
    mkdir $(pwd)/genome_assembly
    reads_directory="${reads}"
        for i in reads_directory
        do 
            gzip $i
        done
        cat *.fasta > all_reads.fasta
        mv all_reads.fasta $(pwd)/reads_assembly
    genome_contamant=$read_removal
    if [[ $read_removal ]]
    then
        blasr all_reads.fasta "${genome_contamant}" --best "${mapping}" --bam \
                                  "${species}".bam --unaligned "${species}".unaligned.fasta
        echo "finished cleaning of the reads"
    fi 
    cp -r "${species}".unaligned.fasta $(pwd)/genome_assembly
    cd $(pwd)/genome_assembly
        echo "working in the present genome assembly directory"
    unaligned_reads_assembly="$(pwd)/genome_assembly/${species}".unaligned.fasta
    canu gridOptions="--time=24:00:00" -corMhapSensitivity="${sensitivity}" 
                                                        -p "{$species}" \ 
                                                        -d "$(pwd)/genome_assembly" \
                                                        -genomeSize="${calibrate}" \
                                                        -minReadLength="${read_selection}" \
                                                        -merylMemory="${merylmemory}" 
                                                        -gnuplotImageFormat=png \
                                                        -ovsThreads="${threads}" \ 
                                                        -ovbThreads="${threads}" \
                                                        -pacbio-raw "${unaligned_reads_assembly}"
    genome_assembly_fasta_file=$(ls *.fasta)
    echo "creating the lastz maps of the assembled genomes to the reference genome"
    if [[ $threshold == "" ]] &&
            [[ $coverage == "" ]] &&
                     [[ $fileformat = "" ]]
    then
        lastz "${compairson}"[multiple] "${genome_assembly_fasta_file}" --identity=75 \
                                                        --coverage=60 --format=sam --ambiguous=iupac \
                                                        --output "${species}".refalignments.sam
    else
    threshold="${threshold}"
    coverage ="${coverage}"
    outputformat="{fileformat}"
       lastz "${compairson}"[multiple] "${genome_assembly_fasta_file}" --identity="${threshold}" \
                                                        --coverage="${coverage}" --format="${outputformat}" \
                                                        --ambiguous=iupac --output "${species}".refalignments.sam
    fi
        echo "genome alignments have been written" 
            echo "perfoming the genome completenes analysis"
            quast "${genome_assembly_fasta_file}" -r ${comparison} -g ${comparisongff} \
                               --fragemented -o "${genome_assembly_fasta_file%}".quast.report \
                                                                 -m 500 -t 40 -e
```



Gaurav Sablok \
ORCID: https://orcid.org/0000-0002-4157-9405 \
WOS: https://www.webofscience.com/wos/author/record/C-5940-2014 \
RubyGems Published: https://rubygems.org/profiles/sablokgaurav \
Python Packages Published : https://pypi.org/user/sablokgaurav/
