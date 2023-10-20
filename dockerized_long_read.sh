#!/usr/bin/env bash
# -*- coding:  utf-8 -*-
# Author: Gaurav Sablok
# date: 2023-10-20
# long read dockerized pipeline easier to execute and also better for the docker application. 
# it will even make your genome browser tracks. This is still in development phase and it will be 
# updated regularly with more options. 

# it supports three long read genome assembler, comprative genome assembly and also the
# update of the previous assembly.

echo                " generating or updating the genome assembly from the long reads such as pacbio"
echo                            "if you have an already assembled genome then use the 
                                update option or else use the genome assembly to update"
read -r -p "do you have all these tools installed such as canu, blasr, lastz, flye, metacatref:": checkavail
if [[ $checkavail == "yes" ]]
then 
   count=0 
   while [[ $count -le 5 ]]
   do
    read -r -p "please provide the path to the canu:": canupath
    read -r -p "please provide the path to the blasr:": blasr
    read -r -p "please provide the path to the lastz:": lastz
    read -r -p "please provide the path to the flye:": flye
    read -r -p "please provide the path to the mecatref:": mecatref
    count+=("count")
    done
        if [[ $count -eq 5 ]]
            then
                break
        fi
else
 conda create -n genomeassembly -y && \
                 conda install -n genomeassembly blasr canu lastz flye mecatref quast -y
 conda clean -t -y
 echo "conda environment has been created"
fi 

canupath="${canupath}"
blasr="${blasr}"
lastz="${lastz}"
flye="${flye}"
mecatref="${mecatref}"
export PATH="${canupath}":$PATH
export PATH="${blasr}":$PATH
export PATH="${lastz}":$PATH
export PATH="${flye}":$PATH
export PATH="${mecatref}":$PATH

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
   read -r -p "please provide the path for the genome assembly reads:" reads
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

echo            "thank you for choosing the workflow options:":
echo        "unless a specified diretory path is provided it will use the $(pwd) as the working directory and all 
                                        analysis will be put int the $(pwd) directory"
closeuptext="You havent selected enough options to run this workflow"
    
if [[ $first == "yes" ]] && 
                    [[ $assembler == "canu" ]] && 
                                [[ $alignment_tracks ]]
then 
    mkdir "$(pwd)/reads_assembly"
    mkdir "$(pwd)/genome_assembly"
    reads_directory="${reads}"
        for i in "${reads_directory}"/*.gz
        do 
            gzip "$i"
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
    cd "$(pwd)/genome_assembly"
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
elif  [[ $first == "yes" ]] && 
                    [[ $assembler == "flye" ]] && 
                              [[ $alignment_tracks ]]
then

    mkdir "$(pwd)/reads_assembly"
    mkdir "$(pwd)/genome_assembly"
    reads_directory="${reads}"
        for i in "${reads_directory}"/*.gz
        do 
            gzip "$i"
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
    cd "$(pwd)/genome_assembly"
        echo "working in the present genome assembly directory"
    unaligned_reads_assembly="$(pwd)/genome_assembly/${species}".unaligned.fasta
    flye --pacbio-raw "${unaligned_reads_assembly}" 
                --genome-size "${calibrate}" \
                --threads "{threads}" \
                --out-dir $(pwd)/genome_assembly 
                --min-overlap "${flyeoverlap}"
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

elif  [[ $first == "yes" ]] && 
                  [[ $assembler == "mecatref" ]] && 
                                  [[ $alignment_tracks ]]
then
    mkdir $(pwd)/reads_assembly
    mkdir $(pwd)/genome_assembly
    reads_directory="${reads}"
        for i in "${reads_directory}"/*.gz
        do 
            gzip "$i"
        done
        cat *.fasta > all_reads.fasta
        mv all_reads.fasta $(pwd)/reads_assembly
    genome_contamant=$read_removal
    if [[ $read_removal ]]
    then
        blasr all_reads.fasta "${genome_contamant}" --best "${mapping}" --bam \
                                  "${species}" --unaligned "${species}".unaligned.fasta
        echo "finished cleaning of the reads"
    fi 
    cp -r "${species}".unaligned.fasta $(pwd)/genome_assembly
    cd $(pwd)/genome_assembly
        echo "working in the present genome assembly directory"
    unaligned_reads_assembly="$(pwd)/genome_assembly/${species}".unaligned.fasta
    reference = "$refgenome"
    mecat2ref -d "${unaligned_reads_assembly}" -r "${reference}" -o "${files%}".reference.sam \
                                            -w "${files%}"_intermediate \
                                                            -t 24 -n 20 -n 50 -m 2 -x 0
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

elif [[ $first == "no" ]] && 
         [[ $update == "yes" ]] &&
                          [[ $assembler == "canu" ]] && 
                                        [[ $alignment_tracks ]]
then
    mkdir $(pwd)/previous_genome_assembly
    mkdir $(pwd)/new_genomic_reads
    mkdir $(pwd)/combined_assembly
    cd $(pwd)/previous_genome_assembly 
    cp -r "${previousassembly}"/* $(pwd)/previous_genome_assembly
    cp -r "${newreadsdirectory}/*" $(pwd)/new_genomic_reads
    cd $(pwd)/new_genomic_reads
    tar zxvf *.gz && cp -r *.fasta $(pwd)/combined_assembly \
                             && cd .. && cd $(pwd)/new_genomic_reads
    for file in *.gz; do gunzip $f; done && \
                               cp -r *.fasta $(pwd)/combined_assembly
    cd $(pwd)/combined_assembly && cat *.fasta >> update_genome_assembly.fasta
    updated_genome_assembly="$(pwd)/combined_assembly/update_genome_assembly.fasta"
    genome_contamant=$read_removal
    if [[ $read_removal ]]
    then
        blasr all_reads.fasta "${genome_contamant}" --best "${mapping}" --bam \
                                  "${updated_genome_assembly}".bam --unaligned \
                                            "${species}".unaligned.fasta
        echo "finished cleaning of the reads"
    fi 
    cp -r "${species}".unaligned.fasta "$(pwd)/genome_assembly"
    mkdir "$(pwd)/refined_genome_assembly"
    canu gridOptions="--time=24:00:00" -corMhapSensitivity="${sensitivity}" 
                                                        -p "{$species}" \ 
                                                        -d "$(pwd)/refined_genome_assembly" \
                                                        -genomeSize="${calibrate}" \
                                                        -minReadLength="${read_selection}" \
                                                        -merylMemory="${merylmemory}" 
                                                        -gnuplotImageFormat=png \
                                                        -ovsThreads="${thread}" \
                                                        -ovbThreads="${thread}" \
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
elif  [[ $first == "no" ]] && 
           [[ $update = "yes" ]] &&
                    [[ $assembler == "flye" ]] && 
                              [[ $alignment_tracks ]]
then

    mkdir "$(pwd)/previous_genome_assembly"
    mkdir "$(pwd)/new_genomic_reads"
    mkdir "$(pwd)/combined_assembly"
    cd "$(pwd)/previous_genome_assembly"
    cp -r "${previousassembly}"/* $(pwd)/previous_genome_assembly
    cp -r "${newreadsdirectory}/*" $(pwd)/new_genomic_reads
    cd "$(pwd)/new_genomic_reads"
    tar zxvf *.gz && cp -r *.fasta $(pwd)/combined_assembly \
                             && cd .. && cd $(pwd)/new_genomic_reads
    for file in *.gz; do gunzip $f; done && \
                               cp -r *.fasta $(pwd)/combined_assembly
    cd "$(pwd)/combined_assembly" && cat *.fasta >> update_genome_assembly.fasta
    updated_genome_assembly="$(pwd)/combined_assembly/update_genome_assembly.fasta"
    genome_contamant=$read_removal
    if [[ $read_removal ]]
    then
        blasr all_reads.fasta "${genome_contamant}" --best "${mapping}" --bam \
                                  "${species}".bam --unaligned "${species}".unaligned.fasta
        echo "finished cleaning of the reads"
    fi 
    cp -r "${species}".unaligned.fasta $(pwd)/genome_assembly
    cd "$(pwd)/genome_assembly"
        echo "working in the present genome assembly directory"
    unaligned_reads_assembly="$(pwd)/genome_assembly/${species}".unaligned.fasta
    flye --pacbio-raw "${unaligned_reads_assembly}" \
                --genome-size "${calibrate}" \
                --threads "{threads}" \
                --out-dir "$(pwd)/genome_assembly" \
                --min-overlap "${flyeoverlap}" 
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

elif  [[ $first == "no" ]] && 
           [[ $update = "yes" ]] &&
                    [[ $assembler == "flye" ]] && 
                              [[ $alignment_tracks ]]
then

    mkdir "$(pwd)/previous_genome_assembly"
    mkdir "$(pwd)/new_genomic_reads"
    mkdir "$(pwd)/combined_assembly"
    cd "$(pwd)/previous_genome_assembly" 
    cp -r "${previousassembly}"/* $(pwd)/previous_genome_assembly
    cp -r "${newreadsdirectory}/*" $(pwd)/new_genomic_reads
    cd "$(pwd)/new_genomic_reads"
    tar zxvf *.gz && cp -r *.fasta $(pwd)/combined_assembly
    cd "$(pwd)/new_genomic_reads"
    for file in *.gz; do gunzip "$f"; done && \
                               cp -r *.fasta "$(pwd)/combined_assembly"
    cd "$(pwd)/combined_assembly" && cat *.fasta >> update_genome_assembly.fasta
    updated_genome_assembly="$(pwd)/combined_assembly/update_genome_assembly.fasta"
    genome_contamant=$read_removal
    if [[ $read_removal ]]
    then
        blasr all_reads.fasta "${genome_contamant}" --best "${mapping}" --bam \
                                  "${species}" --unaligned "${species}".unaligned.fasta
        echo "finished cleaning of the reads"
    fi 
    cp -r "${species}".unaligned.fasta $(pwd)/genome_assembly
    cd "$(pwd)/genome_assembly"
        echo "working in the present genome assembly directory"
    unaligned_reads_assembly="$(pwd)/genome_assembly/${species}".unaligned.fasta
    reference = "$refgenome"
    mecat2ref -d "${unaligned_reads_assembly}" -r "${reference}" -o "${files%}".reference.sam \
                                            -w "${files%}"_intermediate \
                                                            -t 24 -n 20 -n 50 -m 2 -x 0
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
else 
 printf "%s${closeuptext}"   
fi           
