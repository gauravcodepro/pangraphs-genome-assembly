#!/usr/bin/env bash
# -*- coding:  utf-8 -*-
# Universität Potsdam, Germany
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
echo    "it requires canu, bamtools, bowtie2, blasr, lastz, flye, metacatref,pilon, jasper and busco for the processing"
echo    "if you havent downloaded or installed then please press the checkavail as yes"
read -r -p "do you have all these tools installed such as :": checkavail
if [[ $checkavail == "yes" ]]
then
    read -r -p "please provide the path to the canu:": canupath
    read -r -p "please provide the path to the blasr:": blasr
    read -r -p "please provide the path to the lastz:": lastz
    read -r -p "please provide the path to the flye:": flye
    read -r -p "please provide the path to the mecatref:": mecatref
    read -r -p "please provide the path to the bamtools:": bamtools
    read -r -p "please provide the path to the bowtie2:": bowtie2
    read -r -p "please provide the path to the pilon:": pilon
    declare -a path=([1]="${canupath}" [2]="${blasr}" [3]="${lastz}" 
                           [4]="${flye}" [5]="${mecatref}" [6]="${bamtools}" 
                                                                     [7]="${bowtie2}" [8]="${pilon}")
    for  ((i=1; i<=8; i++))
    do
	echo the provided paths are paths: "${path[i]}"
    done
else
    conda create -n genomeassembly -y && \
                 conda install -n genomeassembly blasr canu lastz flye mecatref quast pilon bowtie2 -y
    conda clean -t -y
 echo "conda environment has been created"
fi 

canupath="${canupath}"
blasr="${blasr}"
lastz="${lastz}"
flye="${flye}"
mecatref="${mecatref}"
bamtools="${bamtools}"
bowtie2="${bowtie2}"
pilon="${pilon}"
export PATH="${canupath}":$PATH
export PATH="${blasr}":$PATH
export PATH="${lastz}":$PATH
export PATH="${flye}":$PATH
export PATH="${mecatref}":$PATH
export PATH="${bamtools}":$PATH
export PATH="${bowtie2}":$PATH
export PATH="${pilon}":$PATH

read -r -p "please provide the species name:": species
read -r -p "are you assembling the genome first time:": first
if [[ $first == "yes" ]]
then 
    read -r -p "please provide the path for the genome assembly reads:" reads
    read -r -p "do you have a reference genome or a contaminant file for the read removal:" removal
    if [[ $removal == "yes" ]]
    then 
        read -r -p "please provide the reference genome or a contaminant file:": read_removal
    fi
read -r -p "how many reads you want to keep while blasr mapping:" mapping
read -r -p "do you want to calculate the coverage also:": coverage
read -r -p "do you want to polish the genome also:": polish
if [[ $polish == "yes" ]]
then 
    read -r -p "which polishing way you want to use pilon or the jasper":polishoption
fi
fi
read -r -p "are you updating an existing assembly:": update
if [[ $update == "yes" ]]
then 
    read -r -p "please provide the path for the previous assembly:": previousassembly
    read -r -p "please provide the path for the directory containing the new reads:": newreadsdirectory
    read -r -p "do you want to calculate the coverage also:": coverage
    read -r -p "do you want to polish the genome also:": polish
    read -r -p "how many reads you want to keep while mapping:" mapping
fi
if [[ $coverage == "yes" ]]
then 
    read -r -p "please provide the directory containing the illumina reads:" illuminareads
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
if [[ $alignment_tracks == "yes" ]]
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

echo            "thank you for choosing the workflow options:":
echo        "unless a specified diretory path is provided it will use the $(pwd) as the working directory and all 
                                        analysis will be put int the $(pwd) directory"
closeuptext="You havent selected enough options to run this workflow"
    
if  [[ $species ]]  && [[ $first == "yes" ]] && 
                 [[ $assembler == "canu" ]] && [[ $alignment_tracks == "yes" ]] &&
                 [[ $mapping == "yes" ]] && [[ $coverage == "yes" ]] && [[ $polish == "yes" ]] && [[ $polishoption ]]
then 
    mkdir "$(pwd)/reads_assembly"
    mkdir "$(pwd)/genome_assembly"
    mkdir "$(pwd)/illumina_reads"
    mkdir "$(pwd)/polished_genome"
    mkdir "$(pwd)/sam_file_illumina"
    longreads="$(pwd)/reads_assembly"
    genome_assembly="$(pwd)/genome_assembly"
    illuminareads="$(pwd)/illumina_reads"
    polishedgenome="$(pwd)/polished_genome"
    coverage=$(pwd)/sam_file_illumina
    reads_directory="${reads}"
        for i in "${reads_directory}"/*.gz
        do 
            gzip "$i"
        done
        cat *.fasta > all_reads.fasta
        mv all_reads.fasta "${longreads}"
    illumina_directory="${illuminareads}"
    for file in "${illumina_directory}"/*.gz
    do 
        gunzip "$file"
    done
    for file in "${illumina_directory}"/*.R1.fastq
    do 
        cat *.R1.fastq >> final.R1.fastq
        echo "$file" > illumina.R1.txt
    done
        echo final.R1.fastq >> illumina.R1.final.txt
    for file in "${illumina_directory}"/*.R2.fastq
    do 
        cat *.R2.fastq >> final.R2.fastq
        echo "$file" > illumina.R2.txt
    done
        echo final.R2.fastq >> illumina.R1.final.txt
    cd "${longreads}"
        blasr all_reads.fasta "${read_removal}" --best "${mapping}" --bam \
                                      "${species}".bam --unaligned "${species}".unaligned.fasta
        echo "finished cleaning of the reads"
    cp -r "${species}".unaligned.fasta "${genome_assembly}"
    (cd ..)

    cd "${genome_assembly}"
        echo "working in the present genome assembly directory"
    canu gridOptions="--time=24:00:00" -corMhapSensitivity="${sensitivity}" \
                                                        -p "{$species}" \
                                                        -d "${genomeassembly}" \
                                                        -genomeSize="${calibrate}" \
                                                        -minReadLength="${read_selection}" \
                                                        -merylMemory="${merylmemory}" 
                                                        -gnuplotImageFormat=png \
                                                        -ovsThreads="${threads}" \
                                                        -ovbThreads="${threads}" \
                                                        -pacbio-raw "${species}".unaligned.fasta
    genome_assembly_fasta_file=$(pwd)/$(ls *.fasta)
    echo "creating the lastz maps of the assembled genomes to the reference genome"
    if [[ $threshold == "" ]] &&
            [[ $coverage == "" ]] &&
                     [[ $fileformat = "" ]]
    then
        lastz "${comparison}"[multiple] "${genome_assembly_fasta_file}" --identity=75 \
                                                        --coverage=60 --format=sam --ambiguous=iupac \
                                                                    --output "${species}".refalignments.sam
    else
    threshold="${threshold}"
    coverage ="${coverage}"
    outputformat="{fileformat}"
       lastz "${comparison}"[multiple] "${genome_assembly_fasta_file}" --identity="${threshold}" \
                                                        --coverage="${coverage}" --format="${outputformat}" \
                                                        --ambiguous=iupac --output "${species}".refalignments.sam
    fi
        echo "genome alignments have been written" 
            echo "perfoming the genome completenes analysis"
            quast "${genome_assembly_fasta_file}" -r "${comparison}" -g "${comparisongff}" \
                               --fragemented -o "${genome_assembly_fasta_file%}".quast.report \
                                                                                       -m 500 -t 40 -e
    echo "starting the genome coverage analysis"
    cd "${illumina_directory}"
    cp "${genome_assembly_fasta_file}" "${illumina_directory}"
    paste illumina.R1.txt illumina.R2.txt | while read -r col1 col2; \
    do 
        bowtie2-build "${genome_assembly_fasta_file}" "${species}"
        bowtie2 -t -x "${species}" -p "${threads}" --very-sensitive-local -1 "${col1}" \
                                 -2 "${col2}" -S "${species}".sam --no-unal --al-conc "${species}".aligned.fastq 
    done
    cp -r "${species}".sam "${polishedgenome}"
    (cd ..) 
    cd "${coverage}"
    for f in *.sam; do samtools view -bam "${f}".bam -sam "${f}"; done
    bamtools coverage -in "${f}".bam -out "${species}".bam.coverage.txt
    bamtools stats -in "${f}".bam -insert > "${species}".bam.stats.coverage.illumina.txt
    bamtools sort -in "${f}".bam -out "${f}".sorted.bam -n 500000 -mem 1024 
    cp -r "${f}".sorted.bam "${polishedgenome}"
    (cd ..)

    cd "${polishedgenome}"
    cp -r "${genome_assembly_fasta_file}" "${polishedgenome}"
    if [[ $polishoption == "jasper" ]]
    then 
    wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz 
        if [[ ! -f "jasper-1.0.3.tar.gz" ]]
        then
            echo "jasper file didnt downloaded and retrying"
            wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz
        elif [[ -f "jasper-1.0.3.tar.gz" ]]
        then 
            echo "file downloaded and processing for the installation"
        fi 
        tar zxvf jasper-1.0.3.tar.gz
        cd jasper-1.0.3
        ./configure --prefix=$PWD
        make install
        export PYTHONPATH=$(pwd)
        (cd ..)
    export PATH="${polishedgenome}"/jasper-1.0.3:$PATH
    jasper.sh -t "${threads}" "${genome_assembly_fasta_file}" \
            -r "${illumina_directory}"/final.R1.fastq "${illumina_directory}"/final.R1.fastq
    else 
    java -Xmx100g -jar pilon-1.28.jar --genome "${genome_assembly_fasta_file}" \
                    --frags "${f}".sorted.bam --output "${species}".polished_genome \
                                --threads 60 --changes --fix all --mindepth 10
    fi 
    (cd ..)
    echo "processing finished"
 
elif  [[ $species ]]  && [[ $first == "yes" ]] && 
                 [[ $assembler == "flye" ]] && [[ $alignment_tracks == "yes" ]] &&
                 [[ $mapping == "yes" ]] && [[ $coverage == "yes" ]] && [[ $polish == "yes" ]] && [[ $polishoption ]]

then 
    mkdir "$(pwd)/reads_assembly"
    mkdir "$(pwd)/genome_assembly"
    mkdir "$(pwd)/illumina_reads"
    mkdir "$(pwd)/polished_genome"
    mkdir "$(pwd)/sam_file_illumina"
    longreads="$(pwd)/reads_assembly"
    genome_assembly="$(pwd)/genome_assembly"
    illuminareads="$(pwd)/illumina_reads"
    polishedgenome="$(pwd)/polished_genome"
    coverage=$(pwd)/sam_file_illumina
    reads_directory="${reads}"
        for i in "${reads_directory}"/*.gz
        do 
            gzip "$i"
        done
        cat *.fasta > all_reads.fasta
        mv all_reads.fasta "${longreads}"
    illumina_directory="${illuminareads}"
    for file in "${illumina_directory}"/*.gz
    do 
        gunzip "$file"
    done
    for file in "${illumina_directory}"/*.R1.fastq
    do 
        cat *.R1.fastq >> final.R1.fastq
        echo "$file" > illumina.R1.txt
    done
        echo final.R1.fastq >> illumina.R1.final.txt
    for file in "${illumina_directory}"/*.R2.fastq
    do 
        cat *.R2.fastq >> final.R2.fastq
        echo "$file" > illumina.R2.txt
    done
        echo final.R2.fastq >> illumina.R1.final.txt
    cd "${longreads}"

        blasr all_reads.fasta "${read_removal}" --best "${mapping}" --bam \
                                      "${species}".bam --unaligned "${species}".unaligned.fasta
        echo "finished cleaning of the reads"
    cp -r "${species}".unaligned.fasta "${genome_assembly}"
    (cd ..)

    cd "${genome_assembly}"
        echo "working in the present genome assembly directory"
    flye --pacbio-raw "${species}".unaligned.fasta \
                --genome-size "${calibrate}" \
                --threads "{threads}" \
                --out-dir $(pwd)/genome_assembly \
                --min-overlap "${flyeoverlap}" \
    genome_assembly_fasta_file="$(ls *.fasta)"
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
            quast "${genome_assembly_fasta_file}" -r "${comparison}" -g "${comparisongff}" \
                               --fragemented -o "${genome_assembly_fasta_file%}".quast.report \
                                                                 -m 500 -t 40 -e 
  echo "starting the genome coverage analysis"
    cd "${illumina_directory}"
    cp "${genome_assembly_fasta_file}" "${illumina_directory}"
    paste illumina.R1.txt illumina.R2.txt | while read -r col1 col2; \
    do 
        bowtie2-build "${genome_assembly_fasta_file}" "${species}"
        bowtie2 -t -x "${species}" -p "${threads}" --very-sensitive-local -1 "${col1}" \
                                 -2 "${col2}" -S "${species}".sam --no-unal --al-conc "${species}".aligned.fastq 
    done
    cp -r "${species}".sam "${polishedgenome}"
    (cd ..) 
    cd "${coverage}"
    for f in *.sam; do samtools view -bam "${f}".bam -sam "${f}"; done
    bamtools coverage -in "${f}".bam -out "${species}".bam.coverage.txt
    bamtools stats -in "${f}".bam -insert > "${species}".bam.stats.coverage.illumina.txt
    bamtools sort -in "${f}".bam -out "${f}".sorted.bam -n 500000 -mem 1024 
    cp -r "${f}".sorted.bam "${polishedgenome}"
    (cd ..)

    cd "${polishedgenome}"
    cp -r "${genome_assembly_fasta_file}" "${polishedgenome}"
    if [[ $polishoption == "jasper" ]]
    then 
    wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz 
        if [[ ! -f "jasper-1.0.3.tar.gz" ]]
        then
            echo "jasper file didnt downloaded and retrying"
            wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz
        elif [[ -f "jasper-1.0.3.tar.gz" ]]
        then 
            echo "file downloaded and processing for the installation"
        fi 
        tar zxvf jasper-1.0.3.tar.gz
        cd jasper-1.0.3
        ./configure --prefix=$PWD
        make install
        export PYTHONPATH=$(pwd)
        (cd ..)
    export PATH="${polishedgenome}"/jasper-1.0.3:$PATH
    jasper.sh -t "${threads}" "${genome_assembly_fasta_file}" \
            -r "${illumina_directory}"/final.R1.fastq "${illumina_directory}"/final.R1.fastq
    else 
    java -Xmx100g -jar pilon-1.28.jar --genome "${genome_assembly_fasta_file}" \
                    --frags "${f}".sorted.bam --output "${species}".polished_genome \
                                --threads 60 --changes --fix all --mindepth 10
    fi 
    (cd ..)
    echo "processing finished"

elif  [[ $species ]]  && [[ $first == "yes" ]] && 
                 [[ $assembler == "mecatref" ]] && [[ $alignment_tracks == "yes" ]] &&
                 [[ $mapping == "yes" ]] && [[ $coverage == "yes" ]] && [[ $polish == "yes" ]] && [[ $polishoption ]]
then 
    mkdir "$(pwd)/reads_assembly"
    mkdir "$(pwd)/genome_assembly"
    mkdir "$(pwd)/illumina_reads"
    mkdir "$(pwd)/polished_genome"
    mkdir "$(pwd)/sam_file_illumina"
    longreads="$(pwd)/reads_assembly"
    genome_assembly="$(pwd)/genome_assembly"
    illuminareads="$(pwd)/illumina_reads"
    polishedgenome="$(pwd)/polished_genome"
    coverage=$(pwd)/sam_file_illumina
    reads_directory="${reads}"
        for i in "${reads_directory}"/*.gz
        do 
            gzip "$i"
        done
        cat *.fasta > all_reads.fasta
        mv all_reads.fasta "${longreads}"
    illumina_directory="${illuminareads}"
    for file in "${illumina_directory}"/*.gz
    do 
        gunzip "$file"
    done
   for file in "${illumina_directory}"/*.R1.fastq
    do 
        cat *.R1.fastq >> final.R1.fastq
        echo "$file" > illumina.R1.txt
    done
        echo final.R1.fastq >> illumina.R1.final.txt
    for file in "${illumina_directory}"/*.R2.fastq
    do 
        cat *.R2.fastq >> final.R2.fastq
        echo "$file" > illumina.R2.txt
    done
        echo final.R2.fastq >> illumina.R1.final.txt
    cd "${longreads}"
    genome_contamant=${read_removal}
        blasr all_reads.fasta "${genome_contamant}" --best "${mapping}" --bam \
                                      "${species}".bam --unaligned "${species}".unaligned.fasta
        echo "finished cleaning of the reads"
    cp -r "${species}".unaligned.fasta "${genome_assembly}"
    (cd ..)

    cd "${genome_assembly}"
        echo "working in the present genome assembly directory"
    reference="$refgenome"
    mecat2ref -d "${species}".unaligned.fasta -r "${reference}" -o "${files%}".reference.sam \
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
    echo "starting the genome coverage analysis"
    cd "${illumina_directory}"
    cp "${genome_assembly_fasta_file}" "${illumina_directory}"
    paste illumina.R1.txt illumina.R2.txt | while read -r col1 col2; \
    do 
        bowtie2-build "${genome_assembly_fasta_file}" "${species}"
        bowtie2 -t -x "${species}" -p "${threads}" --very-sensitive-local -1 "${col1}" \
                                 -2 "${col2}" -S "${species}".sam --no-unal --al-conc "${species}".aligned.fastq 
    done
    cp -r "${species}".sam "${polishedgenome}"
    (cd ..) 
    cd "${coverage}"
    for f in *.sam; do samtools view -bam "${f}".bam -sam "${f}"; done
    bamtools coverage -in "${f}".bam -out "${species}".bam.coverage.txt
    bamtools stats -in "${f}".bam -insert > "${species}".bam.stats.coverage.illumina.txt
    bamtools sort -in "${f}".bam -out "${f}".sorted.bam -n 500000 -mem 1024 
    cp -r "${f}".sorted.bam "${polishedgenome}"
    (cd ..)

    cd "${polishedgenome}"
    cp -r "${genome_assembly_fasta_file}" "${polishedgenome}"
    if [[ $polishoption == "jasper" ]]
    then 
    wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz 
        if [[ ! -f "jasper-1.0.3.tar.gz" ]]
        then
            echo "jasper file didnt downloaded and retrying"
            wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz
        elif [[ -f "jasper-1.0.3.tar.gz" ]]
        then 
            echo "file downloaded and processing for the installation"
        fi 
        tar zxvf jasper-1.0.3.tar.gz
        cd jasper-1.0.3
        ./configure --prefix=$PWD
        make install
        export PYTHONPATH=$(pwd)
        (cd ..)
    export PATH="${polishedgenome}"/jasper-1.0.3:$PATH
    jasper.sh -t "${threads}" "${genome_assembly_fasta_file}" \
            -r "${illumina_directory}"/final.R1.fastq "${illumina_directory}"/final.R1.fastq
    else 
    java -Xmx100g -jar pilon-1.28.jar --genome "${genome_assembly_fasta_file}" \
                    --frags "${f}".sorted.bam --output "${species}".polished_genome \
                                --threads 60 --changes --fix all --mindepth 10
    fi 
    (cd ..)
    echo "processing finished"

elif [[ $species ]] && [[ $species ]]  && [[ $first == "no" ]] && [[ $update == "yes" ]]
                 [[ $assembler == "mecatref" ]] && [[ $alignment_tracks == "yes" ]] &&
                 [[ $mapping == "yes" ]] && [[ $coverage == "yes" ]] && [[ $polish == "yes" ]] && [[ $polishoption ]]
then
    mkdir "$(pwd)/new_genomic_reads"
    mkdir "$(pwd)/combined_assembly"
    mkdir "$(pwd)/illumina_reads"
    mkdir "$(pwd)/polished_genome"
    mkdir "$(pwd)/sam_file_illumina"
    newlongreads="$(pwd)/new_genomic_reads"
    combinedassembly="$(pwd)/combined_assembly"
    illuminareads="$(pwd)/illumina_reads"
    polishedgenome="$(pwd)/polished_genome"
    coverage=$(pwd)/sam_file_illumina
    reads_directory="${illuminareads}"
        for i in "${reads_directory}"/*.gz
        do 
            gzip "$i"
        done
        cat *.fasta > all_reads.fasta
        mv all_reads.fasta "${longreads}"
    illumina_directory="${illuminareads}"
    for file in "${illumina_directory}"/*.gz
    do 
        gunzip "$file"
    done
    for file in "${illumina_directory}"/*.R1.fastq
    do 
        cat *.R1.fastq >> final.R1.fastq
        echo "$file" > illumina.R1.txt
    done
        echo final.R1.fastq >> illumina.R1.final.txt
    for file in "${illumina_directory}"/*.R2.fastq
    do 
        cat *.R2.fastq >> final.R2.fastq
        echo "$file" > illumina.R2.txt
    done
        echo final.R2.fastq >> illumina.R1.final.txt
    cp -r "${previousassembly}"/*.fasta "${combinedassembly}"
    cp -r "${newreadsdirectory}/*.gz" "${newlongreads}"
    cd "${newlongreads}"
    tar zxvf *.gz && cp -r *.fasta "${combinedassembly}"
    (cd ..)
    cd "${combinedassembly}"
       cat *.fasta >> update_genome_assembly.fasta         
    updated_genome_assembly=""${combinedassembly}"/update_genome_assembly.fasta"
        blasr all_reads.fasta "${read_removal}" --best "${mapping}" --bam \
                                        "${updated_genome_assembly}".bam --unaligned \
                                                                 "${species}".unaligned.fasta
        canu gridOptions="--time=24:00:00" -corMhapSensitivity="${sensitivity}" \
                                                        -p "{$species}" \
                                                        -d "${combinedassembly}" \
                                                        -genomeSize="${calibrate}" \
                                                        -minReadLength="${read_selection}" \
                                                        -merylMemory="${merylmemory}" 
                                                        -gnuplotImageFormat=png \
                                                        -ovsThreads="${thread}" \
                                                        -ovbThreads="${thread}" \
                                                        -pacbio-raw "${unaligned_reads_assembly}"
    genome_assembly_fasta_file=${(pwd)/$(ls *.fasta)}
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
            quast "${genome_assembly_fasta_file}" -r "${comparison}" -g "${comparisongff}" \
                               --fragemented -o "${genome_assembly_fasta_file%}".quast.report \
                                                                -m 500 -t 40 -e
    echo "starting the genome coverage analysis"
    cd "${illumina_directory}"
    cp "${genome_assembly_fasta_file}" "${illumina_directory}"
    paste illumina.R1.txt illumina.R2.txt | while read -r col1 col2; \
    do 
        bowtie2-build "${genome_assembly_fasta_file}" "${species}"
        bowtie2 -t -x "${species}" -p "${threads}" --very-sensitive-local -1 "${col1}" \
                                 -2 "${col2}" -S "${species}".sam --no-unal --al-conc "${species}".aligned.fastq 
    done
    cp -r "${species}".sam "${polishedgenome}"
    (cd ..) 
    cd "${coverage}"
    for f in *.sam; do samtools view -bam "${f}".bam -sam "${f}"; done
    bamtools coverage -in "${f}".bam -out "${species}".bam.coverage.txt
    bamtools stats -in "${f}".bam -insert > "${species}".bam.stats.coverage.illumina.txt
    bamtools sort -in "${f}".bam -out "${f}".sorted.bam -n 500000 -mem 1024 
    cp -r "${f}".sorted.bam "${polishedgenome}"
    (cd ..)

    cd "${polishedgenome}"
    cp -r "${genome_assembly_fasta_file}" "${polishedgenome}"
    if [[ $polishoption == "jasper" ]]
    then 
    wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz 
        if [[ ! -f "jasper-1.0.3.tar.gz" ]]
        then
            echo "jasper file didnt downloaded and retrying"
            wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz
        elif [[ -f "jasper-1.0.3.tar.gz" ]]
        then 
            echo "file downloaded and processing for the installation"
        fi 
        tar zxvf jasper-1.0.3.tar.gz
        cd jasper-1.0.3
        ./configure --prefix=$PWD
        make install
        export PYTHONPATH=$(pwd)
        (cd ..)
    export PATH="${polishedgenome}"/jasper-1.0.3:$PATH
    jasper.sh -t "${threads}" "${genome_assembly_fasta_file}" \
            -r "${illumina_directory}"/final.R1.fastq "${illumina_directory}"/final.R1.fastq
    else 
    java -Xmx100g -jar pilon-1.28.jar --genome "${genome_assembly_fasta_file}" \
                    --frags "${f}".sorted.bam --output "${species}".polished_genome \
                                --threads 60 --changes --fix all --mindepth 10
    fi 
    (cd ..)
    echo "processing finished"
    echo "processing finished"
    
elif  [[ $species ]] && [[ $species ]]  && [[ $first == "no" ]] && [[ $update == "yes" ]]
                 [[ $assembler == "flye" ]] && [[ $alignment_tracks == "yes" ]] &&
                 [[ $mapping == "yes" ]] && [[ $coverage == "yes" ]] && [[ $polish == "yes" ]] && [[ $polishoption ]]
then
    mkdir "$(pwd)/new_genomic_reads"
    mkdir "$(pwd)/combined_assembly"
    mkdir "$(pwd)/illumina_reads"
    mkdir "$(pwd)/polished_genome"
    mkdir "$(pwd)/sam_file_illumina"
    newlongreads="$(pwd)/new_genomic_reads"
    combinedassembly="$(pwd)/combined_assembly"
    illuminareads="$(pwd)/illumina_reads"
    polishedgenome="$(pwd)/polished_genome"
    coverage=$(pwd)/sam_file_illumina
    reads_directory="${illuminareads}"
        for i in "${reads_directory}"/*.gz
        do 
            gzip "$i"
        done
        cat *.fasta > all_reads.fasta
        mv all_reads.fasta "${longreads}"
    illumina_directory="${illuminareads}"
    for file in "${illumina_directory}"/*.gz
    do 
        gunzip "$file"
    done
   for file in "${illumina_directory}"/*.R1.fastq
    do 
        cat *.R1.fastq >> final.R1.fastq
        echo "$file" > illumina.R1.txt
    done
        echo final.R1.fastq >> illumina.R1.final.txt
    for file in "${illumina_directory}"/*.R2.fastq
    do 
        cat *.R2.fastq >> final.R2.fastq
        echo "$file" > illumina.R2.txt
    done
        echo final.R2.fastq >> illumina.R1.final.txt
    cp -r "${previousassembly}"/*.fasta "${combinedassembly}"
    cp -r "${newreadsdirectory}/*.gz" "${newlongreads}"
    cd "${newlongreads}"
    tar zxvf *.gz && cp -r *.fasta "${combinedassembly}"
    (cd ..)
    cd "${combinedassembly}"
       cat *.fasta >> update_genome_assembly.fasta         
    updated_genome_assembly=""${combinedassembly}"/update_genome_assembly.fasta"
        blasr all_reads.fasta "${read_removal}" --best "${mapping}" --bam \
                                        "${updated_genome_assembly}".bam --unaligned \
                                                                 "${species}".unaligned.fasta
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
            quast "${genome_assembly_fasta_file}" -r "${comparison}" -g "${comparisongff}" \
                               --fragemented -o "${genome_assembly_fasta_file%}".quast.report \
                                                                 -m 500 -t 40 -e 

    echo "starting the genome coverage analysis"
    cd "${illumina_directory}"
    cp "${genome_assembly_fasta_file}" "${illumina_directory}"
    paste illumina.R1.txt illumina.R2.txt | while read -r col1 col2; \
    do 
        bowtie2-build "${genome_assembly_fasta_file}" "${species}"
        bowtie2 -t -x "${species}" -p "${threads}" --very-sensitive-local -1 "${col1}" \
                                 -2 "${col2}" -S "${species}".sam --no-unal --al-conc "${species}".aligned.fastq 
    done
    cp -r "${species}".sam "${polishedgenome}"
    (cd ..) 
    cd "${coverage}"
    for f in *.sam; do samtools view -bam "${f}".bam -sam "${f}"; done
    bamtools coverage -in "${f}".bam -out "${species}".bam.coverage.txt
    bamtools stats -in "${f}".bam -insert > "${species}".bam.stats.coverage.illumina.txt
    bamtools sort -in "${f}".bam -out "${f}".sorted.bam -n 500000 -mem 1024 
    cp -r "${f}".sorted.bam "${polishedgenome}"
    (cd ..)

    cd "${polishedgenome}"
    cp -r "${genome_assembly_fasta_file}" "${polishedgenome}"
    if [[ $polishoption == "jasper" ]]
    then 
    wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz 
        if [[ ! -f "jasper-1.0.3.tar.gz" ]]
        then
            echo "jasper file didnt downloaded and retrying"
            wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz
        elif [[ -f "jasper-1.0.3.tar.gz" ]]
        then 
            echo "file downloaded and processing for the installation"
        fi 
        tar zxvf jasper-1.0.3.tar.gz
        cd jasper-1.0.3
        ./configure --prefix=$PWD
        make install
        export PYTHONPATH=$(pwd)
        (cd ..)
    export PATH="${polishedgenome}"/jasper-1.0.3:$PATH
    jasper.sh -t "${threads}" "${genome_assembly_fasta_file}" \
            -r "${illumina_directory}"/final.R1.fastq "${illumina_directory}"/final.R1.fastq
    else 
    java -Xmx100g -jar pilon-1.28.jar --genome "${genome_assembly_fasta_file}" \
                    --frags "${f}".sorted.bam --output "${species}".polished_genome \
                                --threads 60 --changes --fix all --mindepth 10
    fi 
    (cd ..)
    echo "processing finished"

elif  [[ $species ]] && [[ $species ]]  && [[ $first == "no" ]] && [[ $update == "yes" ]]
                 [[ $assembler == "mecatref" ]] && [[ $alignment_tracks == "yes" ]] &&
                 [[ $mapping == "yes" ]] && [[ $coverage == "yes" ]] && [[ $polish == "yes" ]] && [[ $polishoption ]]
then

    mkdir "$(pwd)/new_genomic_reads"
    mkdir "$(pwd)/combined_assembly"
    mkdir "$(pwd)/illumina_reads"
    mkdir "$(pwd)/polished_genome"
    mkdir "$(pwd)/sam_file_illumina"
    newlongreads="$(pwd)/new_genomic_reads"
    combinedassembly="$(pwd)/combined_assembly"
    illuminareads="$(pwd)/illumina_reads"
    polishedgenome="$(pwd)/polished_genome"
    coverage=$(pwd)/sam_file_illumina
    reads_directory="${illuminareads}"
        for i in "${reads_directory}"/*.gz
        do 
            gzip "$i"
        done
        cat *.fasta > all_reads.fasta
        mv all_reads.fasta "${longreads}"
    illumina_directory="${illuminareads}"
    for file in "${illumina_directory}"/*.gz
    do 
        gunzip "$file"
    done
    for file in "${illumina_directory}"/*.R1.fastq
    do 
        cat *.R1.fastq >> final.R1.fastq
        echo "$file" > illumina.R1.txt
    done
        echo final.R1.fastq >> illumina.R1.final.txt
    for file in "${illumina_directory}"/*.R2.fastq
    do 
        cat *.R2.fastq >> final.R2.fastq
        echo "$file" > illumina.R2.txt
    done
        echo final.R2.fastq >> illumina.R1.final.txt
    cp -r "${previousassembly}"/*.fasta "${combinedassembly}"
    cp -r "${newreadsdirectory}/*.gz" "${newlongreads}"
    cd "${newlongreads}"
    tar zxvf *.gz && cp -r *.fasta "${combinedassembly}"
    (cd ..)
    cd "${combinedassembly}"
       cat *.fasta >> update_genome_assembly.fasta         
    updated_genome_assembly=""${combinedassembly}"/update_genome_assembly.fasta"
    blasr all_reads.fasta "${read_removal}" --best "${mapping}" --bam \
                                        "${updated_genome_assembly}".bam --unaligned \
                                                                 "${species}".unaligned.fasta
    mecat2ref -d "${unaligned_reads_assembly}" -r "${reference}" -o "${files%}".reference.sam \
                                            -w "${files%}"_intermediate -t 24 -n 20 -n 50 -m 2 -x 0
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
    echo "starting the genome coverage analysis"
    cd "${illumina_directory}"
    cp "${genome_assembly_fasta_file}" "${illumina_directory}"
    paste illumina.R1.txt illumina.R2.txt | while read -r col1 col2; \
    do 
        bowtie2-build "${genome_assembly_fasta_file}" "${species}"
        bowtie2 -t -x "${species}" -p "${threads}" --very-sensitive-local -1 "${col1}" \
                                 -2 "${col2}" -S "${species}".sam --no-unal --al-conc "${species}".aligned.fastq 
    done
    cp -r "${species}".sam "${polishedgenome}"
    (cd ..) 
    cd "${coverage}"
    for f in *.sam; do samtools view -bam "${f}".bam -sam "${f}"; done
    bamtools coverage -in "${f}".bam -out "${species}".bam.coverage.txt
    bamtools stats -in "${f}".bam -insert > "${species}".bam.stats.coverage.illumina.txt
    bamtools sort -in "${f}".bam -out "${f}".sorted.bam -n 500000 -mem 1024 
    cp -r "${f}".sorted.bam "${polishedgenome}"
    (cd ..)

    cd "${polishedgenome}"
    cp -r "${genome_assembly_fasta_file}" "${polishedgenome}"
    if [[ $polishoption == "jasper" ]]
    then 
    wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz 
        if [[ ! -f "jasper-1.0.3.tar.gz" ]]
        then
            echo "jasper file didnt downloaded and retrying"
            wget https://github.com/alguoo314/JASPER/releases/download/v1.0.3/jasper-1.0.3.tar.gz
        elif [[ -f "jasper-1.0.3.tar.gz" ]]
        then 
            echo "file downloaded and processing for the installation"
        fi 
        tar zxvf jasper-1.0.3.tar.gz
        cd jasper-1.0.3
        ./configure --prefix=$PWD
        make install
        export PYTHONPATH=$(pwd)
        (cd ..)
    export PATH="${polishedgenome}"/jasper-1.0.3:$PATH
    jasper.sh -t "${threads}" "${genome_assembly_fasta_file}" \
            -r "${illumina_directory}"/final.R1.fastq "${illumina_directory}"/final.R1.fastq
    else 
    java -Xmx100g -jar pilon-1.28.jar --genome "${genome_assembly_fasta_file}" \
                    --frags "${f}".sorted.bam --output "${species}".polished_genome \
                                --threads 60 --changes --fix all --mindepth 10
    fi 
    (cd ..)
    echo "processing finished"
else 
 printf "%s${closeuptext}"   
fi     
