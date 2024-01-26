# long_read_genome_assembly_pangraphs
a complete workflow that can be dockerized for the long read assembly, it allows for the genome assembly update as well as it allows for the assembly from the start. if you have the illumina reads it allows for the genome mapping also. long read dockerized pipeline it will easier to execute and also better for the docker application. it will even make your genome browser tracks. 

You can select the options to make your genome browsers tracks and direct plugins into the jbrowse for the visualization or snp calling tracts for the visualization. 

Finished code first release: 2023-10-20, next update will be with the genome visualization tracts. \
Finished update: 2023-10-21, major changes to the code, added associative arrays and added the support for the multiple genome polishing such as pilon and jasper and also the coverage analysis \
Next update: adding the docker, portainer, jbrowse, multiple visualization and cluster profiling tools directly from the slurm or pbs. \
2024-1-14: \
Another update launch this week: adding the support for the complete docker and the instance based launch. 
Another update launch this week: adding the support for the PacBio HiFi reads also. 
2024-1-26 \
Another update coming with support for the annotations. \
More support for oxford nanopore in the next release. \
Adding support for the decalarative arrays and you will get more information. \
This code is constantly updated and also a docker withe several integrations is in preparation which will allow you to run the dockerized analysis. 
```
# block of code is here for viewing
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
fi
 echo "conda environment has been created"
if  [[ $species ]]  && [[ $first == "yes" ]] && 
                 [[ $assembler == "canu" ]] && [[ $alignment_tracks == "yes" ]] &&
                 [[ $mapping == "yes" ]] && [[ $coverage == "yes" ]] && [[ $polish == "yes" ]]
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
        echo "$file" > illumina.R1.txt
    done
    for file in "${illumina_directory}"/*.R2.fastq
    do 
        echo "$file" > illumina.R2.txt
    done
    cd "${longreads}"
    genome_contamant=${read_removal}
        blasr all_reads.fasta "${genome_contamant}" --best "${mapping}" --bam \
                                      "${species}".bam --unaligned "${species}".unaligned.fasta
        echo "finished cleaning of the reads"
    cp -r "${species}".unaligned.fasta "${genome_assembly}"
    (cd ..)

    cd "${genome_assembly}"
        echo "working in the present genome assembly directory"
    canu gridOptions="--time=24:00:00" -corMhapSensitivity="${sensitivity}" \
                                                        -p "{$species}" \
                                                        -d "$(pwd)/genome_assembly" \
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
    java -Xmx100g -jar pilon-1.28.jar --genome "${genome_assembly_fasta_file}" \
                    --frags "${f}".sorted.bam --output "${species}".polished_genome \
                                --threads 60 --changes --fix all --mindepth 10
    (cd ..)
    echo "processing finished"
fi
```

Gaurav Sablok, \
Academic Staff Member,\
Bioinformatics, \
Institute for Biochemistry and Biology, \
University of Potsdam,Potsdam,\
Germany
