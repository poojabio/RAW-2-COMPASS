#!/bin/bash
#kallisto 0.44.0
#
#Usage: kallisto <CMD> [arguments] ..
#
#Where <CMD> can be one of:
#
#    index         Builds a kallisto index
#    quant         Runs the quantification algorithm
#    pseudo        Runs the pseudoalignment step
#    h5dump        Converts HDF5-formatted results to plaintext
#    inspect       Inspects and gives information about an index
#    version       Prints version information
#    cite          Prints citation information
#
#Running kallisto <CMD> without arguments prints usage information for <CMD>



## ========================== ##

#!/bin/bash
#kallisto 0.44.0

## getting kallisto  + setup to access using just "kallisto"
wget https://github.com/pachterlab/kallisto/releases/download/v0.52.0/kallisto_linux-v0.52.0.tar.gz
tar -zxvf kallisto_linux-v0.52.0.tar.gz
export PATH=$PATH:$(pwd)/kallisto_linux-v0.52.0

##the assembly visits the rest api and looks through retrival of the JSON object matching your species
## for the species you have a match of it will then take trancript assembly
assembly=$(curl -s "https://rest.ensembl.org/info/assembly/${1}?content-type=application/json" \
    | grep -o '"assembly_name":"[^"]*"' | cut -d'"' -f4)

species_cap="$(tr '[:lower:]' '[:upper:]' <<< ${1:0:1})${1:1}"
url="https://ftp.ensembl.org/pub/release-${2}/fasta/${1}/cdna/${species_cap}.${assembly}.cdna.all.fa.gz"

ref_fasta="${1}_cdna.fa.gz"
wget -O "$ref_fasta" "$url"

index_file="${1}_index.idx"
kallisto index -i "$index_file" "$ref_fasta"

##identifying all fastq files
files=($(ls *.fastq.gz | sort))

##for every pair, identify f1 and f2
for (( i=0; i<${#files[@]}; i+=2 )); do
    f1="${files[i]}"
    f2="${files[i+1]}"

    echo "Sample: $f1 , $f2"
    kallisto quant -i "$index_file" -o "${f1}_quant" "$f1" "$f2"
done

echo "All samples complete!"