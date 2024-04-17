INPUT=${1}
ENZYME="NlaIII"
REF=${2}
OUTPUT="all"

pore-c-py digest "${INPUT}" "${ENZYME}" \
    | samtools fastq -T '*' \
    | minimap2 -ay -t 8 -x map-ont "${REF}" - \
    | pore-c-py annotate - "${OUTPUT}" --monomers --stdout --summary --chromunity \
    | tee "${OUTPUT}.ns.bam" \
    | samtools sort --write-index -o "${OUTPUT}.cs.bam" -
samtools index "${OUTPUT}.ns.bam"
