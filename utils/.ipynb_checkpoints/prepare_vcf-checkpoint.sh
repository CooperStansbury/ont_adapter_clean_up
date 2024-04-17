INPUT=${1}
OUTPUT=${2}
THREADS=${3}

echo "INPUT:"  "${INPUT}"
echo "OUTPUT:" "${OUTPUT}"
echo "THREADS:" "${THREADS}"

# bcftools view -r 2 --min-ac 1 -g ^miss -m 2 -M 2 -O z --threads "${THREADS}" "${INPUT}" | bcftools annotate -x INFO,^FORMAT/GT -O z -o "${OUTPUT}"
bcftools view --min-ac 1 -g ^miss -m 2 -M 2 -O z --threads "${THREADS}" "${INPUT}" | bcftools annotate -x INFO,^FORMAT/GT -O z -o "${OUTPUT}"

