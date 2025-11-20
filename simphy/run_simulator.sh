#!/usr/bin/env bash
set -euo pipefail

# Default values
sb="0.000001"
spmin="500000"
spmax="1500000"
out_dir=""
data_base_dir="data"  # Default data directory (backward compatible)
replicates="10"  # Default number of replicates
# Required (no defaults)
taxa_num=""
gene_trees=""

print_usage() {
  cat <<EOF
Usage: $0 -t <taxa_num> -g <gene_trees> [options]

Required:
  -t, --taxa_num    Number of taxa (e.g. 1000)
  -g, --gene_trees  Number of gene trees (e.g. 500)

Optional:
  -o, --out_dir     Output directory (if omitted, auto-built from params)
  -d, --data_dir    Base data directory (default: ${data_base_dir})
      --replicates  Number of replicates (default: ${replicates})
      --sb          Substitution/birthrate parameter (default: ${sb})
      --spmin       Population size minimum (default: ${spmin})
      --spmax       Population size maximum (default: ${spmax})
  -h, --help        Show this help and exit

Example:
  $0 -t 5000 -g 1000 --replicates 5 --sb 0.000001 --spmin 500000 --spmax 1500000
EOF
}

# Use GNU getopt for long options
OPTS=$(getopt -o t:g:o:d:h --long taxa_num:,gene_trees:,out_dir:,data_dir:,replicates:,sb:,spmin:,spmax:,help -n 'run_simulator.sh' -- "$@")
if [ $? != 0 ] ; then
  echo "Failed parsing options." >&2
  exit 1
fi
eval set -- "$OPTS"

while true; do
  case "$1" in
    -t|--taxa_num) taxa_num="$2"; shift 2 ;;
    -g|--gene_trees) gene_trees="$2"; shift 2 ;;
    -o|--out_dir) out_dir="$2"; shift 2 ;;
    -d|--data_dir) data_base_dir="$2"; shift 2 ;;
    --replicates) replicates="$2"; shift 2 ;;
    --sb) sb="$2"; shift 2 ;;
    --spmin) spmin="$2"; shift 2 ;;
    --spmax) spmax="$2"; shift 2 ;;
    -h|--help) print_usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Unexpected argument to run_simulator.sh: $1, Internal error while parsing options!"; exit 1 ;;
  esac
done

# Validate required args
if [ -z "${taxa_num}" ] || [ -z "${gene_trees}" ]; then
  echo "Error: --taxa_num (-t) and --gene_trees (-g) are required."
  print_usage
  exit 1
fi

# Ensure data base directory exists
mkdir -p "${data_base_dir}"

# Construct output folder if not provided
if [ -z "${out_dir}" ]; then
  out_dir="${data_base_dir}/t_${taxa_num}_g_${gene_trees}_sb_${sb}_spmin_${spmin}_spmax_${spmax}"
else
  # If out_dir is provided, use it as-is (could be absolute or relative)
  out_dir="${out_dir%/}"  # strip trailing slash
fi

mkdir -p "${out_dir}"

echo "Running with parameters:"
echo "  taxa_num    = ${taxa_num}"
echo "  gene_trees  = ${gene_trees}"
echo "  replicates  = ${replicates}"
echo "  sb          = ${sb}"
echo "  spmin       = ${spmin}"
echo "  spmax       = ${spmax}"
echo "  data_base_dir = ${data_base_dir}"
echo "  out_dir     = ${out_dir}"
echo ""

# sb = Speciation rate
# spmin Minimum population size
# spmax Maximum population size
# ld Loss rate
# lb Duplication rate
# lt HGT rate
# rs Number of replicates
# rl Number of gene trees
# rg Number of gene trees per replicate
# o Output directory
# sp Effective population size is uniform between spmin and spmax
# su Substitution rate is log-normal with mean -17.27461 and standard deviation 0.6931472
# sg Generation time is fixed at 1
# sl Number of taxa
# st Species tree height is log-normal with mean 16.2 and standard deviation 1
# om Output mode

# v Verbosity level
# od Output directory
# op Output mode
# oc Output mode
# on Output mode
# cs Seed

# Run SimPhy (adjust path to simphy_lnx64 if necessary)
./simphy_lnx64 \
  -sb f:${sb} \
  -ld f:0 \
  -lb f:0 \
  -lt f:0 \
  -rs ${replicates} \
  -rl f:${gene_trees} \
  -rg 1 \
  -o ${out_dir} \
  -sp u:${spmin},${spmax} \
  -su ln:-17.27461,0.6931472 \
  -sg f:1 \
  -sl f:${taxa_num} \
  -st ln:16.2,1 \
  -om 1 \
  -v 2 \
  -od 1 \
  -op 1 \
  -oc 1 \
  -on 1 \
  -cs 42


# Canonicalize directory names by removing leading zeros
echo "Canonicalizing replicate directory names..."
for dir in "${out_dir}"/*/; do
  if [ -d "$dir" ]; then
    dirname=$(basename "$dir")
    # Check if directory name is numeric (with possible leading zeros)
    if [[ "$dirname" =~ ^0*([1-9][0-9]*|0)$ ]]; then
      canonical_name="${BASH_REMATCH[1]}"
      if [ "$dirname" != "$canonical_name" ]; then
        echo "Renaming ${out_dir}/${dirname} to ${out_dir}/${canonical_name}"
        mv "${out_dir}/${dirname}" "${out_dir}/${canonical_name}"
      fi
    fi
  fi
done

# Run concat_gene_trees.py on all replicates
echo "Running concat_gene_trees.py on all replicates..."
for i in $(seq 1 ${replicates}); do
  if [ -d "${out_dir}/${i}" ]; then
    echo "Processing replicate ${i}..."
    python concat_gene_trees.py "${out_dir}/${i}"
  else
    echo "Warning: replicate directory ${out_dir}/${i} not found. Skipping concat_gene_trees.py."
  fi
done

# Reorganize all trees in dataset (run once after all concat operations)
echo "Reorganizing all trees in dataset..."
python reorganize_trees.py "${out_dir}"

echo "Done. Output in ${out_dir}"
