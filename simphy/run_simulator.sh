#!/usr/bin/env bash
set -euo pipefail

# Default values
sb="0.000001"
spmin="500000"
spmax="1500000"
out_dir=""
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
      --sb          Substitution/birthrate parameter (default: ${sb})
      --spmin       Population size minimum (default: ${spmin})
      --spmax       Population size maximum (default: ${spmax})
  -h, --help        Show this help and exit

Example:
  $0 -t 5000 -g 1000 --sb 0.000001 --spmin 500000 --spmax 1500000
EOF
}

# Use GNU getopt for long options
OPTS=$(getopt -o t:g:o:h --long taxa_num:,gene_trees:,out_dir:,sb:,spmin:,spmax:,help -n 'run_simulator.sh' -- "$@")
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
    --sb) sb="$2"; shift 2 ;;
    --spmin) spmin="$2"; shift 2 ;;
    --spmax) spmax="$2"; shift 2 ;;
    -h|--help) print_usage; exit 0 ;;
    --) shift; break ;;
    *) echo "Internal error while parsing options"; exit 1 ;;
  esac
done

# Validate required args
if [ -z "${taxa_num}" ] || [ -z "${gene_trees}" ]; then
  echo "Error: --taxa_num (-t) and --gene_trees (-g) are required."
  print_usage
  exit 1
fi

# Construct output folder if not provided
if [ -z "${out_dir}" ]; then
  out_dir="data/t_${taxa_num}_g_${gene_trees}_sb_${sb}_spmin_${spmin}_spmax_${spmax}"
else
  out_dir="${out_dir%/}"  # strip trailing slash
fi

mkdir -p "${out_dir}"

echo "Running with parameters:"
echo "  taxa_num  = ${taxa_num}"
echo "  gene_trees= ${gene_trees}"
echo "  sb        = ${sb}"
echo "  spmin     = ${spmin}"
echo "  spmax     = ${spmax}"
echo "  out_dir   = ${out_dir}"
echo ""

# Run SimPhy (adjust path to simphy_lnx64 if necessary)
./simphy_lnx64 \
  -sb f:${sb} \
  -ld f:0 \
  -lb f:0 \
  -lt f:0 \
  -rs 1 \
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

# Run concat_gene_trees.py on replicate 1 (keep behavior)
if [ -d "${out_dir}/1" ]; then
  python concat_gene_trees.py "${out_dir}/1"
else
  echo "Warning: replicate directory ${out_dir}/1 not found. Skipping concat_gene_trees.py."
fi

# Reorganize all trees in dataset
python reorganize_trees.py "${out_dir}"

echo "Done. Output in ${out_dir}"
