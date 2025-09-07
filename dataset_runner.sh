#!/bin/bash
# Dataset Runner Script for STELAR-MP
# Runs STELAR-MP on entire dataset folders with validation and error handling
# Based on script.sh template but adapted for STELAR-MP project
#
# Usage: ./dataset_runner.sh [--fresh]
#   --fresh: Overwrite existing outputs and start fresh
#   (default): Skip existing outputs and append to existing CSVs
#
# Example: ./dataset_runner.sh        # Skip existing outputs (default)
#          ./dataset_runner.sh --fresh # Run with fresh outputs

# =============================================================================
# CONFIGURATION SECTION
# =============================================================================

# Define the phylogenetic datasets directory
PHYLO_DIR="/home/aaniksahaa/research/phylogeny"  # UPDATE THIS PATH

# Define computation modes to run - modify as needed
COMPUTATION_MODES=("GPU_PARALLEL")
# COMPUTATION_MODES=("CPU_PARALLEL")  # Use this for single mode testing

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# =============================================================================
# DATASET CONFIGURATION
# =============================================================================

# Define outer folders (taxa numbers) - modify as needed
# folders=("11-taxon" "15-taxon" "37-taxon" "48-taxon" "100-taxon" "200-taxon")

# Define outer folders (taxa numbers)
# folders=("11-taxon")
# folders=("15-taxon")  
# folders=("37-taxon")  
# folders=("48-taxon")
# folders=("100-taxon")
folders=("200-taxon")


# Define inner folder names for each taxa number as space-separated strings
declare -A innerFolderNames
innerFolderNames["11-taxon"]="estimated_Xgenes_strongILS/estimated_5genes_strongILS estimated_Xgenes_strongILS/estimated_15genes_strongILS estimated_Xgenes_strongILS/estimated_25genes_strongILS estimated_Xgenes_strongILS/estimated_50genes_strongILS estimated_Xgenes_strongILS/estimated_100genes_strongILS"
innerFolderNames["15-taxon"]="100gene-100bp/estimated-genetrees 100gene-1000bp/estimated-genetrees 100gene-true 1000gene-100bp/estimated-genetrees 1000gene-1000bp/estimated-genetrees 1000gene-true"
innerFolderNames["37-taxon"]="estimated-genetrees/1X-200-500 estimated-genetrees/1X-200-1000 estimated-genetrees/1X-400-500 estimated-genetrees/1X-400-1000 estimated-genetrees/1X-800-500 estimated-genetrees/1X-800-1000 estimated-genetrees/2X-200-500"
innerFolderNames["48-taxon"]="estimated-genetrees/1X-25-500 estimated-genetrees/1X-50-500 estimated-genetrees/1X-100-500 estimated-genetrees/1X-200-500 estimated-genetrees/1X-500-500 estimated-genetrees/1X-1000-500 estimated-genetrees/2X-1000-500"
innerFolderNames["48-taxon-latest"]="estimated_genetrees/1X-1000-500"
innerFolderNames["100-taxon"]="inner100"
innerFolderNames["200-taxon"]="inner200"
innerFolderNames["500-taxon"]="model.500.2000000.0.000001"
innerFolderNames["biological"]="nuclear"
innerFolderNames["mammalian"]="424genes"
innerFolderNames["amniota"]="aa nt"

# Define number of replicates for each folder
declare -A replicates
replicates["11-taxon"]=20
replicates["15-taxon"]=10
replicates["37-taxon"]=20
replicates["48-taxon"]=20
replicates["48-taxon-latest"]=20
replicates["100-taxon"]=8
replicates["200-taxon"]=10
replicates["500-taxon"]=1
replicates["biological"]=1
replicates["mammalian"]=1
replicates["amniota"]=1

# =============================================================================
# VALIDATION FUNCTIONS
# =============================================================================

validate_phylo_dir() {
    if [ ! -d "$PHYLO_DIR" ]; then
        echo -e "${RED}Error: PHYLO_DIR '$PHYLO_DIR' does not exist.${NC}"
        echo "Please update the PHYLO_DIR variable at the top of this script."
        exit 1
    fi
    
    if [ ! -d "$PHYLO_DIR/datasets" ]; then
        echo -e "${RED}Error: datasets folder not found in '$PHYLO_DIR'.${NC}"
        exit 1
    fi
    
    if [ ! -d "$PHYLO_DIR/RF" ]; then
        echo -e "${RED}Error: RF folder not found in '$PHYLO_DIR'.${NC}"
        exit 1
    fi
}

validate_computation_modes() {
    local valid_modes=("CPU_SINGLE" "CPU_PARALLEL" "GPU_PARALLEL")
    
    for mode in "${COMPUTATION_MODES[@]}"; do
        if [[ ! " ${valid_modes[@]} " =~ " ${mode} " ]]; then
            echo -e "${RED}Error: Invalid computation mode '$mode'${NC}"
            echo "Valid modes: ${valid_modes[*]}"
            exit 1
        fi
    done
}

validate_stelar_binaries() {
    if [ ! -f "target/stelar-mp-1.0-SNAPSHOT.jar" ] || [ ! -d "bin" ]; then
        echo -e "${RED}Error: STELAR-MP binaries not found. Please run build.sh first.${NC}"
        exit 1
    fi
}

validate_cuda_library() {
    # Check if any computation mode requires GPU
    for mode in "${COMPUTATION_MODES[@]}"; do
        if [[ "$mode" == "GPU_PARALLEL" ]]; then
            if [ ! -f "$(pwd)/cuda/libweight_calc.so" ]; then
                echo -e "${RED}Error: CUDA library not found. GPU_PARALLEL mode requires libweight_calc.so${NC}"
                exit 1
            fi
            break
        fi
    done
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

print_header() {
    echo "==============================================="
    echo "=== STELAR-MP Dataset Runner ==="
    echo "==============================================="
    echo "PHYLO_DIR: $PHYLO_DIR"
    echo "Computation modes: ${COMPUTATION_MODES[*]}"
    echo "Fresh run: $fresh"
    echo "Folders to process: ${folders[*]}"
    echo "==============================================="
    echo
}

print_debug_info() {
    echo -e "${YELLOW}Debug Information:${NC}"
    echo "Current directory: $(pwd)"
    echo "Library path: $(pwd)/cuda"
    echo "Library exists: $(if [ -f "$(pwd)/cuda/libweight_calc.so" ]; then echo "Yes"; else echo "No"; fi)"
    echo "Library permissions: $(ls -l "$(pwd)/cuda/libweight_calc.so" 2>/dev/null || echo "Not found")"
    echo
}

run_stelar_mp() {
    local input_file=$1
    local output_file=$2
    local comp_mode=$3
    
    # Create output directory if it doesn't exist
    local output_dir=$(dirname "$output_file")
    if [ ! -d "$output_dir" ]; then
        mkdir -p "$output_dir"
    fi
    
    echo -e "${BLUE}Running STELAR-MP...${NC}"
    echo "Input: $input_file"
    echo "Output: $output_file"
    echo "Mode: $comp_mode"
    
    # Run STELAR-MP with proper library path
    java -Djava.library.path="$(pwd)/cuda" \
         -Djna.debug_load=false \
         -Djna.platform.library.path="$(pwd)/cuda" \
         -cp target/stelar-mp-1.0-SNAPSHOT.jar Main \
         -i "$input_file" -o "$output_file" -m "$comp_mode"
    
    return $?
}

calculate_rf_distance() {
    local estimated_tree=$1
    local true_tree=$2
    
    if [[ -f "$estimated_tree" && -f "$true_tree" ]]; then
        local index_to_extract=2
        local tuple=$(python3 "$PHYLO_DIR/RF/getFpFn.py" -e "$estimated_tree" -t "$true_tree" 2>/dev/null)
        if [ $? -eq 0 ]; then
            local rf_distance=$(echo "$tuple" | cut -d',' -f"$index_to_extract" | tr -d '()')
            echo "$rf_distance"
        else
            echo "ERROR"
        fi
    else
        echo "ERROR"
    fi
}

# =============================================================================
# MAIN PROCESSING FUNCTION
# =============================================================================

process_datasets() {
    local fresh_run=$1
    
    # Create results directory (no timestamp)
    local results_dir="./results"
    mkdir -p "$results_dir"
    
    # Initialize summary files
    local time_csv="$results_dir/execution_times.csv"
    local rf_csv="$results_dir/rf_distances.csv"
    
    # Handle CSV files based on fresh_run flag
    if [[ $fresh_run -eq 1 ]]; then
        # Fresh run: overwrite existing CSV files
        echo "Fresh run: Creating new CSV files"
        printf "%s,%s,%s,%s,%s\n" "folder" "inner_folder" "replicate" "execution_time" "computation_mode" > "$time_csv"
        printf "%s,%s,%s,%s,%s\n" "folder" "inner_folder" "replicate" "rf_distance" "computation_mode" > "$rf_csv"
    else
        # Non-fresh run: append to existing CSV files or create if they don't exist
        if [[ ! -f "$time_csv" ]]; then
            echo "Creating new execution times CSV file"
            printf "%s,%s,%s,%s,%s\n" "folder" "inner_folder" "replicate" "execution_time" "computation_mode" > "$time_csv"
        else
            echo "Appending to existing execution times CSV file"
        fi
        
        if [[ ! -f "$rf_csv" ]]; then
            echo "Creating new RF distances CSV file"
            printf "%s,%s,%s,%s,%s\n" "folder" "inner_folder" "replicate" "rf_distance" "computation_mode" > "$rf_csv"
        else
            echo "Appending to existing RF distances CSV file"
        fi
    fi
    
    echo -e "${GREEN}Results will be saved to: $results_dir${NC}"
    echo
    
    # Loop through each taxa folder
    for folder in "${folders[@]}"; do
        echo -e "${YELLOW}Processing folder: $folder${NC}"
        
        # Check if folder exists
        if [ ! -d "$PHYLO_DIR/datasets/$folder" ]; then
            echo -e "${RED}Warning: Folder '$PHYLO_DIR/datasets/$folder' does not exist, skipping...${NC}"
            continue
        fi
        
        # Get inner folders and number of replicates
        IFS=' ' read -r -a inner_folders <<< "${innerFolderNames[$folder]}"
        local R=${replicates[$folder]}
        
        # Loop through each inner folder
        for inner_folder in "${inner_folders[@]}"; do
            echo "  Processing inner folder: $inner_folder"
            
            # Loop through replicates
            for ((j=1; j<=R; j++)); do
                echo "    Processing replicate R$j"
                
                # Define paths
                local gt_folder="$inner_folder/R$j"
                local input="$PHYLO_DIR/datasets/$folder/$gt_folder/all_gt.tre.rooted"
                local true_tree="$PHYLO_DIR/datasets/$folder/true_tree_trimmed"
                
                # Special handling for certain folders (adapt as needed)
                if [[ "$folder" == "100-taxon" || "$folder" == "200-taxon" ]]; then
                    true_tree="$PHYLO_DIR/datasets/$folder/true-species-trees/R$j/sp-cleaned"
                fi
                
                # Check if input file exists
                if [ ! -f "$input" ]; then
                    echo "    Warning: Input file does not exist: $input"
                    continue
                fi
                
                # Loop through computation modes (innermost loop)
                for comp_mode in "${COMPUTATION_MODES[@]}"; do
                    echo "      Processing computation mode: $comp_mode"
                    
                    # Define output path
                    local output_dir="$PHYLO_DIR/datasets/$folder/$gt_folder/stelar_outputs"
                    local output="$output_dir/stelar_output_${comp_mode}.tre"
                    
                    # Skip if output exists and fresh_run is 0
                    if [[ -f "$output" && $fresh_run -eq 0 ]]; then
                        echo "      Output exists, skipping: $output"
                        continue
                    fi
                    
                    # Run STELAR-MP and measure time
                    local start_time=$(date +%s.%N)
                    if run_stelar_mp "$input" "$output" "$comp_mode"; then
                        local end_time=$(date +%s.%N)
                        local execution_time=$(echo "$end_time - $start_time" | bc)
                        
                        echo "      Execution time: $execution_time seconds"
                        
                        # Calculate RF distance
                        local rf_distance=$(calculate_rf_distance "$output" "$true_tree")
                        if [ "$rf_distance" != "ERROR" ]; then
                            echo "      RF Distance: $rf_distance"
                        else
                            echo "      RF Distance: Could not calculate"
                            rf_distance="N/A"
                        fi
                        
                        # Write to CSV files
                        printf "%s,%s,%s,%.6f,%s\n" "$folder" "$inner_folder" "R$j" "$execution_time" "$comp_mode" >> "$time_csv"
                        printf "%s,%s,%s,%s,%s\n" "$folder" "$inner_folder" "R$j" "$rf_distance" "$comp_mode" >> "$rf_csv"
                        
                        echo -e "      ${GREEN}✓ Completed successfully${NC}"
                    else
                        echo -e "      ${RED}✗ Execution failed${NC}"
                    fi
                done
            done
        done
        echo
    done
}

# =============================================================================
# MAIN SCRIPT
# =============================================================================

# Parse command line arguments
fresh=0  # Default: skip existing outputs (not fresh)
if [[ "$1" == "--fresh" ]]; then
    fresh=1
    echo "Running in fresh mode: will overwrite existing outputs"
else
    echo "Running in resume mode: will skip existing outputs and append to CSVs"
fi

# Print header
print_header

# Validate everything
echo -e "${YELLOW}Performing validation checks...${NC}"
validate_phylo_dir
validate_computation_modes
validate_stelar_binaries
validate_cuda_library
echo -e "${GREEN}✓ All validation checks passed${NC}"
echo

# Print debug information
print_debug_info

# Process datasets
echo -e "${YELLOW}Starting dataset processing...${NC}"
process_datasets "$fresh"

echo -e "${GREEN}Dataset processing complete!${NC}"
echo "Check the results directory for execution times and RF distances."
