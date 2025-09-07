#!/bin/bash
# Script to run Java codebase for inferring species trees from gene trees and compute RF scores.
# Supports multiple methods (base, weighted_2_terminal, weighted_3_terminal) and handles complex Java command.

# JAVA="/home/aaniksahaa/.p2/pool/plugins/org.eclipse.justj.openjdk.hotspot.jre.full.linux.x86_64_21.0.7.v20250502-0916/jre/bin/java"

JAVA="java"

PROJECT_ROOT="/mnt/H/Research/STELAR-extension/STELAR"

JAVA_CMD="${JAVA} \
-Dfile.encoding=UTF-8 \
-Dstdout.encoding=UTF-8 \
-Dstderr.encoding=UTF-8 \
-classpath ${PROJECT_ROOT}/main:${PROJECT_ROOT}/DynaDup/main.jar:${PROJECT_ROOT}/main.jar:${PROJECT_ROOT}/mgd.jar:${PROJECT_ROOT}/miscellaneous-scripts/others/phylonet_v2_4.jar:${PROJECT_ROOT}/STELAR.jar \
-XX:+ShowCodeDetailsInExceptionMessages phylonet.coalescent.MGDInference_DP"


# Define outer folders (taxa numbers)
# folders=("11-taxon")
# folders=("15-taxon")  
# folders=("37-taxon")  
# folders=("48-taxon")
# folders=("100-taxon")
folders=("200-taxon")

fresh=1  # Set to 0 to skip existing output files, 1 to overwrite

# Define methods to run
# methods=("base" "unrooted" "weighted_2_terminal" "weighted_3_terminal")
methods=("unrooted")

# Create dynamic folder structure for summary files
# Join folders array with hyphens
folders_joined=$(IFS='-'; echo "${folders[*]}")

# Join methods array with hyphens
methods_joined=$(IFS='-'; echo "${methods[*]}")

# Create timestamp in hr-min-dd-mm-yy format
timestamp=$(date +"%H-%M-%d-%m-%y")

# Create the summary directory structures
summary_main_dir="./datasets/summary_new/${folders_joined}/${methods_joined}/${timestamp}"
summary_backup_dir="./results/${folders_joined}/${methods_joined}/${timestamp}"

mkdir -p "$summary_main_dir"
mkdir -p "$summary_backup_dir"

echo "Summary results will be saved to: $summary_main_dir"
echo "Backup summary results will be saved to: $summary_backup_dir"

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

# Create global CSV files for summary_new in the new nested structure
timediffcsv_main="${summary_main_dir}/average_diffs.csv"
rf_csv_main="${summary_main_dir}/average_rfs.csv"

# Create backup CSV files
timediffcsv_backup="${summary_backup_dir}/average_diffs.csv"
rf_csv_backup="${summary_backup_dir}/average_rfs.csv"

# Initialize CSV files
printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_diff" > "$timediffcsv_main"
printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_rf" > "$rf_csv_main"
printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_diff" > "$timediffcsv_backup"
printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_rf" > "$rf_csv_backup"

# Loop through each taxa folder
for folder in "${folders[@]}"; do
    echo "Processing folder: $folder"

    # Create local CSV files for each taxa number (keep original location)
    mkdir -p "./datasets/$folder/summary_new"
    timediffcsv="./datasets/$folder/summary_new/average_diffs.csv"
    rf_csv="./datasets/$folder/summary_new/average_rfs.csv"
    printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_diff" > "$timediffcsv"
    printf "%s,%s,%s,%s\n" "summary_new_method" "folder" "inner_folder" "average_rf" > "$rf_csv"

    # Get inner folders and number of replicates
    IFS=' ' read -r -a inner_folders <<< "${innerFolderNames[$folder]}"
    R=${replicates[$folder]}

    # Loop through each inner folder
    for inner_folder in "${inner_folders[@]}"; do
        echo "Processing inner folder: $inner_folder"

        # Initialize arrays to store sums and counts for each method
        declare -A sum_diffs sum_rfs count_diffs
        for method in "${methods[@]}"; do
            sum_diffs[$method]=0
            sum_rfs[$method]=0
            count_diffs[$method]=0
        done

        # Loop through replicates
        for ((j=1; j<=R; j++)); do
            echo "Processing replicate R$j"

            # Define paths
            gt_folder="$inner_folder/R$j"
            input="./datasets/$folder/$gt_folder/all_gt.tre.rooted"
            true_tree="./datasets/$folder/true_tree_trimmed"

            # Special handling for certain folders
            if [[ "$folder" == "100-taxon" || "$folder" == "200-taxon" || "$folder" == "500-taxon" ]]; then
                true_tree="./datasets/$folder/true-species-trees/R$j/sp-cleaned"
            elif [[ "$folder" == "11-taxon-new" ]]; then
                true_tree="./datasets/$folder/higher-ILS/true-speciestrees/R$j.true.tre"
            elif [[ "$folder" == "1000-taxon" ]]; then
                if [[ $j -lt 10 ]]; then
                    gt_folder="$inner_folder/0$j"
                    true_tree="./datasets/$folder/true-species-trees/0$j/s_tree.trees"
                else
                    gt_folder="$inner_folder/$j"
                    true_tree="./datasets/$folder/true-species-trees/$j/s_tree.trees"
                fi
                input="./datasets/$folder/$gt_folder/stelar_inputs/stelar_input.tre"
            fi

            # Create output directory (keep original location)
            mkdir -p "./datasets/$folder/$gt_folder/stelar_outputs"

            # Loop through each method (innermost loop)
            for method in "${methods[@]}"; do
                echo "Processing method: $method for replicate R$j"

                # Determine input file and Java method based on method type
                if [[ "$method" == stelar_* ]]; then
                    # For stelar_x methods, use base as Java method and stelar_input_cleaned_x.tre as input
                    stelar_suffix="${method#stelar_}"  # Extract x from stelar_x
                    
                    # Handle different file naming conventions for 100-taxon folder
                    if [[ "$folder" == "100-taxon" ]]; then
                        method_input="./datasets/$folder/$gt_folder/stelar_inputs/stelar_input_sanitized${stelar_suffix}.tre"
                    else
                        method_input="./datasets/$folder/$gt_folder/stelar_inputs/stelar_input_cleaned_${stelar_suffix}.tre"
                    fi
                    java_method="base"
                else
                    # For regular methods, use the original input and method
                    method_input="$input"
                    java_method="$method"
                fi

                # Special fallback handling for unrooted method
                if [[ "$method" == "unrooted" && ! -f "$method_input" ]]; then
                    # Handle different file naming conventions for 100-taxon folder in fallback
                    if [[ "$folder" == "100-taxon" ]]; then
                        fallback_input="./datasets/$folder/$gt_folder/stelar_inputs/stelar_input_sanitizedMV.tre"
                    else
                        fallback_input="./datasets/$folder/$gt_folder/stelar_inputs/stelar_input_cleaned_MV.tre"
                    fi
                    if [[ -f "$fallback_input" ]]; then
                        echo "Original input not found for unrooted method, using fallback: $fallback_input"
                        method_input="$fallback_input"
                    fi
                fi

                # Define output path (keep original location)
                output="./datasets/$folder/$gt_folder/stelar_outputs/stelar_output_${method}.tre"

                # Skip if output exists and fresh is 0
                if [[ -f "$output" && $fresh -eq 0 ]]; then
                    echo "Output exists, skipping: $output"
                    continue
                fi

                # Ensure method input file exists
                if [[ ! -f "$method_input" ]]; then
                    echo "Error: Method input file does not exist: $method_input"
                    continue
                fi

                # Print the input file path before running
                echo "Running on input file: $method_input"

                # Run Java codebase with method and measure time
                START=$(date +%s.%N)
                $JAVA_CMD -m "$java_method" -i "$method_input" -o "$output"
                END=$(date +%s.%N)
                DIFF=$(echo "$END - $START" | bc)
                echo "Time taken for method $method: $DIFF seconds"

                # Calculate RF distance
                if [[ -f "$output" && -f "$true_tree" ]]; then
                    index_to_extract=2
                    tuple=$(python3 ./RF/getFpFn.py -e "$output" -t "$true_tree")
                    RFdistance=$(echo "$tuple" | cut -d',' -f"$index_to_extract" | tr -d '()')
                    echo "RF Distance for method $method: $RFdistance"

                    # Update sums and count for this method
                    sum_diffs[$method]=$(echo "${sum_diffs[$method]} + $DIFF" | bc)
                    sum_rfs[$method]=$(echo "${sum_rfs[$method]} + $RFdistance" | bc)
                    ((count_diffs[$method]++))
                else
                    echo "Error: Output or true tree file missing for RF calculation"
                fi
            done
        done

        # Calculate and store averages for each method
        for method in "${methods[@]}"; do
            if [[ ${count_diffs[$method]} -gt 0 ]]; then
                average_diff=$(echo "${sum_diffs[$method]} / ${count_diffs[$method]}" | bc -l)
                average_rf=$(echo "${sum_rfs[$method]} / ${count_diffs[$method]}" | bc -l)
                printf "Average DIFF for %s (%s): %.6f\n" "$inner_folder" "$method" "$average_diff"
                printf "Average RF for %s (%s): %.6f\n" "$inner_folder" "$method" "$average_rf"

                # Write to local CSV (keep original location)
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_diff" >> "$timediffcsv"
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_rf" >> "$rf_csv"

                # Write to global CSV (new nested location)
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_diff" >> "$timediffcsv_main"
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_rf" >> "$rf_csv_main"

                # Write to backup CSV
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_diff" >> "$timediffcsv_backup"
                printf "%s,%s,%s,%.6f\n" "$method" "$folder" "$inner_folder" "$average_rf" >> "$rf_csv_backup"
            fi
        done

        # Add newline to CSVs
        printf "\n" >> "$timediffcsv"
        printf "\n" >> "$rf_csv"
        printf "\n" >> "$timediffcsv_main"
        printf "\n" >> "$rf_csv_main"
        printf "\n" >> "$timediffcsv_backup"
        printf "\n" >> "$rf_csv_backup"
    done
done

echo "Processing complete."
echo "Individual results saved in ./datasets/"
echo "Summary results saved in: $summary_main_dir"
echo "Backup summary results saved in: $summary_backup_dir"


