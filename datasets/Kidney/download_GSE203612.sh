#!/bin/bash

# Download script for GSE203612 - Pan-cancer single-cell RNA-seq with endometrial cancer samples
# Dataset: Cancer cell states recur across tumor types and form specific interactions with the tumor microenvironment
# Publication: Nature Genetics (2022)
# BioProject: PRJNA841591
# Samples: 31 experiments across 15 cancer types (including endometrial cancer)
# Platform: Illumina NextSeq 500

set -e  # Exit on any error

# Configuration
BASE_DIR=/gpfs/scratch/nk4167/
cd "$BASE_DIR"
PROJECT_ID="PRJNA841591"
GSE_ID="GSE203612"
OUTPUT_DIR="./GSE203612_data"
TEMP_DIR="./temp_${GSE_ID}"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_DIR"

echo "=== Downloading GSE203612 Pan-cancer Dataset ==="
echo "Focus: Endometrial cancer samples from pan-cancer study"
echo "BioProject: $PROJECT_ID"
echo "GEO Series: $GSE_ID"
echo "Output directory: $OUTPUT_DIR"
echo "Temporary directory: $TEMP_DIR"

# Step 1: Define run accessions for GSE203612/PRJNA841591
echo "Step 1: Setting up run accessions for $PROJECT_ID..."

# All 31 SRA run accessions from PRJNA841591
# Note: This pan-cancer study includes endometrial cancer among 15 cancer types
cat > "$TEMP_DIR/run_accessions.txt" << 'EOF'
SRR19360476
SRR19360477
SRR19360478
SRR19360481
SRR19360484
SRR19360485
SRR19360486
EOF

echo "Run accessions prepared (31 total):"
cat "$TEMP_DIR/run_accessions.txt"



RUN_COUNT=$(cat "$TEMP_DIR/run_accessions.txt" | wc -l)
echo "Total runs to download: $RUN_COUNT"

# Step 4: Download FASTQ files
echo "Step 4: Downloading FASTQ files..."
cd "$OUTPUT_DIR"

counter=1
while IFS= read -r run_accession; do
    echo "[$counter/$RUN_COUNT] Downloading $run_accession..."
    
    # Download using fasterq-dump
    fasterq-dump "$run_accession" --outdir . --threads 4 --progress
    
    # Compress FASTQ files to save space
    if [ -f "${run_accession}_1.fastq" ]; then
        echo "  Compressing ${run_accession}_1.fastq..."
        gzip "${run_accession}_1.fastq"
    fi
    
    if [ -f "${run_accession}_2.fastq" ]; then
        echo "  Compressing ${run_accession}_2.fastq..."
        gzip "${run_accession}_2.fastq"
    fi
    
    # Handle single-end reads
    if [ -f "${run_accession}.fastq" ]; then
        echo "  Compressing ${run_accession}.fastq..."
        gzip "${run_accession}.fastq"
    fi
    
    ((counter++))
done < "../$TEMP_DIR/run_accessions.txt"

cd ..

# Step 5: Copy metadata and supplementary files
echo "Step 5: Organizing metadata and supplementary files..."
if [ -f "$TEMP_DIR/GSE203612_series_matrix.txt" ]; then
    cp "$TEMP_DIR/GSE203612_series_matrix.txt" "$OUTPUT_DIR/"
fi

# Copy any extracted supplementary files
if [ -d "$TEMP_DIR" ]; then
    find "$TEMP_DIR" -name "*.gz" -o -name "*.txt" -o -name "*.csv" | while read file; do
        if [[ "$(basename "$file")" != "run_accessions.txt" && "$(basename "$file")" != "GSE203612_series_matrix.txt" ]]; then
            cp "$file" "$OUTPUT_DIR/"
        fi
    done
fi

# Step 6: Generate metadata summary
echo "Step 6: Generating download summary..."
echo "GSE203612 Download Summary" > "$OUTPUT_DIR/download_summary.txt"
echo "Pan-cancer single-cell RNA-seq with endometrial cancer samples" >> "$OUTPUT_DIR/download_summary.txt"
echo "Date: $(date)" >> "$OUTPUT_DIR/download_summary.txt"
echo "BioProject: $PROJECT_ID" >> "$OUTPUT_DIR/download_summary.txt"
echo "GEO Series: $GSE_ID" >> "$OUTPUT_DIR/download_summary.txt"
echo "Total runs downloaded: $RUN_COUNT" >> "$OUTPUT_DIR/download_summary.txt"
echo "" >> "$OUTPUT_DIR/download_summary.txt"
echo "Study details:" >> "$OUTPUT_DIR/download_summary.txt"
echo "- Pan-cancer analysis across 15 cancer types" >> "$OUTPUT_DIR/download_summary.txt"
echo "- Includes endometrial cancer samples" >> "$OUTPUT_DIR/download_summary.txt"
echo "- Single-cell RNA sequencing data" >> "$OUTPUT_DIR/download_summary.txt"
echo "- Published in Nature Genetics (2022)" >> "$OUTPUT_DIR/download_summary.txt"
echo "" >> "$OUTPUT_DIR/download_summary.txt"
echo "Files downloaded:" >> "$OUTPUT_DIR/download_summary.txt"
ls -lh "$OUTPUT_DIR"/*.fastq.gz >> "$OUTPUT_DIR/download_summary.txt" 2>/dev/null || echo "No FASTQ files found" >> "$OUTPUT_DIR/download_summary.txt"



# Step 8: Prepare for Cell Ranger processing
echo "Step 8: Preparing Cell Ranger directory structure..."
mkdir -p "$OUTPUT_DIR/cellranger_input"

# Create sample sheet for Cell Ranger
echo "Creating sample sheet for Cell Ranger processing..."
echo "Sample,R1,R2" > "$OUTPUT_DIR/cellranger_input/sample_sheet.csv"

# Group files by sample
for r1_file in "$OUTPUT_DIR"/*_1.fastq.gz; do
    if [ -f "$r1_file" ]; then
        base_name=$(basename "$r1_file" _1.fastq.gz)
        r2_file="$OUTPUT_DIR/${base_name}_2.fastq.gz"
        
        if [ -f "$r2_file" ]; then
            echo "$base_name,$r1_file,$r2_file" >> "$OUTPUT_DIR/cellranger_input/sample_sheet.csv"
        fi
    fi
done

# Create Cell Ranger processing script
cat > "$OUTPUT_DIR/run_cellranger.sh" << 'EOF'
#!/bin/bash

# Cell Ranger processing script for GSE203612
REFERENCE_GENOME="/path/to/refdata-gex-GRCh38-2020-A"  # Update this path
OUTPUT_DIR="./cellranger_output"
SAMPLE_SHEET="./cellranger_input/sample_sheet.csv"

mkdir -p "$OUTPUT_DIR"

# Process each sample
tail -n +2 "$SAMPLE_SHEET" | while IFS=',' read -r sample r1 r2; do
    echo "Processing sample: $sample"
    
    cellranger count \
        --id="$sample" \
        --transcriptome="$REFERENCE_GENOME" \
        --fastqs="$(dirname "$r1")" \
        --sample="$sample" \
        --localcores=8 \
        --localmem=32
        
    mv "$sample" "$OUTPUT_DIR/"
done

echo "Cell Ranger processing complete. Results in: $OUTPUT_DIR"
EOF

chmod +x "$OUTPUT_DIR/run_cellranger.sh"

# Cleanup
echo "Step 9: Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo ""
echo "=== Download Complete ==="
echo "Dataset: GSE203612 (Pan-cancer with endometrial samples)"
echo "Output directory: $OUTPUT_DIR"
echo "Files downloaded: $RUN_COUNT runs"
echo ""
echo "Next steps:"
echo "1. Run sample identification: cd $OUTPUT_DIR && python identify_endometrial_samples.py"
echo "2. Update reference genome path in: $OUTPUT_DIR/run_cellranger.sh"
echo "3. Run Cell Ranger processing: cd $OUTPUT_DIR && ./run_cellranger.sh"
echo "4. Review download summary: cat $OUTPUT_DIR/download_summary.txt"
echo ""
echo "Note: This is a pan-cancer study. Use the identification script to find endometrial samples."