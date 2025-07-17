#!/bin/bash

# Download script for SRP349751 - Endometrial Carcinoma scRNA-seq Dataset
# Dataset: Single-cell transcriptomic analysis of human endometrioid endometrial carcinoma
# Publication: Nature Communications (2022)
# Samples: 15 experiments (6 EEC, 5 AEH, 4 Normal endometrial)
# Platform: Illumina NovaSeq 6000, Paired-end

set -e  # Exit on any error

# Configuration
BASE_DIR=/gpfs/scratch/nk4167/
cd "$BASE_DIR"
PROJECT_ID="SRP349751"
OUTPUT_DIR="./SRP349751_data"
TEMP_DIR="./temp_${PROJECT_ID}"

# Create directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_DIR"

echo "=== Downloading SRP349751 Endometrial Carcinoma Dataset ==="
echo "Output directory: $OUTPUT_DIR"
echo "Temporary directory: $TEMP_DIR"

# Step 1: Get run information using SRA toolkit
echo "Step 1: Fetching run information for $PROJECT_ID..."
esearch -db sra -query "$PROJECT_ID" | efetch -format runinfo > "$TEMP_DIR/runinfo.csv"

# Check if we got the run info
if [ ! -s "$TEMP_DIR/runinfo.csv" ]; then
    echo "Error: Failed to fetch run information for $PROJECT_ID"
    exit 1
fi

echo "Found the following runs:"
cat "$TEMP_DIR/runinfo.csv" | cut -d',' -f1,12,30 | head -10

# Step 2: Extract run accessions
echo "Step 2: Extracting run accessions..."
tail -n +2 "$TEMP_DIR/runinfo.csv" | cut -d',' -f1 > "$TEMP_DIR/run_accessions.txt"

RUN_COUNT=$(cat "$TEMP_DIR/run_accessions.txt" | wc -l)
echo "Total runs to download: $RUN_COUNT"

# Step 3: Download FASTQ files
echo "Step 3: Downloading FASTQ files..."
cd "$OUTPUT_DIR"

counter=1
while IFS= read -r run_accession; do
    echo "[$counter/$RUN_COUNT] Downloading $run_accession..."
    
    # Download using fasterq-dump (faster than fastq-dump)
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

# Step 4: Generate metadata summary
echo "Step 4: Generating metadata summary..."
echo "SRP349751 Download Summary" > "$OUTPUT_DIR/download_summary.txt"
echo "Date: $(date)" >> "$OUTPUT_DIR/download_summary.txt"
echo "Total runs downloaded: $RUN_COUNT" >> "$OUTPUT_DIR/download_summary.txt"
echo "" >> "$OUTPUT_DIR/download_summary.txt"
echo "Files downloaded:" >> "$OUTPUT_DIR/download_summary.txt"
ls -lh "$OUTPUT_DIR"/*.fastq.gz >> "$OUTPUT_DIR/download_summary.txt" 2>/dev/null || echo "No FASTQ files found" >> "$OUTPUT_DIR/download_summary.txt"

# Copy run info for reference
cp "$TEMP_DIR/runinfo.csv" "$OUTPUT_DIR/"

# Step 5: Prepare for Cell Ranger processing
echo "Step 5: Preparing Cell Ranger directory structure..."
mkdir -p "$OUTPUT_DIR/cellranger_input"

# Create sample sheet for Cell Ranger
echo "Creating sample sheet for Cell Ranger processing..."
echo "Sample,R1,R2" > "$OUTPUT_DIR/cellranger_input/sample_sheet.csv"

# Group files by sample (assuming standard 10X naming)
for r1_file in "$OUTPUT_DIR"/*_1.fastq.gz; do
    if [ -f "$r1_file" ]; then
        base_name=$(basename "$r1_file" _1.fastq.gz)
        r2_file="$OUTPUT_DIR/${base_name}_2.fastq.gz"
        
        if [ -f "$r2_file" ]; then
            echo "$base_name,$r1_file,$r2_file" >> "$OUTPUT_DIR/cellranger_input/sample_sheet.csv"
        fi
    fi
done

# Step 6: Create Cell Ranger processing script template
cat > "$OUTPUT_DIR/run_cellranger.sh" << 'EOF'
#!/bin/bash

# Cell Ranger processing script for SRP349751
# Modify the reference genome path and other parameters as needed

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
        
    # Move output to organized directory
    mv "$sample" "$OUTPUT_DIR/"
done

echo "Cell Ranger processing complete. Results in: $OUTPUT_DIR"
EOF

chmod +x "$OUTPUT_DIR/run_cellranger.sh"

# Cleanup
echo "Step 7: Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo ""
echo "=== Download Complete ==="
echo "Dataset: SRP349751 (Endometrial Carcinoma)"
echo "Output directory: $OUTPUT_DIR"
echo "Files downloaded: $RUN_COUNT runs"
echo ""
echo "Next steps:"
echo "1. Update the reference genome path in: $OUTPUT_DIR/run_cellranger.sh"
echo "2. Run Cell Ranger processing: cd $OUTPUT_DIR && ./run_cellranger.sh"
echo "3. Review download summary: cat $OUTPUT_DIR/download_summary.txt"
echo ""
echo "Dataset ready for analysis!"