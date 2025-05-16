#!/usr/bin/env python3
"""
check_interval_inputs.py - A utility script to prepare input files for Trackplot pipeline.

This script checks if GTF and interval files are properly compressed with bgzip and 
indexed with tabix. If not, it processes them accordingly. This allows for trackplot to immediately parallelise across input intervals when running the snakemake pipeline

Usage:
    python check_inputs.py -g <gtf_file> -i <intervals_file> -o <output_path>

Requirements:
    - Python 3.9+
    - samtools/bgzip/tabix installed and in PATH
"""

import argparse
import os
import subprocess
import sys
import tempfile
from typing import List, Optional, Tuple, Dict, Any

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Check and prepare input files for Trackplot")
    parser.add_argument('-g', '--gtf', required=True, help='Path to the GTF file')
    parser.add_argument('-i', '--intervals', required=False, help='Path to the intervals.txt file')
    parser.add_argument('-o', '--output_path', default='intervals.validated.txt', help='Output path for re-processed/validated intervals.txt file')
    return parser.parse_args()

def is_bgzipped(file_path: str) -> bool:
    """Check if a file is bgzipped by checking file extension (.gz or .bgz)."""
    return file_path.endswith('.gz') or file_path.endswith('.bgz')

def is_tabix_indexed(file_path: str) -> bool:
    """Check if a tabix index exists for the given file."""
    index_path = file_path + '.tbi'
    return os.path.exists(index_path)

def check_file_status(file_path: str) -> Tuple[bool, bool]:
    """
    Check if a file is bgzipped and indexed, return status.
    
    Args:
        file_path: Path to the file to check
        
    Returns:
        Tuple of (is_compressed, is_indexed)
        
    Raises:
        FileNotFoundError: If the file does not exist
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")
    
    is_compressed = is_bgzipped(file_path)
    is_indexed = False
    
    if is_compressed:
        is_indexed = is_tabix_indexed(file_path)
    
    return is_compressed, is_indexed

def run_command(cmd: List[str], description: str = "") -> bool:
    """
    Run a shell command and handle errors.
    
    Args:
        cmd: Command to run as a list of strings
        description: Description of the command for error reporting
        
    Returns:
        True if command succeeded, False otherwise
    """
    try:
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error in {description}: {e}", file=sys.stderr)
        print(f"STDERR: {e.stderr.decode()}", file=sys.stderr)
        return False

def determine_tabix_preset(file_path: str) -> str:
    """
    Determine the appropriate tabix preset based on file extension.
    
    Args:
        file_path: Path to the file
        
    Returns:
        Tabix preset: 'bed' for bed files, 'gff' for others
    """
    if ".bed" in file_path.lower():
        return "bed"
    else:
        return "gff"


def process_file(file_path: str, tabix_preset: str) -> str:
    """
    Process a file: bgzip and tabix if needed.
    
    Args:
        file_path: Path to the file to process
        tabix_preset: Tabix preset to use ('bed' or 'gff')
        
    Returns:
        Path to the processed file
        
    Raises:
        RuntimeError: If processing fails at any step
    """
    is_compressed, is_indexed = check_file_status(file_path)
    
    if is_compressed and is_indexed:
        print(f"File {file_path} is already bgzipped and indexed.")
        return file_path
    
    # Determine output path
    output_path = file_path
    if not is_compressed:
        # For uncompressed files
        base_path = file_path
        output_path = base_path + '.bgz'

    # determine sort command based on file suffix
    if tabix_preset == "bed":
        sort_cmd_prefix = ["sort", "-k1V", "-k2n", "-k3n"]
    else:
        sort_cmd_prefix = ["sort", "-k1,1", "-k4,4n"]
    
    # Process file if needed
    if not is_compressed:
        # Sort and bgzip directly
        with tempfile.NamedTemporaryFile(suffix='.sorted') as temp_file:
            # Sort the file
            sort_cmd = sort_cmd_prefix + ["-o", temp_file.name, file_path]
            if not run_command(sort_cmd, f"sorting {tabix_preset} file"):
                raise RuntimeError(f"Failed to sort file: {file_path}")
            
            # Bgzip the sorted file
            bgzip_cmd = ["bgzip", "-c", temp_file.name]
            try:
                with open(output_path, 'wb') as out_file:
                    subprocess.run(bgzip_cmd, check=True, stdout=out_file)
                print(f"Created bgzipped file: {output_path}")
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"Failed to bgzip file {file_path}: {e}")
    
    # At this point, we have a compressed file (either the original or newly compressed)
    
    # First attempt to create tabix index directly
    if not is_indexed:
        tabix_cmd = ["tabix", "-p", tabix_preset, output_path]
        indexing_success = run_command(tabix_cmd, f"indexing {tabix_preset} file")
        
        # If indexing failed, we need to decompress and recompress
        if not indexing_success:
            print("Direct indexing failed. The file might not be properly bgzipped. Attempting to decompress and re-bgzip...")
            
            # Determine new output filename for the re-bgzipped file
            base_path = output_path
            if base_path.endswith('.gz'):
                base_path = base_path[:-3]  # Remove .gz extension
            elif base_path.endswith('.bgz'):
                base_path = base_path[:-4]  # Remove .bgz extension
            rebgzipped_output = base_path + '.bgz'
            
            # First decompress to a temporary file
            with tempfile.NamedTemporaryFile(suffix='.uncompressed', delete=False) as temp_uncompressed:
                temp_uncompressed_path = temp_uncompressed.name
            
            try:
                # Use zcat or gzip -dc to decompress
                decompress_cmd = ["zcat", output_path]
                try:
                    with open(temp_uncompressed_path, 'wb') as temp_file:
                        subprocess.run(decompress_cmd, check=True, stdout=temp_file)
                except subprocess.CalledProcessError as e:
                    os.unlink(temp_uncompressed_path)
                    raise RuntimeError(f"Failed to decompress file {output_path}: {e}")
                
                # Sort the decompressed file
                with tempfile.NamedTemporaryFile(suffix='.sorted', delete=False) as temp_sorted:
                    temp_sorted_path = temp_sorted.name
                
                sort_cmd = sort_cmd_prefix + ["-o", temp_sorted_path, temp_uncompressed_path]
                if not run_command(sort_cmd, f"sorting decompressed file"):
                    os.unlink(temp_uncompressed_path)
                    os.unlink(temp_sorted_path)
                    raise RuntimeError(f"Failed to sort decompressed file {temp_uncompressed_path}")
                
                # Clean up the uncompressed file, we don't need it anymore
                os.unlink(temp_uncompressed_path)
                
                # Bgzip the sorted file
                bgzip_cmd = ["bgzip", "-c", temp_sorted_path]
                try:
                    with open(rebgzipped_output, 'wb') as out_file:
                        subprocess.run(bgzip_cmd, check=True, stdout=out_file)
                    print(f"Created properly bgzipped file at: {rebgzipped_output}")
                except subprocess.CalledProcessError as e:
                    os.unlink(temp_sorted_path)
                    raise RuntimeError(f"Failed to bgzip sorted file {temp_sorted_path}: {e}")
                
                # Clean up the sorted file
                os.unlink(temp_sorted_path)
                
                # Try tabix on the re-bgzipped file
                tabix_cmd = ["tabix", "-p", tabix_preset, rebgzipped_output]
                if run_command(tabix_cmd, f"indexing re-bgzipped {tabix_preset} file"):
                    output_path = rebgzipped_output
                    print(f"Successfully indexed re-bgzipped file: {output_path}")
                else:
                    raise RuntimeError(f"Failed to create tabix index even after re-bgzipping: {rebgzipped_output}")
            except Exception as e:
                # Clean up any temporary files
                if os.path.exists(temp_uncompressed_path):
                    os.unlink(temp_uncompressed_path)
                if 'temp_sorted_path' in locals() and os.path.exists(temp_sorted_path):
                    os.unlink(temp_sorted_path)
                    
                # Propagate the error
                if isinstance(e, RuntimeError):
                    raise
                else:
                    raise RuntimeError(f"Error during fallback processing of {file_path}: {e}")
        else:
            print(f"Successfully indexed file: {output_path}")
    
    return output_path

def process_gtf(gtf_path: str) -> str:
    """
    Process GTF file: check if bgzipped and tabix indexed, process if needed.
    
    Args:
        gtf_path: Path to the GTF file
        
    Returns:
        Path to the processed GTF file
        
    Raises:
        RuntimeError: If processing fails
        FileNotFoundError: If the file doesn't exist
    """
    try:
        return process_file(gtf_path, "gff")
    except RuntimeError as e:
        print(f"Fatal error processing GTF file: {e}", file=sys.stderr)
        raise

def process_intervals_file(intervals_path: str, output_path: str) -> Optional[str]:
    """
    Process intervals file, update file paths if needed, and write a new intervals file.
    
    Args:
        intervals_path: Path to the intervals file
        output_dir: Output directory for processed files
        
    Returns:
        Path to the processed intervals file, or None if no intervals file was provided
    """
    if not intervals_path or not os.path.exists(intervals_path):
        print("No intervals file provided or file not found.")
        return None
    
    output_intervals = output_path
    
    # Read and process interval file
    updated_lines: List[str] = []
    updated = False
    
    with open(intervals_path, 'r') as infile:
        for line in infile:
            line = line.strip()
            if not line or line.startswith('#'):
                updated_lines.append(line)
                continue
            
            parts = line.split()
            if len(parts) < 2:
                updated_lines.append(line)  # Preserve malformed lines
                continue
            
            file_path = parts[0]
            try:
                is_compressed, is_indexed = check_file_status(file_path)
                
                # Process file if needed
                if not (is_compressed and is_indexed):
                    # Determine tabix preset based on file extension
                    tabix_preset = determine_tabix_preset(file_path)
                    
                    # Process file and get new path
                    proc_path = process_file(file_path, tabix_preset)
                    if proc_path != file_path:
                        updated = True
                        parts[0] = proc_path
            except FileNotFoundError:
                print(f"Warning: File specified in intervals not found: {file_path}", file=sys.stderr)
            
            updated_lines.append('\t'.join(parts))
    
    # Write updated intervals file
    if updated:
        with open(output_intervals, 'w') as outfile:
            for line in updated_lines:
                outfile.write(f"{line}\n")
        print(f"Updated intervals file saved to: {output_intervals}")
        return output_intervals
    else:
        print("No updates needed for intervals file.")
        return intervals_path

def main() -> None:
    """Main function."""
    args = parse_args()
    
    # Process GTF file
    try:
        processed_gtf = process_gtf(args.gtf)
        print(f"GTF file: {processed_gtf}")
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except RuntimeError as e:
        print(f"Error processing GTF file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Process intervals file if provided
    if args.intervals:
        try:
            processed_intervals = process_intervals_file(args.intervals, args.output_path)
            print(f"Intervals file: {processed_intervals}")
        except FileNotFoundError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        except RuntimeError as e:
            print(f"Error processing intervals file: {e}", file=sys.stderr)
            sys.exit(1)
    
    print("\nAll input files have now been validated or re-processed")

if __name__ == "__main__":
    main()