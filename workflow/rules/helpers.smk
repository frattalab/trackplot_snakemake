from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import os
import sys

# Simple store for view and focus coordinates of each input event
@dataclass
class CoordinatePair:
    view: str = None
    focus: str = None

# Parse BED files matched by Name field into interval strings ready for trackplot
def parse_bed_files(view_bed_path: str, focus_bed_path: str, max_width: int) -> Dict[str, CoordinatePair]:
    '''
    Parse BED files and create coordinate pairs, adjusting intervals that exceed max_width.
    
    Parameters:
    -----------
    view_bed_path : str
        Path to the BED file containing view coordinates (full plotting window)
    focus_bed_path : str
        Path to the BED file containing focus coordinates (window to highlight)
    max_width : int
        Maximum width allowed for intervals. If an interval exceeds this width,
        it will be adjusted to center around the focus region with the max_width.
    
    Returns:
    --------
    Dict[str, CoordinatePair]
        Dictionary of coordinate pairs keyed by name
    '''
    # Dictionary to store the coordinate pairs, keyed by the name field
    coord_dict = {}
    
    # Parse the view BED file
    with open(view_bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue  # Skip comments and empty lines
                
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue  # Skip lines without enough fields
                
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            strand = "+" if len(fields) < 6 else fields[5]
            
            # Create the view string in the format chr:start-end:strand
            view_str = f"{chrom}:{start}-{end}:{strand}"
            
            # Initialize the dictionary entry with the view string
            coord_dict[name] = CoordinatePair(view=view_str)
    
    # Parse the focus BED file
    with open(focus_bed_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue  # Skip comments and empty lines
                
            fields = line.strip().split('\t')
            if len(fields) < 4:
                continue  # Skip lines without enough fields
                
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            
            # Create the focus string in the format start-end
            focus_str = f"{start}-{end}"
            
            # Update the dictionary entry if the name exists
            if name in coord_dict:
                coord_dict[name].focus = focus_str
                
                # Extract the original view data
                view_parts = coord_dict[name].view.split(':')
                view_chrom = view_parts[0]
                view_coords = view_parts[1].split('-')
                view_start = int(view_coords[0])
                view_end = int(view_coords[1])
                view_strand = view_parts[2]
                
                # Calculate current width
                width = view_end - view_start
                
                # If width exceeds max_width, adjust it around the focus
                if width > max_width:
                    # Output warning to stderr
                    print(f"Warning: Interval '{name}' exceeds max_width of {max_width}bp (current width: {width}bp). Adjusting to center around focus region.", file=sys.stderr)
                    
                    # Calculate center of focus region
                    focus_center = (start + end) // 2
                    
                    # Set new view start and end centered on focus
                    new_start = max(0, focus_center - max_width // 2)
                    new_end = new_start + max_width
                    
                    # Update the view string
                    new_view_str = f"{view_chrom}:{new_start}-{new_end}:{view_strand}"
                    
                    # Add information about the change to the warning
                    print(f"  - Original view: {coord_dict[name].view}", file=sys.stderr)
                    print(f"  - Updated view: {new_view_str}", file=sys.stderr)
                    
                    coord_dict[name].view = new_view_str
    
    return coord_dict


def is_bgzipped_and_indexed(file_path: str) -> Tuple[bool, bool]:
    """
    Check if a file is bgzipped and has a tabix index.
    
    Parameters:
    -----------
    file_path : str
        Path to the file to be checked
        
    Returns:
    --------
    Tuple[bool, bool]
        A tuple of (is_bgzipped, is_indexed) where:
        - is_bgzipped: True if the file exists and has a .gz or .bgz extension
        - is_indexed: True if a tabix index (.tbi) file exists for this file
    """
    # Check if the file exists first
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"file_path doesn't exist - {file_path}")
    
    # Check if the file has a .gz or .bgz extension
    is_bgzipped = file_path.endswith('.gz') or file_path.endswith('.bgz')
    
    # Check if a tabix index (.tbi) file exists for this file
    is_indexed = False
    if is_bgzipped:
        index_path = file_path + '.tbi'
        is_indexed = os.path.isfile(index_path)
    
    return is_bgzipped, is_indexed


# def target_gtf_file(gtf_path: str, output_directory: str) -> str:
#     """
#     Determine the appropriate GTF file to use for trackplot.
    
#     Parameters:
#     -----------
#     gtf_path : str
#         Path to the input GTF file
#     output_directory : str
#         Directory where a processed (bgzipped) GTF file should be created if needed
        
#     Returns:
#     --------
#     str
#         Path to the GTF file to use (either input file if already bgzipped
#         or path to where bgzipped file will be created)
#     """
#     # Check if the input file exists
#     if not os.path.isfile(gtf_path):
#         raise FileNotFoundError(f"GTF file not found: {gtf_path}")
    
#     # Check if the file is already bgzipped
#     is_bgzipped, _ = is_bgzipped_and_indexed(gtf_path)
    
#     # If the file is already properly bgzipped, return it
#     if is_bgzipped:
#         return gtf_path
    
#     # File needs bgzipping, determine the output path
#     input_filename = os.path.basename(gtf_path)
#     output_filename = input_filename + '.bgz'
    
#     # Create the full output path
#     os.makedirs(output_directory, exist_ok=True)
#     output_path = os.path.join(output_directory, output_filename)
    
#     return output_path


# def target_tabix_file(file_path: str, output_directory: str) -> str:
#     """
#     Determine the path to the tabix index file for a given input file.
#     If the input file is not already bgzipped, determines the path to
#     where the bgzipped version would be created first.
    
#     Parameters:
#     -----------
#     file_path : str
#         Path to the input file
#     output_directory : str
#         Directory where processed files should be created if needed
        
#     Returns:
#     --------
#     str
#         Path to the tabix index file (either existing .tbi or where it should be created)
#     """
#     # Check if the input file exists
#     if not os.path.isfile(file_path):
#         raise FileNotFoundError(f"File not found: {file_path}")
    
#     # Check if the file is already bgzipped and indexed
#     is_bgzipped, is_indexed = is_bgzipped_and_indexed(file_path)
    
#     # If the file is already properly bgzipped and indexed, return the index path
#     if is_bgzipped and is_indexed:
#         return file_path + '.tbi'
    
#     # If the file is bgzipped but not indexed, we just need to create the index
#     if is_bgzipped:
#         return file_path + '.tbi'
    
#     # File needs bgzipping first, then indexing
#     input_filename = os.path.basename(file_path)
#     output_filename = input_filename + '.bgz'
    
#     # Create the full output path for the bgzipped file
#     os.makedirs(output_directory, exist_ok=True)
#     output_path = os.path.join(output_directory, output_filename)
    
#     # Return the path to the tabix index that will be created
#     return output_path + '.tbi'


# def check_interval_files(interval_list_path: str, output_directory: str) -> Tuple[Dict[str, str], str, List[str]]:
#     """
#     Checks the interval files listed in the interval list file to identify which ones 
#     need bgzip compression and tabix indexing.
    
#     Parameters:
#     -----------
#     interval_list_path : str
#         Path to the interval list file
#     output_directory : str
#         Directory where processed files will be stored
        
#     Returns:
#     --------
#     Tuple[Dict[str, str], str, List[str]]
#         (
#             interval_files_to_process: Dict mapping original file paths to target paths,
#             target_interval_list_path: Path to the updated interval list file,
#             updated_intervals: List of lines for the updated interval file
#         )
        
#     Raises:
#     -------
#     FileNotFoundError
#         If the interval list file doesn't exist
#     """
#     if not os.path.exists(interval_list_path):
#         raise FileNotFoundError(f"Interval list file not found: {interval_list_path}")
    
#     # Ensure output directory exists
#     os.makedirs(output_directory, exist_ok=True)
    
#     interval_files_to_process: Dict[str, str] = {}
#     updated_intervals: List[str] = []

#     # Parse the interval list file
#     with open(interval_list_path, 'r') as f:
#         lines = f.readlines()
        
#     # First line should be the header
#     header = lines[0] if lines else "#file_location\tlabel_of_file\n"
#     updated_intervals.append(header)
    
#     # Process each file in the interval list
#     for line in lines[1:]:
#         if line.strip() and not line.startswith('#'):
#             fields = line.strip().split('\t')
#             if len(fields) >= 2:
#                 file_path = fields[0]
#                 label = fields[1]
                
#                 # Check if the file exists
#                 if not os.path.exists(file_path):
#                     raise FileNotFoundError(f"Interval file referenced in {interval_list_path} not found: {file_path}")
                
#                 # Check if the file is already bgzipped and indexed
#                 is_bgzipped, is_indexed = is_bgzipped_and_indexed(file_path)
                
#                 if is_bgzipped and is_indexed:
#                     # File is already processed correctly
#                     target_file = file_path
#                 else:
#                     # File needs processing
#                     if not is_bgzipped:
#                         # Need to bgzip
#                         target_file = target_file + ".bgz"
#                         if not (target_file.endswith('.gz') or target_file.endswith('.bgz')):
#                             # Add .bgz extension if neither .gz nor .bgz exists
#                             target_file = target_file + '.bgz'
#                     else:
#                         # Already bgzipped, but needs index
#                         target_file = file_path
                    
#                     # Add to the list of files to process
#                     interval_files_to_process[file_path] = target_file
                
#                 # Update the intervals file with the new path
#                 updated_line = f"{target_file}\t{label}\n"
#                 updated_intervals.append(updated_line)
    
#     # Create the target interval list path
#     base_name = os.path.basename(interval_list_path)
#     target_interval_list_path = os.path.join(output_directory, f"processed_{base_name}")
    
#     return interval_files_to_process, target_interval_list_path, updated_intervals

# def target_intervals_file(interval_list_path: str, output_directory: str) -> str:
#     """
#     Creates an updated interval list file with paths to processed (bgzipped and indexed) files.
    
#     Parameters:
#     -----------
#     interval_list_path : str
#         Path to the original interval list file
#     output_directory : str
#         Directory where processed files will be stored
        
#     Returns:
#     --------
#     str
#         Path to the updated interval list file
        
#     Raises:
#     -------
#     FileNotFoundError
#         If the interval list file doesn't exist
#     """
#     # Get the files to process and the updated intervals content
#     _, target_path, updated_intervals = check_interval_files(interval_list_path, output_directory)
    
#     # Write the updated interval list file
#     os.makedirs(os.path.dirname(target_path), exist_ok=True)
#     with open(target_path, 'w') as f:
#         f.writelines(updated_intervals)
    
#     return target_path

# def get_processed_interval_files(interval_list_path: str, output_directory: str) -> List[str]:
#     """
#     Returns a list of paths to the processed interval files that need to be created.
    
#     Parameters:
#     -----------
#     interval_list_path : str
#         Path to the interval list file
#     output_directory : str
#         Directory where processed files will be stored
        
#     Returns:
#     --------
#     List[str]
#         List of paths to the processed files (both bgzipped and indexed)
        
#     Raises:
#     -------
#     FileNotFoundError
#         If the interval list file doesn't exist
#     """
#     # Get the files that need processing
#     files_to_process, _, _ = check_interval_files(interval_list_path, output_directory)
    
#     # Generate a list of all output files (both bgzipped and indexed)
#     output_files: List[str] = []
#     for original_path, target_path in files_to_process.items():
#         output_files.append(target_path)  # Bgzipped file
#         output_files.append(target_path + ".tbi")  # Tabix index
    
#     return output_files



def test_is_bgzipped_and_indexed():
    """
    Test the is_bgzipped_and_indexed function with different cases.
    """
    # Test cases
    test_cases = [
        # Case 1: File exists but doesn't have .gz or .bgz extension
        {'file_path': 'data/intervals.txt', 'expected': (False, False), 'should_raise': False},
        
        # Case 2: File with .gz extension but no .tbi index
        {'file_path': 'data/example.gz', 'expected': (True, False), 'should_raise': False},
        
        # Case 3: File with .gz extension and .tbi index (positive case)
        {'file_path': 'data/example_with_index.gz', 'expected': (True, True), 'should_raise': False},
        
        # Case 4: Non-existent file (should raise FileNotFoundError)
        {'file_path': 'non_existent_file.gz', 'expected': None, 'should_raise': True}
        ]
    
    # Before running the test, create test files for cases 2, 3, and 4
    import os
    
    # Ensure data directory exists
    os.makedirs('data', exist_ok=True)
    
    # For case 2
    with open('data/intervals.txt', 'w') as f:
        f.write('test content')
    
    # For case 3
    with open('data/example.gz', 'w') as f:
        f.write('test content')
    
    # For case 4
    with open('data/example_with_index.gz', 'w') as f:
        f.write('test content')
    with open('data/example_with_index.gz.tbi', 'w') as f:
        f.write('test index content')
    
    # Run tests
    for i, case in enumerate(test_cases):
        file_path = case['file_path']
        expected = case['expected']
        should_raise = case['should_raise']
        
        try:
            result = is_bgzipped_and_indexed(file_path)
            if should_raise:
                assert False, f"Test case {i+1} failed: {file_path}, expected FileNotFoundError but got result {result}"
            else:
                assert result == expected, f"Test case {i+1} failed: {file_path}, expected {expected}, got {result}"
        except FileNotFoundError:
            if not should_raise:
                assert False, f"Test case {i+1} failed: {file_path}, unexpected FileNotFoundError"
            # Otherwise, this is the expected behavior for case 4
        
    # Clean up test files
    if os.path.exists('data/intervals.txt'):
        os.remove('data/intervals.txt')
    if os.path.exists('data/example.gz'):
        os.remove('data/example.gz')
    if os.path.exists('data/example_with_index.gz'):
        os.remove('data/example_with_index.gz')
    if os.path.exists('data/example_with_index.gz.tbi'):
        os.remove('data/example_with_index.gz.tbi')
    
    print("All tests passed!")

# Run the test
if __name__ == "__main__":
    test_is_bgzipped_and_indexed()