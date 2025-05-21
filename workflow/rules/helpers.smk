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