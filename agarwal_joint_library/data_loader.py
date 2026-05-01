'''Handle Loading and Validating Evaluator Input/Request Data'''

import os
import json
from collections import Counter
import functools
import pandas as pd
from config import EVALUATOR_INPUT_PATH, PLASMID_BACKBONE_INPUT_PATH

class DuplicateKeysError(ValueError):
    """Raised when duplicate keys are found in a JSON object."""
    pass

# Internal helper function to detect duplicates during JSON parsing
def _detect_duplicates(pairs, duplicate_keys_state):

    """
    Detects duplicate keys during JSON parsing and counts occurrences of each key.

    This function intercepts the key-value pairs provided by `json.loads` and ensures that
    duplicate keys are flagged. It constructs the dictionary normally but counts how often
    each key appears, recording any keys that occur more than once.

    Args:
        pairs (list of tuple): A list of key-value pairs at the current level of the JSON.
        duplicate_keys_state (dict): The dictionary to update with any duplicates found.

    Returns:
        result_dict: A dictionary created from the key-value pairs.
    """

    # Use a local Counter to count occurrences of keys at this level
    local_counts = Counter()
    result_dict = {}
    for key, value in pairs:
        # Increment the count for each key
        local_counts[key] += 1
        # If the key is a duplicate, record it in the duplicate_keys dictionary
        if local_counts[key] > 1:
            duplicate_keys_state[key] = local_counts[key]
        # Add the key-value pair to the resulting dictionary
        result_dict[key] = value
    return result_dict

def _process_results(data, duplicate_keys):
    """
    Checks the duplicate_keys dictionary and prints a report.

    Args:
        data (dict): The dictionary of parsed data. 
        duplicate_keys (dict): The dictionary of duplicates.
    
    Raises:
        DuplicateKeysError: If duplicate keys are found in the JSON structure.
        
    Returns:
        data: The parsed data if no errors or duplicates are found.
    """
    # Report duplicates if any were found
    if duplicate_keys:
        print("Duplicate keys found:")
        error_messages = [f"Key: '{key}', Count: {count}" for key, count in duplicate_keys.items()]
        raise DuplicateKeysError(f"Duplicate keys found:\n" + "\n".join(error_messages))
    else:
        print("No duplicates found.")
        return data # Return the parsed data if no duplicates.


# Function to check for duplicate keys in JSON object

def check_duplicates_from_string(json_string):

    """
    Parses a JSON string to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.loads()` to intercept key-value pairs
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts. Keys reused in separate objects within arrays (e.g. lists) 
    are not considered duplicates.

    Args:
        json_string (str): The JSON content as a string to parse and check for duplicates.

    Raises:
        json.JSONDecodeError: If the string is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict or list: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}
    
    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Parse the JSON string using the helper to track duplicates
    data = json.loads(json_string, object_pairs_hook=hook)
    
    return _process_results(data, duplicate_keys)
    
# Function for check for duplicate keys if input file is in JSON format

def check_duplicates_from_json(json_file_path):
    """
    Parses a JSON file to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.load()` to intercept key-value pairs 
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts and paths. Keys reused in separate objects within arrays 
    (e.g. lists) are not considered duplicates.
    
    Args:
        json_file_path (str): The path to the JSON file to parse and check for duplicates.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        json.JSONDecodeError: If the file content is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict or list: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}
    
    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Open and parse the JSON file, using the helper to track duplicates
    with open(json_file_path, 'r') as file:
        data = json.load(file, object_pairs_hook=hook)
        
    return _process_results(data, duplicate_keys)


def create_json_from_xlsx():
    """
    Loads and validates the input file specified in `config.py`.
    Parses an Excel file, extracts to Pandas DataFrame, and creates a JSON object for the Predictor.
    
    This function checks for file existence and output directory, handles JSON or MsgPack formats,
    and runs duplicate key validation.

    Returns:
        evaluator_dict (dict): The constructed JSON object containing request payload.

    Raises:
        FileNotFoundError: If the input file specified in `EVALUATOR_INPUT_PATH` 
                           does not exist.
        ValueError: If the data is malformed (e.g. invalid JSON/MsgPack),
                    or if duplicate keys are found during validation.

    """

    # Validate evaluator input file exists
    if not os.path.exists(EVALUATOR_INPUT_PATH):
        print(f"ERROR: Evaluator input file '{EVALUATOR_INPUT_PATH}' not found.")
        raise FileNotFoundError(f"Evaluator input file not found: {EVALUATOR_INPUT_PATH}")

    try:
        # Read the Excel file, treating the second row as the header (skipping the empty first row)
        df = pd.read_excel(EVALUATOR_INPUT_PATH, header=0)

        # Load the backbone
        backbone = pd.read_csv(PLASMID_BACKBONE_INPUT_PATH, header=0, sep= '\t', index_col=0)
        print(backbone)

        # get sequences from backbone
        upstream_seq = backbone.loc["upstream_padded", "sequence"]
        downstream_seq = backbone.loc["downstream", "sequence"]
        promoter_coordinates = json.loads(backbone.loc["promoter_coordinates", "sequence"])
        
        # print(f"DEBUG: promoter_coordinates type: {type(promoter_coordinates)}")
        # print(f"DEBUG: promoter_coordinates value: {promoter_coordinates}")
        
        # Extract the first column (seq_id (names)) and the last column (sequences -- "230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)")
        names = df.iloc[:, 0]      # sequence ID
        sequences = df.iloc[:, -1] # sequences
        
        # Create a dictionary, mapping each sequence name to its corresponding sequence
        sequence_dict = dict(zip(names, sequences))
        
        # Add predictions ranges that span the gene (69bp)
        prediction_ranges = {name: promoter_coordinates for name in sequence_dict.keys()}

        # Define the prediction tasks as a separate variable
        prediction_tasks_str = """
        [
            {
                "name": "agarwal_joint_lib_wtc11",
                "type": "expression",
                "cell_type": "induced pluripotent stem cells (iPS cells; WTC11)",
                "scale": "log",
                "species": "homo_sapiens"
            },
            {
                "name": "agarwal_joint_lib_k562",
                "type": "expression",
                "cell_type": "lymphoblasts (K562)",
                "scale": "log",
                "species": "homo_sapiens"
            },
            {
                "name": "agarwal_joint_lib_hepg2",
                "type": "expression",
                "cell_type": "human hepatocytes (HepG2)",
                "scale": "log",
                "species": "homo_sapiens"
            }
        ]
        """
        # print(f"prediction_tasks_str type: {type(prediction_tasks_str)}")
        
        # Check for duplicate keys in prediction_tasks
        prediction_tasks = check_duplicates_from_string(prediction_tasks_str)
        # print(f"prediction_tasks type: {type(prediction_tasks)}")
        
        # Build the JSON evaluator object
        evaluator_dict = {
            "readout": "point",
            "prediction_tasks": prediction_tasks,
            "upstream_seq": upstream_seq,
            "downstream_seq": downstream_seq,
            "sequences": sequence_dict,
            "prediction_ranges": prediction_ranges
        }
        
        # Convert the dictionary to a JSON string with indentation for readability
        json_string = json.dumps(evaluator_dict, indent=4)
        
        # Final validation check
        check_duplicates_from_string(json_string)
        print("Input data loaded and validated successfully.")
        
        return evaluator_dict
        
    except (json.JSONDecodeError, DuplicateKeysError, KeyError) as e:
        # Raise a general ValueError that the main script's handler
        # will catch and report cleanly
        raise ValueError(f"Input data is invalid.\nDetails: {e}") from e

if __name__ == '__main__':
    # ---- TEST SUITE (Remove this for other Evaluators) ----
    try:
        evaluator_dict = create_json_from_xlsx()
        
        # Basic structure checks
        expected_keys = {"readout", "prediction_tasks", "upstream_seq", "downstream_seq", "sequences", "prediction_ranges"}
        actual_keys = set(evaluator_dict.keys())
        missing = expected_keys - actual_keys
        extra = actual_keys - expected_keys
        
        if missing:
            print(f"FAIL: Missing keys: {missing}")
        if extra:
            print(f"INFO: Extra keys (not necessarily wrong): {extra}")
        if not missing:
            print("PASS: All expected keys present.")
        
        # Type checks
        assert isinstance(evaluator_dict["sequences"], dict), "FAIL: 'sequences' is not a dict"
        assert isinstance(evaluator_dict["prediction_tasks"], list), "FAIL: 'prediction_tasks' is not a list"
        assert isinstance(evaluator_dict["prediction_ranges"], dict), "FAIL: 'prediction_ranges' is not a dict"
        assert isinstance(evaluator_dict["upstream_seq"], str), "FAIL: 'upstream_seq' is not a string"
        assert isinstance(evaluator_dict["downstream_seq"], str), "FAIL: 'downstream_seq' is not a string"
        print("PASS: All type checks passed.")
        
        # Count checks
        num_seqs = len(evaluator_dict["sequences"])
        num_ranges = len(evaluator_dict["prediction_ranges"])
        num_tasks = len(evaluator_dict["prediction_tasks"])
        print(f"Sequences: {num_seqs}, Prediction ranges: {num_ranges}, Tasks: {num_tasks}")
        
        assert num_seqs == num_ranges, f"FAIL: sequences ({num_seqs}) != prediction_ranges ({num_ranges})"
        assert num_seqs > 0, "FAIL: No sequences loaded"
        assert num_tasks > 0, "FAIL: No prediction tasks loaded"
        print("PASS: Count checks passed.")
        
        # Sample output
        first_seq_id = next(iter(evaluator_dict["sequences"]))
        print(f"\nSample sequence ID: {first_seq_id}")
        print(f"Sample sequence (first 50 chars): {evaluator_dict['sequences'][first_seq_id][:50]}...")
        print(f"Sample prediction_range: {evaluator_dict['prediction_ranges'][first_seq_id]}")
        print(f"Upstream seq length: {len(evaluator_dict['upstream_seq'])}")
        print(f"Downstream seq length: {len(evaluator_dict['downstream_seq'])}")
        
        print("\n=== ALL CHECKS PASSED ===")
        
    except Exception as e:
        print(f"\n=== TEST FAILED ===\n{e}")
        import traceback
        traceback.print_exc()