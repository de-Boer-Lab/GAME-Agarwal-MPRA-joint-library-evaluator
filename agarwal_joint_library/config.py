'''Configuration Script for Evaluator Name, Input File, and Preferred Data Format'''

import os

# --- Core Evaluator Settings ---
# Evaluator name for predictions file and metrics CSV
EVALUATOR_NAME = "agarwal_2025_joint_lib_56k"
# Name of the file being used for predictions
input_file = "56k_measured_final_sequence_file.xlsx"

# --- Directory Settings ---
# Get the absolute path of the script's directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Determine if running inside a container or not
if os.path.exists("/.singularity.d"):
    # Running inside the container
    EVALUATOR_DATA_DIR = "/evaluator_data/2023-03-03628-s5"
else:
    # Running outside the container
    EVALUATOR_DATA_DIR = os.path.join(SCRIPT_DIR, "evaluator_data", "2023-03-03628-s5")

EVALUATOR_INPUT_PATH = os.path.join(EVALUATOR_DATA_DIR, input_file)

output_filename_base = f'{EVALUATOR_NAME}_predictions_{input_file.replace(".xlsx", "")}'

# Debug logs for validation
print(f"Using input file: {EVALUATOR_INPUT_PATH}")

# --- API Communication Settings ---
REQUEST_FORMAT = "application/json"
REQUEST_FORMAT = REQUEST_FORMAT.lower()

RESPONSE_FORMAT = "application/msgpack"
RESPONSE_FORMAT = RESPONSE_FORMAT.lower()

# HTTP request retry
MAX_RETRIES = 50
RETRY_INTERVAL = 30 # Seconds to wait between each retry attempt