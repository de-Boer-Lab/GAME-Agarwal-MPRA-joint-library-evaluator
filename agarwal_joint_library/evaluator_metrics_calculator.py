'''Calculate and save the final evaluation metrics.

This module handles the comparison between Predictor outputs (JSON) and Measured Data (Excel/CSV).
It calculates:
1. Pearson Correlation (per task)
2. Cell-Type Specificity (metrics across specific cell-type pairs)
'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

import os
import sys
import json
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from datetime import datetime, timezone

from config import EVALUATOR_NAME, EVALUATOR_DATA_DIR

def calculate_task_correlation(
    measured_df: pd.DataFrame,
    single_task_data: dict,
    measured_value_column: str,
    seq_id_column: str,
    chromosome_column: str = None,
    chromosomes_to_filter: list = None
    ):
    
    """
    Calculates Pearson r and extracts metadata for single prediction task,
    using pre-loaded measured_df and a single task data dictionary.
    
    Args:
        measured_df (pd.DataFrame): DataFrame containing the measured data.
        single_task_data (dict): Dictionary containing metadata and predictions for a single task.
        measured_value_column (str): The column name in measured_df that contains the measured values to correlate against.
        seq_id_column (str): The column name in measured_df that contains the sequence identifiers to merge on.
        chromosome_column (str, optional): The column name in measured_df that contains chromosome information (if filtering by chromosome is needed).
        chromosomes_to_filter (list, optional): List of chromosome names to filter by (if filtering by chromosome is needed).
        
    Returns:
        correlation_details (dict): Dictionary containing 'pearson_r' (float or None) and task metadata ("task_name", "task_type", "cell_type_actual")
    """
    
    # Extract metadata from single task data
    print(f"\n--- Extracting prediction_task metadata ---")
    task_name = single_task_data.get("name")
    task_type_actual = single_task_data.get("type_actual")
    cell_type_actual = single_task_data.get("cell_type_actual")
    predictions_dict =  single_task_data.get("predictions")
    scale_prediction_actual = single_task_data.get("scale_prediction_actual")
    scale_prediction_requested = single_task_data.get("scale_prediction_requested")
    
    pearson_r_value = None # If there's an error, default to None
    
    # --- Early exit checks ---
    # Error in predictions
    if "error" in predictions_dict:
        print("No predictions were returned for this task -> Skipping evaluation calculation")
        correlation_details = {
        'task_name': task_name, 
        'task_type': task_type_actual,
        'cell_type_actual': cell_type_actual,
        'pearson_r': pearson_r_value
        }
        return correlation_details

    # Scale mismatch
    if scale_prediction_actual != scale_prediction_requested:
        print("Predictions scale does not match requested scale (log): Skipping evaluation calculation")
        correlation_details = {
        'task_name': task_name, 
        'task_type': task_type_actual,
        'cell_type_actual': cell_type_actual,
        'pearson_r': pearson_r_value
        }
        return correlation_details

    # --- Data Validation and processing ---
    print("--- Validating data ---")
    if not isinstance(predictions_dict, dict):
        print(f"WARNING: 'predictions' in task: '{task_name}'\
            \nCell type: {cell_type_actual}\
            \nType: {task_type_actual}\
            \nis not a valid dictionary or is missing.")
    elif not predictions_dict:
        print(f"WARNING: 'predictions' dictionary is empty in task: '{task_name}'\
            \nCell type: {cell_type_actual}\
            \nType: {task_type_actual}")
    elif seq_id_column not in measured_df.columns:
        print(f"ERROR: Sequence ID column '{seq_id_column}' not found in measured_df.\
            \nCannot merge for task: {task_name}.")
    elif measured_value_column not in measured_df.columns:
         print(f"ERROR: Measured value column '{measured_value_column}'\
             \nfor cell type '{cell_type_actual}' not found in measured_df.\
             \nCannot correlate task: '{task_name}'.") 
         # NOTE: More checks can be added.
    else:
        # Proceed with calculation if checks pass
        # Create DataFrame from Predictions
        # The keys in 'predictions_dict' must match those in measured_df[seq_id_column]
        print("--- Creating predictions_df ---")
        predictions_df = pd.DataFrame(list(predictions_dict.items()), columns=[seq_id_column, 'Predicted_Value'])
        predictions_df['Predicted_Value'] = predictions_df['Predicted_Value'].apply(
            lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
            )

        # check if there are any NaN values in the prediction values before merging with measured_df
        # if there are NaNs in prediction values, skip the evaluation for this task
        na_rows = predictions_df[predictions_df['Predicted_Value'].isna()]
        if not na_rows.empty:
            print("NA values were found in the predictions, skipping evaluation")
            print(na_rows)
            correlation_details = {
            'task_name': task_name, 
            'task_type': task_type_actual,
            'cell_type_actual': cell_type_actual,
            'pearson_r': pearson_r_value
            }
            return correlation_details
        
        # Now select only the necessary columns of measured_df
        columns_to_keep = [seq_id_column, measured_value_column]
        if (chromosome_column and chromosomes_to_filter and (chromosome_column in measured_df)):
            if chromosome_column not in columns_to_keep:
                columns_to_keep.append(chromosome_column)
                
        measured_df_subset = measured_df[columns_to_keep].copy()
        # Use inner join to allow for subset testing.
        # This ensures we only evaluate on sequences present in BOTH the measured file and the predictions.
        # Validation that the Predictor returned the correct number of sequences happens in evaluator_RestAPI.py.
        merged_df = pd.merge(measured_df_subset, predictions_df, on=seq_id_column, how="inner")
        
        # Filter by chromosome (if needed)
        if (chromosomes_to_filter and chromosome_column and (chromosome_column in merged_df.columns)):
            print(f"Filtering chromosomes: {chromosomes_to_filter}...")
            merged_df[chromosome_column] = merged_df[chromosome_column].astype(str)
            # In order to handle None
            str_chromosomes_to_filter = [str(c) for c in chromosomes_to_filter]
            filtered_df = merged_df[merged_df[chromosome_column].isin(str_chromosomes_to_filter)]
        else:
            print("No chromosomes to filter.")
            filtered_df = merged_df # No chromosome filter
        
        # Drop any rows that have NaNs for either column, any NAs in the Predicted values should have already been caught
        print("Original size of the measurement file is:")
        print(filtered_df.shape)
        na_rows = filtered_df[filtered_df.isna().any(axis=1)]
        if not na_rows.empty:
            print("Rows with NaN values in any of measured value column (will be dropped):")
            print(na_rows)
    
        final_df = filtered_df.dropna() # Drop rows where measured data is NaN
        print("Size of data that will be used to calculate pearson r")
        print(final_df.shape)
        
        # Sanitize the final_df in case values are non-numeric
        if not final_df.empty:
            print("Sanitizing final_df in case values are non-numeric for correlation...")
            final_df.loc[:, 'Predicted_Value'] = pd.to_numeric(final_df['Predicted_Value'], errors='coerce')
            final_df.loc[:, measured_value_column] = pd.to_numeric(final_df[measured_value_column], errors='coerce')

            rows_before = len(final_df)
            final_df = final_df.dropna(subset=['Predicted_Value', measured_value_column])
            rows_dropped = rows_before - len(final_df)
            if rows_dropped > 0:
                print(f"Warning: Dropped {rows_dropped} rows with non-numeric values after coercion.")
                
            # Check for 0 variance in either column (all values identical)
            std_predicted = final_df['Predicted_Value'].std()
            std_measured = final_df[measured_value_column].std()
            
            # If predictions have 0 variance, there is no deviation, correlation is undefined -> score 0
            # If measured values have 0 variance, correlation is undefined -> score 0
            if std_predicted == 0 or std_measured == 0:
                print(f"WARNING: Zero variance detected for task '{task_name}'. Setting Pearson r to 0 to reflect no correlation.")
                pearson_r_value = 0.0
            
            else:
                # Calculate pearson r
                try:
                    r, _ = pearsonr(final_df['Predicted_Value'], final_df[measured_value_column])
                    pearson_r_value = 0.0 if np.isnan(r) else float(r)
                    print(f"Calculated Pearson r for {task_name}: {pearson_r_value}")
                except Exception as e:
                    print(f"Error during Pearson correlation calculation for task: '{task_name}': {e}")
        else:
            print(f"DataFrame is empty after numeric conversion and NaN drop for task: '{task_name}'")
            
    correlation_details = {
        'task_name': task_name, 
        'task_type': task_type_actual,
        'cell_type_actual': cell_type_actual,
        'pearson_r': pearson_r_value
    }
    return correlation_details

def calculate_cell_type_specificity(
    predictions_file_content: dict,
    measured_df: pd.DataFrame,
    seq_id_column: str,
    evaluator_name: str,
    predictor_name: str,
    cell_type_pairs: list,
    measured_col_map: dict
    ):
    
    """
    Calculates cell-type specificity by correlating the difference in expression
    between pairs of cell types.
    
    Handles cases where predictions may have identical values (0 variance) by assigning
    a Pearson r of 0 in such scenarios.
    
    Args:
        predictions_file_content (dict): The loaded JSON content from the predictor.
        measured_df (pd.DataFrame): DataFrame containing the measured data.
        seq_id_column (str): The column name in measured_df that contains the sequence identifiers to merge on.
        evaluator_name (str): The name of the evaluator.
        predictor_name (str): The name of the predictor.
        cell_type_pairs (list): A list of tuples, each containing cell types to be compared.  
        measured_col_map (dict): A dictionary mapping cell type names to their corresponding
                                 column names in the measured_df.
    
    Returns:
        evaluation_output (pd.DataFrame): DataFrame containing the cell-type specificity results.
    """
    
    print("\n----- Calculating Cell-Type Specificity -----")
    
    predictions_data = {} # Consolidate predictions into a dictionary first
    
    for task in predictions_file_content.get("prediction_tasks", []):
        req_cell_type = task.get("cell_type_requested")
        predictions = task.get("predictions")
        
        scale_prediction_actual = task.get("scale_prediction_actual")
        scale_prediction_requested = task.get("scale_prediction_requested")
        
        # Skip if scale mismatch
        if scale_prediction_actual != scale_prediction_requested:
            print("Predictions scale does not match requested scale (log): Skipping calculation")
            return None
        
        # Skip if predictions are missing or contain an error
        if not predictions or "error" in predictions:
            print(f"Skipping cell type '{req_cell_type}' due to missing predictions or error.")
            continue
        
        # Convert predictions to DataFrame
        df = pd.DataFrame(list(predictions.items()), columns=[seq_id_column, f'prediction_{req_cell_type}'])
        # Handle cases where predictions are lists vs scalars
        df[f'prediction_{req_cell_type}'] = df[f'prediction_{req_cell_type}'].apply(
            lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
        )
        
        # Set index to 'name' for easier merging later
        predictions_data[f'prediction_{req_cell_type}'] = df.set_index(seq_id_column)[f'prediction_{req_cell_type}']
        
    if not predictions_data:
        print("No valid predictions found for any cell type. Skipping specificity calculation.")
        return None
    
    # Merge all predictions into a single DataFrame
    predictions_df = pd.DataFrame(predictions_data) 
    
    # Check for NA values in predictions, if any exist, skip specificity calculation
    if predictions_df.isna().any().any():  # any column, any row
        print("NA values were found in the predictions, skipping specificity evaluation")
        return None
    
    # Merge measured data with predictions at ONCE
    predictions_df = predictions_df.reset_index().rename(columns={'index': seq_id_column}) # Reset index to merge
    # Use inner join to allow for subset testing.
    # This ensures we only evaluate on sequences present in BOTH the measured file and the predictions.
    # Validation that the Predictor returned the correct number of sequences happens in evaluator_RestAPI.py.
    merged_df = pd.merge(measured_df, predictions_df, left_on=seq_id_column, right_on=seq_id_column, how="inner") # Use inner join to ensure alignment
    
    if merged_df.empty:
        print("Merged DataFrame is empty after merging measured and predicted data. Skipping specificity calculation.")
        return None
    
    final_df = merged_df.dropna().copy() # This is to drop rows with any NaNs in measured or predicted values
    # Primarily if the measured data has more sequences than predictions
    
    evaluation_output = []
    
    prediction_task_data_no_predictions = {}
    for task in predictions_file_content["prediction_tasks"]:
        cell_type = task["cell_type_requested"]
        # Store the metadata for this cell type without the predictions to include in the output
        prediction_task_data_no_predictions[cell_type] = [
            {k: v for k, v in task.items() if k != "predictions"}
        ]
    
    # Iterate through each cell type pair to calculate specificity
    for cell1, cell2 in cell_type_pairs:
        # Define the expected column names
        pred_col1 = f'prediction_{cell1}'
        pred_col2 = f'prediction_{cell2}'

        # Check for the actual column names
        if pred_col1 not in final_df.columns or pred_col2 not in final_df.columns:
            print(f"One or both of the prediction columns '{pred_col1}' or '{pred_col2}' not found. Skipping this pair.")
            continue
        
        measured_col1 = measured_col_map.get(cell1)
        measured_col2 = measured_col_map.get(cell2)
        
        # Calculate differences
        final_df[f'{cell1}-{cell2}_measured'] = final_df[measured_col1] - final_df[measured_col2]
        final_df[f'{cell1}-{cell2}_predicted'] = final_df[f'prediction_{cell1}'] - final_df[f'prediction_{cell2}']
        
        # Check for 0 variance in either column
        std_measured_diff = final_df[f'{cell1}-{cell2}_measured'].std()
        std_predicted_diff = final_df[f'{cell1}-{cell2}_predicted'].std()
        
        if std_measured_diff == 0 or std_predicted_diff == 0:
            print(f"WARNING: Zero variance detected for cell type pair '{cell1}-{cell2}'. Setting Pearson r to 0 to reflect no correlation.")
            pearson_r_val = 0.0
        else:
            try:
                r, _ = pearsonr(final_df[f'{cell1}-{cell2}_measured'], final_df[f'{cell1}-{cell2}_predicted'])
                pearson_r_val = 0.0 if np.isnan(r) else float(r)
            except Exception as e:
                print(f"Error calculating Pearson r for cell type pair '{cell1}-{cell2}': {e}")
                pearson_r_val = None
                
        # Always append, converting None to "NaN"
        val_str = "NaN" if pearson_r_val is None else str(pearson_r_val)
        
        # Get UTC timestamp for this calculation to match correlation schema
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
        
        new_result = {
            'evaluator_name': evaluator_name,
            'description': f'Cell type specific expression ({cell1} - {cell2})',
            'predictor_name': predictor_name,
            'time_stamp': timestamp,
            'metric': 'pearson_r',
            'value': val_str,
            'prediction_task(s)_data': str(prediction_task_data_no_predictions[cell1]) + " - " + str(prediction_task_data_no_predictions[cell2])
        }
        evaluation_output.append(new_result)
            
    return pd.DataFrame(evaluation_output)

def calculate_and_save_metrics(saved_predictions_path, output_dir):
    """
    Calculates custom evaluation metrics and saves them to CSV files.
    This is the primary function to customize for a new evaluator.
    """
    print("----- Starting Evaluation Calculation and Saving as CSV -----")
    
    
    MEASURED_DATA_PATH = os.path.join(EVALUATOR_DATA_DIR, "56k_measured_file.xlsx") # NOTE: This may not be the same for other evaluators
    print(f"Using measured data from: {MEASURED_DATA_PATH}")
    print(f"Using predictions from: {saved_predictions_path}")
    print(f"Correlation metadata will be saved in {output_dir}")
 
    # Define output paths
    evaluation_metrics_filepath = os.path.join(output_dir, f"evaluation_summary_{EVALUATOR_NAME}.csv") # NOTE: This is now evaluation summary, not correlation_summary
    
    # Initialize an empty list to get summary for all tasks
    all_task_correlation_results = []
    # Initialize cell-type specificity result container
    specificity_results_df = None
    
    try:
        try:
            # Load measured data file and predictions file ONCE (not with every function call).
            # NOTE: Evaluator builders: If measured_file_path is not an excel file,
            # this line (pd.read_excel) will need to be adjusted or replaced with the
            # appropriate pandas read function (e.g. pd.read_excel, pd.read_csv with different sep)
            # or custom loading logic (e.g. for .npy files).
            measured_df = pd.read_excel(MEASURED_DATA_PATH, header=0)
            
            # Now load predictions
            with open(saved_predictions_path, 'r') as f:
                predictions_file_content = json.load(f)
            
            # Extract Predictor Name
            predictor_name_base = predictions_file_content.get("predictor_name", "UnknownPredictor") # Resort UnknownPredictor if predictor name is not available
            predictor_name = predictor_name_base.replace(" ", "_").replace("/", "_") # Sanitize spaces and slashes for filenames, add underscores

        except Exception as e:
            print(f"Error loading data files: {e}")
            sys.exit(1) # Exit if essential data can't be loaded
            
        try:
            seq_column = "name" # This can change depending on data
            measured_value_columns_map = {
                "lymphoblasts (K562)": "K562 [log2(rna/dna)]",
                "human hepatocytes (HepG2)": "HepG2 [log2(rna/dna)]", 
                "induced pluripotent stem cells (iPS cells; WTC11)": "WTC11 [log2(rna/dna)]"
            }
            
            #These will be used to calculate cell type specificity predictions
            cell_type_pairs_for_specificity = [
                ("human hepatocytes (HepG2)", "induced pluripotent stem cells (iPS cells; WTC11)"),
                ("human hepatocytes (HepG2)", "lymphoblasts (K562)"),
                ("induced pluripotent stem cells (iPS cells; WTC11)", "lymphoblasts (K562)")
            ]
            
            if (
                "prediction_tasks" not in predictions_file_content or
                # Also flag cases in case prediction_tasks key is returned empty
                not predictions_file_content["prediction_tasks"] or
                # And flag if any 'predictions' keys are empty
                any(not key.get("predictions") for key in predictions_file_content["prediction_tasks"])
            ):
                print("WARNING: 'prediction_tasks' key missing, empty, or one of the tasks has empty predictions.")
            else:
                # Loop through each prediction_task from Predictor
                # Calculate the correlation for each task separately
                for task_index, single_task_data_dict in enumerate(predictions_file_content["prediction_tasks"]):
                    if not isinstance(single_task_data_dict, dict):
                        print(f"WARNING: Task item at index {task_index} is not a dictionary. Skipping!")
                        continue
                    
                    # Extract metadata from this task
                    task_type_actual = single_task_data_dict.get("type_actual")
                    predicted_cell_type = single_task_data_dict.get("cell_type_actual")
                    # We also want to extract the cell_type_requested to map it to measured_value_columns_map
                    requested_cell_type = single_task_data_dict.get("cell_type_requested")

                    # Find the corresponding measured data column from the map
                    measured_col_for_task = measured_value_columns_map.get(requested_cell_type)
                    
                    print(f"\nProcessing task {task_index+1} (Cell Type: {predicted_cell_type}). Correlating against measured column '{measured_col_for_task}'\
                    for requested cell type '{requested_cell_type}'.")
                    prediction_task_data_no_predictions = [{k: v for k, v in single_task_data_dict.items() if k != "predictions"}]

                    # Call the correlation calculation function
                    task_correlation_dict = calculate_task_correlation(
                        measured_df=measured_df,
                        single_task_data=single_task_data_dict,
                        measured_value_column=measured_col_for_task,
                        seq_id_column=seq_column # chromosome_column and chromosomes_to_filter_list can be added as arguments
                    )
                    
                    if task_correlation_dict:
                        pearson_r_value = task_correlation_dict.get('pearson_r')
                        
                        # Convert None (failed/disqualified tasks) to "NaN" for the CSV
                        val_str = "NaN" if pearson_r_value is None else str(pearson_r_value)
                        
                        # Get UTC timestamp for this calculation
                        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
                        
                        description = f"Agarwal Joint MPRA ({requested_cell_type})"
                        all_task_correlation_results.append({
                            "evaluator_name": EVALUATOR_NAME,
                            "description": description,
                            "predictor_name": predictor_name,
                            "time_stamp": timestamp,
                            'metric': 'pearson_r',
                            'value': val_str,
                            'prediction_task(s)_data': prediction_task_data_no_predictions
                        })
            
        except Exception as e:
            print(f"An error occurred during correlation calculation: {e}")
            
        try:
            # Calculate cell-type specificity
            specificity_results_df = calculate_cell_type_specificity(
                predictions_file_content=predictions_file_content,
                measured_df=measured_df,
                seq_id_column=seq_column,
                evaluator_name=EVALUATOR_NAME,
                predictor_name=predictor_name,
                cell_type_pairs=cell_type_pairs_for_specificity,
                measured_col_map=measured_value_columns_map
            )
        except Exception as e:
            print(f"An error occurred during cell-type specificity calculation: {e}")
            
        # Once all the data is received, save them all into a summary CSV
        # print(all_task_correlation_results)
        all_results = []
        
        # Add correlation results
        if all_task_correlation_results:
            all_results.extend(all_task_correlation_results)
            
        # Add specificity results
        if specificity_results_df is not None and not specificity_results_df.empty:
            # Convert the DataFrame back to a list of dicts to easily combine them
            all_results.extend(specificity_results_df.to_dict('records'))

        # Save to a single CSV
        if all_results:
            summary_df = pd.DataFrame(all_results)
            csv_file_exists: bool = os.path.isfile(evaluation_metrics_filepath)
            
            try:
                summary_df.to_csv(evaluation_metrics_filepath, mode='a',
                                  sep='\t', header=(not csv_file_exists), index=False)
                if csv_file_exists:
                    print(f"Appended merged metrics to existing CSV file: {evaluation_metrics_filepath}")
                else:
                    print(f"Created a new merged metrics CSV file: {evaluation_metrics_filepath}")
            except IOError as e:
                print(f"\nFATAL: No metrics results were saved! {e}", file=sys.stderr)
        else:
            print("\nNo valid metrics were calculated to save.")


    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        