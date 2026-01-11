'''Calculate and save the final evaluation metrics.'''

# NOTE: Every evaluator will do this slightly differently depending on how the data is presented

import os
import sys
import json
import pandas as pd
import numpy as np
import itertools
from scipy.stats import pearsonr
from datetime import datetime, timezone

from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH, EVALUATOR_DATA_DIR

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
    usign pre-loaded measured_df and a single task data dictionary.
    
    Args:
        measured_df (pd.DataFrame): The dataframe of measured values loaded once by the caller (Evaluator __main__ block).
        single_task_data (dict): The dictionary with predictions and metadata for a single task from the `prediction_tasks` 
                                 list of a predictions JSON.
        measured_value_column (str): Column name in measured_df for correlation.
        id_column (str): Common identifier for sequences column (IDs or sequences).
        chromosome_column (str, optional): Chromosome column name in measured_df. Defaults to 'None'.
        chromosomes_to_filter (list, optional): List of chromosomes for filtering. Defaults to 'None'. 

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
    if "error" in predictions_dict:
        print("No predictions were returned for this task -> Skipping evaluation calculation")
        correlation_details = {
        'task_name': task_name, 
        'task_type': task_type_actual,
        'cell_type_actual': cell_type_actual,
        'pearson_r': pearson_r_value
        }
        return correlation_details

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
        print("--- Creating predictions_df ---")
        predictions_df = pd.DataFrame(list(predictions_dict.items()), columns=[seq_id_column, 'Predicted_Value'])
        predictions_df['Predicted_Value'] = predictions_df['Predicted_Value'].apply(
            lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
            )

        #check here is there is NA is any of the prediction values
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
        merged_df = pd.merge(measured_df_subset, predictions_df, on=seq_id_column, how="left")
        
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
        print("Original size of the measurment file is:")
        print(filtered_df.shape)
        na_rows = filtered_df[filtered_df.isna().any(axis=1)]
        if not na_rows.empty:
            print("Rows with NaN values in any of measured value column (will be dropped):")
            print(na_rows)
    
        final_df = filtered_df.dropna()
        print("Size of data that will be used to calculate peasron r")
        print(final_df.shape)
        
        # Sanitize the final_df in case values are non-numeric
        if not final_df.empty:
            print("Sanitizing final_df in case values are non-numeric for correlation...")
            final_df.loc[:, 'Predicted_Value'] = pd.to_numeric(final_df['Predicted_Value'], errors='coerce')
            final_df.loc[:, measured_value_column] = pd.to_numeric(final_df[measured_value_column], errors='coerce')

            # Calculate pearson r
            try:
                r, _ = pearsonr(final_df['Predicted_Value'], final_df[measured_value_column])
                print(f"Calculated Pearson r for {task_name}: {r}") 
                if np.isnan(r):
                    print(f"WARNING: Pearson r is NaN for task '{task_name}'")
                    pearson_r_value = None
                else:
                    pearson_r_value = float(r)
            except ValueError as e:
                print(f"ValueError during Pearson correlation calculation for task: '{task_name}': {e}")
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
    
    Args:
        predictions_file_content (dict): The loaded JSON content from the predictor.
        measured_df (pd.DataFrame): DataFrame containing the measured data.
        evaluator_name (str): The name of the evaluator.
        predictor_name (str): The name of the predictor.
        cell_type_pairs (list): A list of tuples, each containing cell types to be compared.  
        measured_col_map (dict): A dictionary mapping cell type names to their corresponding
                                 column names in the measured_df.
    
    Returns:
        list: A list of dictionaries, where each dict is a specificity result.
    """
    
    print("\n----- Calculating Cell-Type Specificity -----")
    
    predictions_df = None
    
    for task_index, single_task_data_dict in enumerate(predictions_file_content["prediction_tasks"]):
        
        if not isinstance(single_task_data_dict, dict):
            print(f"WARNING: Task item at index {task_index} is not a dictionary. Skipping!")
            continue
        # Extract metadata from this task
        task_type_actual = single_task_data_dict.get("type_actual")
        predicted_cell_type = single_task_data_dict.get("cell_type_actual")
        # We also want to extract the cell_type_requested to map it to measured_value_columns_map
        requested_cell_type = single_task_data_dict.get("cell_type_requested")
        predictions_dict =  single_task_data_dict.get("predictions")
        scale_prediction_actual = single_task_data_dict.get("scale_prediction_actual")
        scale_prediction_requested = single_task_data_dict.get("scale_prediction_requested")
        if scale_prediction_actual != scale_prediction_requested:
            print("Predictions scale does not match requested scale (log): Skipping calculation")
            return None

        df = pd.DataFrame(list(predictions_dict.items()), columns=['name', f'prediction_{requested_cell_type}'])
        df[f'prediction_{requested_cell_type}'] = df[f'prediction_{requested_cell_type}'].apply(
            lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
        )
        # Merge into predictions_df
        if predictions_df is None:
            predictions_df = df
        else:
            # Merge on 'name' to ensure alignment
            predictions_df = pd.merge(predictions_df, df, on=seq_id_column, how='inner')
    
    #Check that none of the predicted values have NA in them, if they do skip specificity calculation
    if predictions_df.isna().any().any():  # any column, any row
        print("NA values were found in the predictions, skipping specificity evaluation")
        return None

    merged_df = pd.merge(measured_df, predictions_df, on="name", how="left")
    print(merged_df)
    final_df = merged_df.dropna().copy()
    print(final_df)
    #check size of final data frame here after any NA rows have been dropped for the Measured values
    #same thing here if the Predictions have NA might be a problem
    evaluation_output = pd.DataFrame()
    
    for cell1, cell2 in cell_type_pairs:
        measured_col1 = measured_col_map[cell1]

        measured_col2 = measured_col_map[cell2]

        final_df[f'{cell1}-{cell2}_measured'] = final_df[measured_col1] - final_df[measured_col2]
        final_df[f'{cell1}-{cell2}_predicted'] = final_df[f'prediction_{cell1}'] - final_df[f'prediction_{cell2}']

        r, _ = pearsonr(final_df[f'{cell1}-{cell2}_measured'], final_df[f'{cell1}-{cell2}_predicted'])
        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")

        prediction_task_data_nopredictions = {}
        for task in predictions_file_content["prediction_tasks"]:
            cell_type = task["cell_type_requested"]
            # create an entry only if you want to remove predictions
            prediction_task_data_nopredictions[cell_type] = [
                {k: v for k, v in task.items() if k != "predictions"}
            ]
        #print(prediction_task_data_nopredictions)
        new_row = pd.DataFrame([{
            'Evaluator': evaluator_name,
            "Description": f'Cell type specific expression ({cell1} - {cell2})',
            'Predictor_name': predictor_name,
            'Metric': 'pearson_r',
            'Value': str(r),
            'Prediction_task(s)_data': str(prediction_task_data_nopredictions[cell1]) + " - " + str(prediction_task_data_nopredictions[cell2])

        }])
        evaluation_output = pd.concat([evaluation_output, new_row], ignore_index=True)
    return evaluation_output

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
    correlation_filepath = os.path.join(output_dir, f"correlation_summary_{EVALUATOR_NAME}.csv")
    specificity_filepath = os.path.join(output_dir, f"cell_type_specific_expression_{EVALUATOR_NAME}.csv")
    # Initialize an empty list to get summary for all tasks
    all_task_correlation_results = []
    # Initialize cell-type specificity result container
    specificity_results_df = None
    
    try:
        try:
            # Load measured data file and predictions file ONCE (not with every function call).
            # NOTE: Evaluator builders: If measured_file_path is not a tab-separated file,
            # this line (pd.read_csv) will need to be adjusted or replaced with the
            # appropriate pandas read function (e.g., pd.read_excel, pd.read_csv with different sep)
            # or custom loading logic (e.g., for .npy files).
            measured_df = pd.read_excel(MEASURED_DATA_PATH, header=0)
            
            # Now load predictions
            with open(saved_predictions_path, 'r') as f:
                predictions_file_content = json.load(f)
            
            # Extract Predictor Name
            predictor_name_base = predictions_file_content.get("predictor_name", None) # Resort to None if predictor name is not available
            predictor_name = predictor_name_base.replace(" ", "_").replace("/", "_")

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
            cell_type_pairs_for_specificity = [("human hepatocytes (HepG2)", "induced pluripotent stem cells (iPS cells; WTC11)"), ("human hepatocytes (HepG2)", "lymphoblasts (K562)"), ("induced pluripotent stem cells (iPS cells; WTC11)", "lymphoblasts (K562)")]
            
            chromosome_column = None # If provided
            chromosomes_to_filter_list = None 
            
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
                # Calculate the correlation for each task seperately
                for task_index, single_task_data_dict in enumerate(predictions_file_content["prediction_tasks"]):
                    if not isinstance(single_task_data_dict, dict):
                        print(f"WARNING: Task item at index {task_index} is not a dictionary. Skipping!")
                        continue
                    
                    # Extract metadata from this task
                    task_type_actual = single_task_data_dict.get("type_actual")
                    predicted_cell_type = single_task_data_dict.get("cell_type_actual")
                    # We also want to extract the cell_type_requested to map it to measured_value_columns_map
                    requested_cell_type = single_task_data_dict.get("cell_type_requested")

                    # Find the correspoding measured data column from the map
                    measured_col_for_task = measured_value_columns_map.get(requested_cell_type)
                    
                    print(f"\nProcessing task {task_index+1} (Cell Type: {predicted_cell_type}). Correlating against measured column '{measured_col_for_task}'\
                    for requested cell type '{requested_cell_type}'.")
                    prediction_task_data_nopredictions = [{k: v for k, v in single_task_data_dict.items() if k != "predictions"}]

                    # Call the correlation calculation function
                    task_correlation_dict = calculate_task_correlation(
                        measured_df=measured_df,
                        single_task_data=single_task_data_dict,
                        measured_value_column=measured_col_for_task,
                        seq_id_column=seq_column # chromosome_column and chromosomes_to_filter_list can be added as arguments
                    )
                    
                    if task_correlation_dict:
                        pearson_r_value = task_correlation_dict.get('pearson_r')
                        
                        # Get UTC timestamp for predictor_name
                        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")
                        # And append it to the predictor_name
                        predictor_identifier = f"{predictor_name_base}_{task_index}_{timestamp}" if predictor_name_base else f"UnknownPredictor_{task_index}_{timestamp}"
                        description = f"Agarwal Joint MPRA ({requested_cell_type})"
                        all_task_correlation_results.append({
                            "Evaluator": EVALUATOR_NAME,
                            "Description": description,
                            "Predictor_name": predictor_name,
                            "Time_stamp": timestamp,
                            'Metric': 'pearson_r',
                            'Value': str(pearson_r_value),
                            'Prediction_task(s)_data': prediction_task_data_nopredictions
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
        if all_task_correlation_results:
            summary_df = pd.DataFrame(all_task_correlation_results)
            csv_file_exists: bool = os.path.isfile(correlation_filepath)
            try:
                summary_df.to_csv(correlation_filepath, mode='a',
                                  sep='\t', header=(not csv_file_exists), index=False)
                if csv_file_exists:
                    print("Appended to existing correlation summary CSV file")
                else:
                    print("Created a new correlation summary CSV file")
                print(f"Saved correlation summary to {correlation_filepath}!")
            except IOError as e:
                print("\nNo correlation resuls were saved!")
                
        if specificity_results_df is not None and not specificity_results_df.empty:
            csv_file_exists: bool = os.path.isfile(specificity_filepath)
            try:
                specificity_results_df.to_csv(specificity_filepath, mode='a',
                                              sep='\t', header=(not csv_file_exists), index=False)
                if csv_file_exists:
                    print("Appended to existing cell-type specificity CSV file")
                else:
                    print("Created a new cell-type specificity CSV file")
                print(f"Saved cell-type specificity to {specificity_filepath}!")
            except IOError as e:
                print("\nNo cell-type specificity resuls were saved!")


    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
