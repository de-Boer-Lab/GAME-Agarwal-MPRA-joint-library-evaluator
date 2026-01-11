#script to make downsampled data from original
import json
import pandas as pd

#Read in measurement file from paper
measured_values = pd.read_excel("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/game_dev/Evaluators/Agarwal_Joint_Library/evaluator_data/2023-03-03628-s5/all_cell_type_measured.xlsx")
print(measured_values)
sequences = pd.read_excel("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/game_dev/Evaluators/Agarwal_Joint_Library/evaluator_data/2023-03-03628-s5/2023-03-03628C-Table_S10-joint_lib_design_56k_measured.xlsx", header=0)
print(sequences)

#Drop any rows that don't have a value for all the cell types
measured_values = measured_values.dropna()
print(measured_values)
#Filter both to only keep sequence that have measurements
filtered_sequences = sequences[sequences["name"].isin(measured_values["name"])]
print(filtered_sequences)
filtered_sequences.to_excel(
    "/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/game_dev/Evaluators/Agarwal_Joint_Library/evaluator_data/2023-03-03628-s5/56k_measured_final_sequence_file.xlsx", 
    index=False  # prevents pandas from writing the row index as a column
)

#Also filter the measured data to only keep measurments that we have the sequences for
filtered_measurements = measured_values[measured_values["name"].isin(sequences["name"])]
print(filtered_measurements)
filtered_measurements.to_excel(
    "/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/game_dev/Evaluators/Agarwal_Joint_Library/evaluator_data/2023-03-03628-s5/56k_measured_file.xlsx", 
    index=False  # prevents pandas from writing the row index as a column
)
has_na = filtered_measurements.isna().any().any()
print(has_na)  # True


has_na = filtered_sequences.isna().any().any()
print(has_na)  # True

print(filtered_sequences.columns)
rows_with_na = filtered_sequences[filtered_sequences["230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)"].isna()]
print(rows_with_na)