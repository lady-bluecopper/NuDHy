# Format the data

Run the script `format_data.py` to extract distributions and other
representations for each dataset:

`python format_data.py <dataset_name> <data_path>`

# AME solutions

To obtain the AME solution, you need to define the experience parameters in a
configuration file.
The following command finds the solutions of the AMEs related to configuration file `lyon_exp1_ame.json` in the `dat` directory and it will save the results
in the same directory using pickle format:

`python ame_sol.py --dir dat/ --conf lyon_exp1_ame`