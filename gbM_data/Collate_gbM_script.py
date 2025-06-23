import csv
import os

def extract_gbm_genes(file_path, gbm_genes_set):
    """
    Reads a CSV file, extracts gene IDs classified as 'gbM', and adds them to a set.
    It attempts to identify the correct gene ID and classification columns based on common names.

    Args:
        file_path (str): The path to the input CSV file.
        gbm_genes_set (set): A set to store unique gbM gene IDs.
    """
    # Define possible column names for gene ID and classification
    possible_gene_id_cols = ['Gene', 'gene_ID']
    possible_classification_cols = ['Classification', 'Methylation_status']

    try:
        with open(file_path, mode='r', newline='', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)

            # Identify the actual column names from the header of the current file
            header = reader.fieldnames
            gene_id_col = None
            classification_col = None

            # Find the gene ID column
            for col in possible_gene_id_cols:
                if col in header:
                    gene_id_col = col
                    break
            # Find the classification column
            for col in possible_classification_cols:
                if col in header:
                    classification_col = col
                    break

            # Check if required columns were found
            if not gene_id_col:
                print(f"Warning: Could not find a suitable gene ID column (e.g., 'Gene', 'gene_ID') in {file_path}. Skipping file.")
                return
            if not classification_col:
                print(f"Warning: Could not find a suitable classification column (e.g., 'Classification', 'Methylation_status') in {file_path}. Skipping file.")
                return

            print(f"Processing {file_path} using '{gene_id_col}' for gene ID and '{classification_col}' for classification.")

            # Iterate through each row and extract gbM genes
            for row in reader:
                gene_id = row.get(gene_id_col)
                classification = row.get(classification_col)

                # Check if gene_id and classification exist and if classification is 'gbM' (case-insensitive)
                if gene_id and classification and classification.strip().lower() == 'gbm':
                    # Remove any version numbers (e.g., .1, .2) from gene IDs for consistency
                    # Example: "AT1G01040.2" becomes "AT1G01040"
                    gene_id_stripped = gene_id.split('.')[0]
                    gbm_genes_set.add(gene_id_stripped)

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}. Please ensure the path is correct.")
    except Exception as e:
        print(f"An error occurred while processing {file_path}: {e}.")

def create_gbm_gene_list_csv(output_file_path, gbm_genes_set):
    """
    Writes the unique gbM gene IDs to a single-column CSV file.

    Args:
        output_file_path (str): The path to the output CSV file.
        gbm_genes_set (set): A set containing unique gbM gene IDs.
    """
    try:
        with open(output_file_path, mode='w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(['Gene']) # Write the header for the single column

            # Sort the genes alphabetically for consistent output
            sorted_genes = sorted(list(gbm_genes_set))
            for gene in sorted_genes:
                writer.writerow([gene])
        print(f"Successfully created '{output_file_path}' with {len(gbm_genes_set)} unique gbM genes.")
    except Exception as e:
        print(f"An error occurred while writing to {output_file_path}: {e}.")

def main():
    """
    Main function to orchestrate the process of extracting and writing gbM genes.
    """
    # --- Configuration ---
    # The base directory where the CSV files are located.
    base_directory = '/group/sms029/GBM_Data'

    # List of input CSV filenames.
    csv_filenames = ['Bewick_et_al_2016.csv', 'Cahn_et_al_2024.csv']

    # Construct full paths for input files.
    input_files = [os.path.join(base_directory, filename) for filename in csv_filenames] # [4, 5, 6, 7]

    # The name of the output CSV file that will contain the unique gbM genes.
    output_file = 'unique_gbm_genes.csv'
    # --- End Configuration ---

    # Initialize a set to store unique gene IDs
    all_gbm_genes = set()

    # Process each input file
    for file_path in input_files:
        if os.path.exists(file_path): # 8
            extract_gbm_genes(file_path, all_gbm_genes)
        else:
            print(f"Input file not found: {file_path}. Please ensure it exists in the specified location.")

    # Write the collected unique gbM genes to the output CSV file
    if all_gbm_genes:
        # The output file will be created in the current working directory where the script is run
        create_gbm_gene_list_csv(output_file, all_gbm_genes)
    else:
        print("No gbM genes were found across all specified input files. No output CSV will be created.")

if __name__ == "__main__":
    main()