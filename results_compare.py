#!/usr/bin/env python3

import os
import argparse
import filecmp

def get_text_and_tsv_files(folder_path):
    # 1. Initialize an empty list to store the files we find
    matching_files = []
    
    # 2. Get a list of every item inside the directory
    try:
        all_items = os.listdir(folder_path)
    except FileNotFoundError:
        return "Folder not found. Please check the path."

    # 3. Loop through every item one by one
    for item_name in all_items:
        
        # 4. Check if the item is a file (and not a sub-folder)
        # We join the folder path and item name to get the full path
        full_path = os.path.join(folder_path, item_name)
        
        if os.path.isfile(full_path):
            
            # 5. Check the file extension
            # We convert to lowercase so we don't miss .TXT or .TSV
            if item_name.lower().endswith('.txt'):
                matching_files.append(item_name)
                
            elif item_name.lower().endswith('.tsv'):
                matching_files.append(item_name)

    # 6. Return the final list
    return matching_files

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--folders', 
        nargs=2, 
        metavar=('FOLDER1', 'FOLDER2'),
        required=True,
        help='Specify the paths to two folders'
    )
    
    args = parser.parse_args()

    #get list of folder contents excluding the metadata file
    contents_A = [f for f in get_text_and_tsv_files(args.folders[0]) if "metadata" not in f.lower()]
    contents_B = [f for f in get_text_and_tsv_files(args.folders[1]) if "metadata" not in f.lower()]

    print(f"Filenames of both folders are the same (excluding metadata): {contents_A == contents_B}")

    for file in contents_A:
        path_a = os.path.join(args.folders[0], file)
        path_b = os.path.join(args.folders[1], file)

        if file in contents_B:
            if not filecmp.cmp(path_a, path_b, shallow=False):
                print(f"Discrepancy in content for: {file}")
            else:
                print(f"Match for: {file}")
        else:
            print(f"Skipping comparison: {file} not found in second folder.")

    # if args.folders:
    #     folder_one = args.folders[0]
    #     folder_two = args.folders[1]
        
    #     print(f"First folder: {folder_one}")
    #     print(f"Second folder: {folder_two}")
    # else:
    #     print("No folders provided.")
    
    

if __name__ == "__main__":
    main()

