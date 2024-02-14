#!/bin/sh
## Change $CG_AA_PATH to the directories containing the binary executables (ABC2A)
export CG_AA_PATH=/Users/s3000762/Documents/Asia/ABC2A/Program
fragment_path=$CG_AA_PATH/fragment/
file_folder="RNA_CG"
echo "$file_folder"
mkdir rebuild_AA/ 
mkdir rmsd_detail
cd $file_folder
for file in "."/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        if [[ "$filename" == *.pdb ]]; then
           RNA_CG="${filename%.pdb}"
           # RNA_CG="${RNA_CG_1%CG_}"
           echo "$RNA_CG"
           $CG_AA_PATH/ABC2A $RNA_CG $fragment_path
           mv ${RNA_CG}_AA.pdb ../rebuild_AA/
           mv ${RNA_CG}_rmsd.txt ../rmsd_detail/
        fi
    fi
done


