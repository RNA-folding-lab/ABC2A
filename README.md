# ABC2A
## Description 
A straightforward and fast method for accurate backmapping of RNA coarse-grained models to all-atom structures
## Requirement
1. gcc
2. shell (option)
## Usage
1. Put the source file ABC2A.c and fragment file into a user-defined directory (e.g., /home/ABC2A/), then compile it using gcc to generate the executable file ABC2A (e.g., gcc -Wall -o3 ABC2A.c -o ABC2A -lm).
2. Create a folder named "RNA_CG" in any working directory (e.g., /home/demo/), and place the coarse-grained model(s) to be reconstructed (in PDB format) into this folder.
3. Copy the shell script "ABC2A.sh" to the working directory(i.e., /home/demo/), and modify the paths of the executable file and fragments within it before saving.
   e.g., export CG_AA_PATH=/home/ABC2A/   &&   fragment_path=$CG_AA_PATH/fragment/
4. Open the terminal, navigate to the working directory, and execute the shell script (e.g., bash ABC2A.sh).
## Outputs
1. "rebuild_AA" (which contains the full atomic reconstructed structures (PDB format) corresponding to the coarse-grained model)
2. "rmsd_detail" (where the RMSD of each nucleotide, reconstructed with different fragments, is stored. This helps us easily identify which fragment was used to reconstruct each nucleotide in the final structure).

Note: The current version of ABC2A is only compatible with the 3-bead coarse-grained model (i.e., P, C4', N1, or N9; see Refs. 1-3). Users can modify the source file (C program) according to their own needs. If assistance is needed, please contact us (Email: yzshi@wtu.edu.cn). Of course, we are currently expanding this method to accommodate any coarse-grained model. Please stay tuned for updates. Thanks for your attention!
## Reference
1. Shi YZ, Wang FH, Wu YY, Tan ZJ. A coarse-grained model with implicit salt for RNAs: predicting 3D structure, stability and salt effect. J Chem Phys. 2014; 141(10):105102. 
2. Shi YZ, Jin L, Feng CJ, Tan YL, Tan ZJ. Predicting 3D structure and stability of RNA pseudoknots in monovalent and divalent ion solutions. PLoS Comput Biol. 2018;14(6):e1006222
3. Wang X, Tan YL, Yu S, Shi YZ, Tan ZJ. Predicting 3D structures and stabilities for complex RNA pseudoknots in ion solu-tions. Biophys J. 2023; 122(8):1503-1516.
4. Perry ZR, Pyle AM, Zhang C. Arena: Rapid and accurate reconstruction of full atomic RNA structures from coarse-grained models. J Mol Biol. 2023; 435(18):168210. 
