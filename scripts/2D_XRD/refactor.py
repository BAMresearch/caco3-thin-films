# -*- coding: utf-8 -*-
"""
Code Refactoring Helper
=======================
Author: Tomasz Stawski, tomasz.stawski@bam.de
Date: 2026-06-19
Version: 0.1

Description:
A utility script to refactor and clean up the 2D-XRD processing code structure.
"""
import os

with open('processing.py', 'r') as f:
    lines = f.readlines()

# 1. Identify where imports end and argparse starts
start_argparse = 0
for i, line in enumerate(lines):
    if line.startswith('parser = argparse.ArgumentParser'):
        start_argparse = i - 2 # including `import argparse`
        break

header_lines = lines[:start_argparse]

# 2. Extract load_reference_peaks
load_ref_start = 0
load_ref_end = 0
for i, line in enumerate(lines):
    if line.strip().startswith('def load_reference_peaks('):
        load_ref_start = i
    if load_ref_start > 0 and line.strip() == 'return ref_peaks':
        load_ref_end = i + 1
        break

load_ref_lines = [line.lstrip() for line in lines[load_ref_start:load_ref_end]]
for i in range(1, len(load_ref_lines)):
    load_ref_lines[i] = "    " + load_ref_lines[i] # fix indentation

# 3. Identify core processing start
core_start = 0
for i, line in enumerate(lines):
    if line.startswith('converted_h5 = '):
        core_start = i
        break

core_lines = lines[core_start:]

# Remove the inline load_reference_peaks from core_lines
new_core_lines = []
skip = False
for line in core_lines:
    if line.strip().startswith('def load_reference_peaks('):
        skip = True
    if skip and line.strip() == 'return ref_peaks':
        skip = False
        continue
    if not skip:
        new_core_lines.append(line)

core_lines = new_core_lines

# Replace plt.show() with plt.close()
for i in range(len(core_lines)):
    if 'plt.show()' in core_lines[i]:
        core_lines[i] = core_lines[i].replace('plt.show()', 'plt.close()')

# Indent core lines
for i in range(len(core_lines)):
    if core_lines[i] != '\n':
        core_lines[i] = '    ' + core_lines[i]

new_script = "".join(header_lines) + """
import argparse
from pathlib import Path
import glob

def load_reference_peaks(filepath, min_intensity=2.0):
""" + "".join(load_ref_lines[1:]) + """

def process_single_file(input_gfrm, timestamp, calcite_refs, vaterite_refs):
""" + "".join(core_lines)

# We need to return the data we want to compare.
# Find the end of core_lines
# Replace print(...) at the end with a return statement.
# We also need to get the sample_name, theta_lin, corrected_profile, integrated_profile, phase_result, orientation_result.
# Let's add the return statement at the very end of process_single_file.
new_script += """
    return {
        'sample_name': sample_name,
        'theta_lin': theta_lin,
        'corrected_profile': corrected_profile,
        'integrated_profile': integrated_profile,
        'phase_result': phase_result,
        'orientation_result': orientation_result
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process .gfrm files.')
    parser.add_argument('input_path', nargs='?', default='dataXRD', help='Path to .gfrm file or directory')
    args = parser.parse_args()

    input_path = Path(args.input_path)
    if input_path.is_file():
        gfrm_files = [str(input_path)]
    else:
        gfrm_files = [str(p) for p in input_path.rglob('*.gfrm')]

    if not gfrm_files:
        print(f"No .gfrm files found in {args.input_path}")
        exit(0)

    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    crystal_data_dir = os.path.join(script_dir, 'crystaldata')
    calcite_refs = load_reference_peaks(os.path.join(crystal_data_dir, 'Calcite__0000985.txt'), min_intensity=5.0)
    vaterite_refs = load_reference_peaks(os.path.join(crystal_data_dir, 'Vaterite__0004854.txt'), min_intensity=5.0)

    comparison_data = []

    print(f"Found {len(gfrm_files)} .gfrm files. Starting batch processing...")
    for f in gfrm_files:
        print(f"Processing {f}...")
        try:
            data = process_single_file(f, timestamp, calcite_refs, vaterite_refs)
            comparison_data.append(data)
        except Exception as e:
            print(f"Error processing {f}: {e}")

    if comparison_data:
        comp_dir = os.path.join("output", f"comparison_{timestamp}")
        os.makedirs(comp_dir, exist_ok=True)
        
        # Plot comparison of all corrected profiles
        plt.figure(figsize=(12, 8))
        for d in comparison_data:
            plt.plot(d['theta_lin'], d['corrected_profile'], label=d['sample_name'], alpha=0.7)
        
        plt.xlabel('2θ (degrees)')
        plt.ylabel('Intensity')
        plt.title('Comparison of Corrected Profiles')
        # Only show legend if not too many
        if len(comparison_data) <= 15:
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(comp_dir, 'comparison_profiles.png'), dpi=150)
        plt.close()

        # Save comparison data to excel
        comp_excel = os.path.join(comp_dir, 'comparison_summary.xlsx')
        with pd.ExcelWriter(comp_excel, engine='xlsxwriter') as writer:
            # Summary sheet
            summary_df = pd.DataFrame([
                {'Sample Name': d['sample_name'], 
                 'Detected Phase': d['phase_result'], 
                 'Orientation': d['orientation_result']} 
                 for d in comparison_data
            ])
            summary_df.to_excel(writer, sheet_name='Summary', index=False)

            # Combined Corrected Profiles
            profiles_df = pd.DataFrame({d['sample_name']: d['corrected_profile'] for d in comparison_data})
            # Assume theta_lin is same for all (or very similar), use the first one as index
            profiles_df.insert(0, '2Theta', comparison_data[0]['theta_lin'])
            profiles_df.to_excel(writer, sheet_name='Corrected_Profiles', index=False)

        print(f"\\nBatch processing complete! Comparison saved to: {comp_dir}")
"""

with open('processing_new.py', 'w') as f:
    f.write(new_script)
