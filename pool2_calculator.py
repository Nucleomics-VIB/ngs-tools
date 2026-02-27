#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NGS Pool2 Calculator - Equal Depth Repooling
============================================

PURPOSE:
This script calculates optimal volumes for a second pooling round (pool2) to 
achieve equal target read depth across all samples when combining reads from 
sequencing run1 (pool1) + run2 (pool2). Uses INDIVIDUAL sample reads/volume 
ratios from pool1 for precise per-sample planning.

INPUT CSV FORMAT - REQUIRED COLUMNS (exact names):
sample                    Sample identifier (e.g. 'SampleA', 'CTRL1')
used_volume               Volume used in pool1 in µL (e.g. 2.0)
obtained_read_count       Read count from pool1 sequencing (e.g. 1250000)
volume_left               Remaining volume after pool1 in µL (e.g. 8.0)

OUTPUT CSV COLUMNS:
sample,used_volume,obtained_read_count,volume_left,
pool2_volume,expected_read_count,expected_total_reads,status

USAGE:
    python pool2_calculator.py input.csv output.csv [--target_reads 4000000]

Authors: Stephane Plaisance <stephane.plaisance@vib.be> + Perplexity AI
Affiliation: VIB - Nucleomics Core
Version: 1.0 - Feb 2026 (Status: icon-only, double-warning <80%)
License: GNU GPLv3
"""

import argparse
import pandas as pd
import numpy as np

def main(input_csv, output_csv, target_reads=4000000):
    df = pd.read_csv(input_csv)
    
    required_cols = ['sample', 'used_volume', 'obtained_read_count', 'volume_left']
    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    print("Input data loaded:")
    print(df[required_cols].round(1))
    print()
    
    # Calculate INDIVIDUAL reads per volume ratio for each sample
    df['reads_per_volume'] = df['obtained_read_count'] / df['used_volume']
    
    # Additional reads needed (0 if already at/above target)
    df['additional_reads_needed'] = np.maximum(target_reads - df['obtained_read_count'], 0)
    
    # Volume needed for pool2 based on SAMPLE-SPECIFIC reads/µL
    df['pool2_volume_needed'] = df['additional_reads_needed'] / df['reads_per_volume']
    
    # Clip to available volume, then intelligently round
    df['pool2_volume'] = np.minimum(df['pool2_volume_needed'], df['volume_left'])
    
    for idx in df.index:
        v_needed = df.loc[idx, 'pool2_volume_needed']
        v_available = df.loc[idx, 'volume_left']
        v_calc = df.loc[idx, 'pool2_volume']
        rpv = df.loc[idx, 'reads_per_volume']  # Sample-specific ratio
        
        v_rounded = round(v_calc, 1)
        
        reads_after_rounding = v_rounded * rpv
        additional_needed = df.loc[idx, 'additional_reads_needed']
        
        if (reads_after_rounding < additional_needed and 
            v_rounded < v_available):
            v_rounded = min(v_available, v_rounded + 0.1)
        
        df.loc[idx, 'pool2_volume'] = v_rounded
    
    # Calculate expected reads from final adjusted volumes (sample-specific)
    df['expected_read_count'] = (df['pool2_volume'] * df['reads_per_volume']).round().astype(int)
    df['expected_total_reads'] = df['obtained_read_count'] + df['expected_read_count']
    
    # Status column: icon-only, double warning if <80% target
    df['status'] = '✅'
    df.loc[df['pool2_volume_needed'] > df['volume_left'], 'status'] = '⚠️'
    df.loc[df['expected_total_reads'] < target_reads * 0.8, 'status'] = '⚠️⚠️'
    
    # Output columns
    output_cols = ['sample', 'used_volume', 'obtained_read_count', 
                   'volume_left', 'pool2_volume', 'expected_read_count', 
                   'expected_total_reads', 'status']
    
    df_output = df[output_cols].copy()
    df_output['pool2_volume'] = df_output['pool2_volume'].round(1)
    df_output.to_csv(output_csv, index=False)
    
    print(f"✓ Pool2 plan saved: {output_csv}")
    print(f"  Target reads/sample: {target_reads:,}")
    print("\nPool2 Volumes (rounded to 1 decimal, per-sample optimized):")
    print(df_output[['sample', 'volume_left', 'pool2_volume', 'expected_read_count', 'expected_total_reads', 'status']])
    
    # Volume-limited warning
    limited = df[df['pool2_volume_needed'] > df['volume_left']]
    if len(limited) > 0:
        print(f"\n⚠️  {len(limited)} samples limited by remaining volume")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate pool2 volumes for equal sequencing depth using per-sample reads/volume ratios",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pool2_calculator.py pool1_results.csv pool2_plan.csv
  python pool2_calculator.py data.csv plan.csv --target_reads 4000000
        """
    )
    parser.add_argument("input_csv", help="Input CSV with pool1+seq1 results")
    parser.add_argument("output_csv", help="Output CSV with pool2 plan")
    parser.add_argument("--target_reads", type=int, default=4000000, 
                       help="Target total reads per sample (default: 4M)")
    args = parser.parse_args()
    
    main(args.input_csv, args.output_csv, args.target_reads)
