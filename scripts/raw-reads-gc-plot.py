#!/usr/env python

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks


data_file = r"C:\Users\cbe45\OneDrive\Desktop\Masters\raw-reads-gc-plot\dc1-gc-reads-data.tsv"

# Load data (GC content in first column, count in second column, tab-separated)
gc_content = []
counts = []
with open(data_file, 'r') as f:
    for line in f:
        if line.strip():
            parts = line.strip().split('\t')
            gc = float(parts[0])
            count = float(parts[1])
            gc_content.append(gc)
            counts.append(count)

# Sort data by GC content for smooth plotting
sorted_indices = np.argsort(gc_content)
gc_content_sorted = np.array(gc_content)[sorted_indices]
counts_sorted = np.array(counts)[sorted_indices]

# Find peaks (local maxima)
peaks, _ = find_peaks(counts_sorted)

# Find the peak closest to 50% GC
fifty_idx = np.argmin(np.abs(gc_content_sorted[peaks] - 50))
fifty_peak = peaks[fifty_idx]

# Find the peak closest to 17% GC
seventeen_idx = np.argmin(np.abs(gc_content_sorted[peaks] - 17))
seventeen_peak = peaks[seventeen_idx]

# Find the peak closest to 16% GC (move left circle to 16%)
sixteen_idx = np.argmin(np.abs(gc_content_sorted[peaks] - 16))
sixteen_peak = peaks[sixteen_idx]

# Plot as a smooth line
plt.plot(gc_content_sorted, counts_sorted, color='blue')

# Add dashed vertical line at 28% GC content
plt.axvline(x=29, color='black', linestyle='--', linewidth=1, label='28% GC cutoff')

# Add circles at the 50% peak (green), 17% peak and 16% peak (red)
plt.scatter(gc_content_sorted[fifty_peak], counts_sorted[fifty_peak], color='lightgreen', marker='o', s=100, label='Expected GC Peak')
plt.scatter([gc_content_sorted[seventeen_peak], gc_content_sorted[sixteen_peak]],
            [counts_sorted[seventeen_peak], counts_sorted[sixteen_peak]],
            color='red', marker='o', s=100, label='Anomalous AT-rich Peaks')

plt.xlabel('GC Content (%)')
plt.ylabel('Number of Sequences')
plt.title('GC Content Distribution of Raw Reads')
plt.tight_layout()
plt.legend()
plt.show()
