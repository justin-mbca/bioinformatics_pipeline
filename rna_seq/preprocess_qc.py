import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# Usage: python preprocess_qc.py <input_counts.csv> <output_counts.csv> <output_plot.png>
if __name__ == "__main__":
	# Use the global 'snakemake' variable if present (set by Snakemake), else fallback to CLI args
	if 'snakemake' in globals():
		in_file = snakemake.input[0]
		out_counts = snakemake.output[0]
		out_plot = snakemake.output[1]
	else:
		if len(sys.argv) != 4:
			print("Usage: python preprocess_qc.py <input_counts.csv> <output_counts.csv> <output_plot.png>")
			sys.exit(1)
		in_file, out_counts, out_plot = sys.argv[1:4]

	df = pd.read_csv(in_file, index_col=0)
	print('Count matrix shape:', df.shape)
	print('First 5 rows:')
	print(df.head())

	# Summary statistics
	print('\nSummary statistics:')
	print(df.describe())

	# Check for missing values
	print('\nMissing values per sample:')
	print(df.isnull().sum())

	# Library size per sample
	library_sizes = df.sum(axis=0)
	print('\nLibrary sizes:')
	print(library_sizes)

	# Plot library sizes
	plt.figure(figsize=(6,4))
	sns.barplot(x=library_sizes.index, y=library_sizes.values, palette='Blues_d')
	plt.ylabel('Total Counts')
	plt.title('Library Size per Sample')
	plt.tight_layout()
	plt.savefig(out_plot)
	# plt.show()  # Commented out for non-interactive runs

	# Save cleaned matrix for DESeq2
	cleaned = df.fillna(0).astype(int)
	cleaned.to_csv(out_counts)
	print(f'\nCleaned count matrix saved as {out_counts}')
