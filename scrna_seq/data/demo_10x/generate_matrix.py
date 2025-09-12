import numpy as np
import random

num_genes = 500
num_cells = 50
sparsity = 0.92  # ~92% zeros
max_count = 10

entries = []
for gene in range(1, num_genes + 1):
    for cell in range(1, num_cells + 1):
        if random.random() > sparsity:
            count = random.randint(1, max_count)
            entries.append((gene, cell, count))

with open("matrix.mtx", "w") as f:
    f.write("%%MatrixMarket matrix coordinate integer general\n%\n")
    f.write(f"{num_genes} {num_cells} {len(entries)}\n")
    for gene, cell, count in entries:
        f.write(f"{gene} {cell} {count}\n")

print(f"Wrote matrix.mtx with {len(entries)} nonzero entries.")
