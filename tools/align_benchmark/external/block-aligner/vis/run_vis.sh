set -e
jupyter nbconvert --to notebook --inplace --execute block_aligner_vis.ipynb
jupyter trust block_aligner_vis.ipynb
