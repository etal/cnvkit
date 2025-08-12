# HMM Segmentation in CNVkit

HMM (Hidden Markov Model) segmentation is an optional feature in CNVkit that provides advanced segmentation methods using machine learning.

## Installation

### With pip
```bash
# Install CNVkit with HMM support
pip install cnvkit[hmm]

# Or install HMM dependencies separately
pip install cnvkit
pip install pomegranate>=1.0.0
```

### With conda
```bash
# Use the HMM-enabled environment
conda env create -f conda-env-hmm.yml
conda activate cnvkit-hmm

# Or add to existing environment
conda install pomegranate>=1.0.0 pytorch>=1.13.0
```

## Usage

HMM segmentation provides three methods:
- `hmm` - 3-state model with flexible means (loss, neutral, gain)
- `hmm-tumor` - 5-state model for tumor samples (deletion, loss, neutral, gain, amplification)  
- `hmm-germline` - 3-state model with fixed means for germline samples

```bash
# Use HMM segmentation
cnvkit.py segment sample.cnr -m hmm -o sample.cns
cnvkit.py segment sample.cnr -m hmm-tumor -o sample.cns
cnvkit.py segment sample.cnr -m hmm-germline -o sample.cns
```

## Dependencies

HMM segmentation requires:
- pomegranate >= 1.0.0 (machine learning library)
- pytorch >= 1.13.0 (deep learning backend)

These are heavy dependencies (~500MB+) which is why HMM functionality is optional.

## Performance Notes

CNVkit v0.9.12+ uses pomegranate v1.0.0 which has a PyTorch backend. While more powerful, this can be slower than previous versions. Performance has been optimized by:
- Reducing max iterations from 100,000 to 100
- Adding early stopping tolerance (1e-4)
- Using float32 tensors for better memory efficiency

If you need the old faster HMM implementation, use CNVkit v0.9.11 or earlier with pomegranate v0.14.8.