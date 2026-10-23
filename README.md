# TMA Core Splitter — QPTIFF Prototype

A Python prototype for detecting tissue microarray (TMA) cores in QPTIFF images, reviewing candidate regions, and exporting individual cores while preserving the original channel pixel data.

## Installation

Use **Python 3.11 or later**. Install the required dependencies:

```bash
pip install numpy tifffile imagecodecs zarr scipy scikit-image pillow
```

Run commands from the Core-Splitter project directory with the appropriate Python environment activated.

## 1. Detect Core Candidates

Run the detector using a named channel, such as DAPI:

```bash
python split_tma_qptiff.py slide.qptiff core_output --channel DAPI
```

Replace `slide.qptiff` with the path to your input image. Enclose paths in quotation marks if they contain spaces.

### Generated Files

| File | Contents |
|---|---|
| `core_output/cores.csv` | Candidate core coordinates and review information |
| `core_output/core_preview.png` | Preview of the detected candidates |
| `core_output/core_labels.tif` | Labelled detection mask |

Review the CSV and preview before exporting cores.

### Review the CSV

The CSV contains full-resolution, exclusive-upper-bound crop coordinates, circularity, solidity, area, and flags. Flags identify potential issues; they do not automatically discard samples.

Correct the review CSV as needed:

- Set `export=false` for debris or false detections.
- Correct crop coordinates for incorrectly positioned boxes.
- Merge fragmented detections of the same core into one enclosing box, removing redundant rows.
- Add missing cores manually.

**Detection IDs are not patient IDs and do not guarantee TMA grid order.**

## 2. Review and Export Cores

After reviewing `cores.csv`, export the selected cores while preserving all channels:

```bash
python split_tma_qptiff.py slide.qptiff core_output \
  --review-csv core_output/cores.csv \
  --export
```

- `--review-csv` identifies the reviewed CSV.
- `--export` enables core export.

### Existing Output Folders

The script protects an existing review CSV from accidental overwriting. Running detection in a folder that already contains `cores.csv` stops with a `FileExistsError`.

Use a new output folder for a fresh detection run, or use `--review-csv` to export from an existing reviewed CSV.

Do not delete or overwrite a reviewed CSV unless you have backed it up. Existing exported core files are also protected from overwriting.

## 3. Run Detection and Export Together

For a quick proof of concept:

```bash
python split_tma_qptiff.py slide.qptiff core_output --channel DAPI --export
```

This tests the pipeline without manually reviewing the CSV first. Inspect the resulting detections and exported crops before using them for downstream analysis.

For a fresh run, choose an output folder that does not already contain a protected review CSV.

## 4. Export All Biomarker Channels as OME-TIFF

The separate script, `split_tma_all_biomarkers.py`, exports each core as an **OME-TIFF** with explicitly named biomarker channels.

Both scripts preserve all full-resolution channel pixels. The alternative script uses OME-TIFF channel metadata to make the channels identifiable in compatible viewers.

### Detect Candidates

```bash
python split_tma_all_biomarkers.py slide.qptiff core_output --channel DAPI
```

Review `core_output/cores.csv` and `core_output/core_preview.png`.

### Export Reviewed Candidates

```bash
python split_tma_all_biomarkers.py slide.qptiff core_output \
  --review-csv core_output/cores.csv \
  --export
```

### Detect and Export Together

```bash
python split_tma_all_biomarkers.py slide.qptiff core_output_export \
  --channel DAPI \
  --export
```

### Output Behaviour

- `--channel DAPI` selects the channel used for detection only.
- Every full-resolution biomarker channel is exported.
- Channel names, original pixel intensities and data types are retained.
- Channel colours are transferred when the source contains valid RGB metadata.
- Each output is named `core_XXX.ome.tif`.
- Every exported channel is checked against the corresponding source crop.
- The detection preview displays one channel; it is not the multichannel export.

Open the `.ome.tif` files in an OME-aware viewer, such as QuPath or Fiji/Bio-Formats.

**The OME-TIFF exporter does not transfer the complete vendor-specific QPTIFF metadata or generate an image pyramid.** Check compatibility and performance on an actual slide before processing a large batch.

## What the QPTIFF Exporter Preserves

The following applies to `split_tma_qptiff.py`.

### Pixels and Metadata

- Full-resolution grayscale channel pixels and data types.
- Channel XML metadata, including names and exposure/normalization metadata.
- TIFF resolution and physical crop origin.
- Consistent rectangular crop coordinates across all channels.

Source images are never modified. Background pixels are retained: masks locate cores but do not erase pixels.

### Output Structure

The output uses:

- Lossless LZW compression.
- 512-pixel tiles for large images.
- Baseline channel IFDs.
- An RGB thumbnail.
- Level-major downsampled channel IFDs.
- A crop overview.
- A new GUID for each crop.

### Pixel Verification

Before a `.partial` output is promoted to `.qptiff`, every baseline pixel in every channel is compared with the corresponding source region.

## Limits and Quality Control

### Vendor Compatibility

This is a derived QPTIFF writer prototype, **not an Akoya-certified exporter**.

- The vendor `ValidationCode` is removed because it cannot authenticate modified images.
- Vendor image-integrity checks should not be expected to pass.
- Acquisition/channel XML is retained, except for integrity fields and non-baseline `ScanProfile` entries.
- A JSON manifest records crop provenance.
- Original slide label/macro images are not copied; the overview shows the crop instead.
- Pyramids are rebuilt using nearest-neighbour sampling rather than copied from the source.

### Supported Input Layouts

The exporter supports flat QPI TIFF IFD layouts with one grayscale page per channel and resolution.

It rejects:

- Missing required metadata.
- RGB baseline data.
- Inconsistent pyramid channels.
- Layouts it cannot identify.

Actual PhenoCycler files may require adaptation. Compatibility with QuPath and Phenochart must be checked on the actual input files before batch processing.

**Pixel equality alone does not establish compatibility with vendor software.**

## Detection Method

Detection uses one named low-resolution channel, DAPI by default, with:

1. Gaussian smoothing.
2. Otsu thresholding, unless a manual threshold is supplied.
3. Morphological closing.
4. Hole filling.
5. Connected-component analysis.

Solidity and circularity describe the processed detection mask, not untouched tissue anatomy.

### Known Limitations

- Fragmented cores may be detected as multiple candidates.
- Touching cores may merge into one candidate.
- Overlapping crop boxes are flagged and should be corrected so a crop cannot include a neighbouring core.
- Missing tissue cannot be reconstructed automatically.
- Merged samples are not automatically separated.

## Parameter Tuning

Tune the following options using the detection preview:

| Option | Purpose | Units / behaviour |
|---|---|---|
| `--min-area` | Minimum candidate area | Preview pixels |
| `--max-area` | Flag candidates above this area | Preview pixels; `0` disables the flag |
| `--close-radius` | Morphological closing radius | Preview pixels |
| `--sigma` | Gaussian smoothing | Preview pixels |
| `--threshold` | Manual detection threshold | Original preview intensity units |
| `--margin` | Crop margin | Full-resolution pixels |
| `--min-circularity` | Flag candidates below this circularity | Default: `0.65` |
| `--min-solidity` | Flag candidates below this solidity | Default: `0.85` |

The shape thresholds are provisional flags, not validated biological cutoffs.

### Memory and Image Reading

One channel crop is read at a time, and TIFF/Zarr reads intersecting tiles.

For QPTIFF export, small pyramid arrays are held across channels.

## Reference

[Akoya Phenochart 2.2.0 User Guide, Appendix B](https://akoyabio.boldfocus.com/wp-content/uploads/PhenochartUserManual_2_2_0_rev0.pdf)
