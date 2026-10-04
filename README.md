# TMA core splitter — QPTIFF prototype

Install Python 3.11+ dependencies:

```bash
pip install numpy tifffile imagecodecs zarr scipy scikit-image pillow
```

Detect core candidates first:

```bash
python split_tma_qptiff.py slide.qptiff core_output --channel DAPI
```

Check `core_preview.png` and `cores.csv`. The CSV contains full-resolution,
exclusive-upper-bound crop coordinates, circularity, solidity, area, and flags.
Flags do not discard samples. Set `export=false` to omit debris or false detections.
Correct boxes, merge fragmented-core rows into one enclosing box, or add missing
cores manually. IDs are detection IDs, not patient IDs or guaranteed TMA grid order.

Export all selected cores, preserving every channel:

```bash
python split_tma_qptiff.py slide.qptiff core_output --review-csv core_output/cores.csv --export
```

For an immediate proof of concept, detection and export can run together:

```bash
python split_tma_qptiff.py slide.qptiff core_output --channel DAPI --export
```

## What is preserved

All full-resolution grayscale channel pixels and dtypes, channel XML (including
names, exposure/normalization metadata), TIFF resolution, and physical crop origin.
All channels use the same rectangle. Source images are never altered. Background
is kept: masks locate cores, they do not erase pixels. Output uses lossless LZW,
512-pixel tiles for large images, baseline channel IFDs, an RGB thumbnail,
level-major downsampled channel IFDs, and a crop overview. Each crop gets a new GUID.
Every baseline pixel in every channel is compared with the source before the
`.partial` output is promoted to `.qptiff`.

## Limits and quality control

This is a derived QPTIFF writer prototype, not an Akoya-certified exporter.
Vendor ValidationCode is removed because it cannot authenticate changed images.
Do not expect vendor image-integrity checks to pass. Original XML is otherwise
retained, including acquisition history; the JSON manifest records crop provenance.
Original slide label/macro images are not copied; the overview shows the crop.
Pyramids are rebuilt with nearest-neighbour sampling, not copied from the source.

Supports flat QPI TIFF IFD layouts with one grayscale page per channel/resolution.
It rejects missing metadata, RGB baseline data, inconsistent pyramid channels, and
layouts it cannot identify. Actual PhenoCycler files may require adaptation.
Compatibility with QuPath/Phenochart must be checked on your actual file before
batch use; pixel equality alone does not establish vendor-software compatibility.

Detection uses one named low-resolution channel (DAPI by default), Gaussian
smoothing, Otsu threshold, closing, hole filling, and connected components.
Solidity/circularity describe that processed mask, not untouched tissue anatomy.
Fragmented cores can split; touching cores can merge. Overlapping crop boxes are
flagged: adjust them so a crop cannot include a neighbour. There is no automatic
reconstruction of missing tissue or automatic separation of merged samples.

Tune `--min-area`, `--max-area`, `--close-radius`, `--sigma`, and `--threshold` on the
preview. These are preview-pixel units; `--margin` is in full-resolution pixels.
Shape thresholds (`--min-circularity .65 --min-solidity .85`) are provisional flags,
not validated biological cutoffs. One channel crop is read at a time; TIFF/Zarr
reads intersecting tiles. Small pyramid arrays are held across channels.

Specification reference: Akoya Phenochart 2.2.0 user guide, Appendix B:
https://akoyabio.boldfocus.com/wp-content/uploads/PhenochartUserManual_2_2_0_rev0.pdf
