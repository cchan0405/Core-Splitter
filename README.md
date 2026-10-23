TMA Core Splitter — QPTIFF Prototype

A Python prototype for detecting tissue microarray (TMA) cores in QPTIFF images, reviewing candidate regions, and exporting individual cores while preserving the original channel pixel data.

Installation

Use Python 3.11 or later. Install the required dependencies:

pip install numpy tifffile imagecodecs zarr scipy scikit-image pillow


Run commands from the Core-Splitter project directory with the appropriate Python environment activated.

1. Detect core candidates

Run the detector using a named channel, such as DAPI:

python split_tma_qptiff.py slide.qptiff core_output --channel DAPI


Replace slide.qptiff with the path to your input image. Enclose paths in quotation marks if they contain spaces.

The script generates:

core_output/cores.csv — candidate core coordinates and review information.
core_output/core_preview.png — a preview of the detected candidates.

Review both files before exporting cores.

The CSV contains full-resolution, exclusive-upper-bound crop coordinates, circularity, solidity, area, and flags. Flags identify potential issues; they do not automatically discard samples.

Correct the review CSV as needed:

Set export=false for debris or false detections that should not be exported.
Correct crop coordinates for incorrectly positioned boxes.
Merge rows corresponding to fragmented detections of the same core into one enclosing box, removing redundant rows as appropriate.
Add missing cores manually.

Detection IDs are not patient IDs and do not guarantee TMA grid order.

2. Review and export cores

After reviewing cores.csv, export the selected cores while preserving all channels:

python split_tma_qptiff.py slide.qptiff core_output --review-csv core_output/cores.csv --export


The --review-csv option identifies the reviewed CSV. The --export option enables core export.

Important: The script protects an existing review CSV from accidental overwriting. If the output folder already contains cores.csv, the command may stop with a FileExistsError. Use a new output folder for a fresh run, or follow the script's supported review-CSV workflow. Do not delete or overwrite a reviewed CSV unless you have backed it up.

3. Run detection and export together

For a quick proof of concept, detection and export can run in one command:

python split_tma_qptiff.py slide.qptiff core_output --channel DAPI --export


This is suitable when you want to test the pipeline without manually reviewing the CSV first. Inspect the resulting detections and exported crops before using the results for downstream analysis.

For a fresh run, choose an output folder that does not already contain a protected review CSV.

4. Processing all biomarker channels

The project also includes a separate script, split_tma_all_biomarkers.py, for the all-biomarkers workflow.

Example:

python split_tma_all_biomarkers.py slide.qptiff core_output --channel DAPI


To request export:

python split_tma_all_biomarkers.py slide.qptiff core_output_export --channel DAPI --export


Here, DAPI specifies the channel used for core detection. The script reports detected candidate counts and generates review outputs, depending on its implementation and configuration.

Use a new output directory for a fresh run if an existing review CSV is protected. Review the generated CSV and preview before relying on the exported crops.

Note: The exact export behaviour, supported review-CSV options, and handling of all biomarker channels depend on the implementation of split_tma_all_biomarkers.py. Confirm these behaviours against the script before running a large batch.

What is preserved

The QPTIFF exporter is designed to preserve:

Full-resolution grayscale channel pixels and data types.
Channel XML metadata, including names and exposure/normalization metadata.
TIFF resolution and physical crop origin.
Consistent rectangular crop coordinates across all channels.

Source images are never modified. Background pixels are retained: masks locate cores but do not erase pixels.

The output uses lossless LZW compression, 512-pixel tiles for large images, baseline channel IFDs, an RGB thumbnail, level-major downsampled channel IFDs, and a crop overview. Each crop receives a new GUID.

Before a .partial output is promoted to .qptiff, every baseline pixel in every channel is compared with the corresponding source region.

Limits and quality control

This is a derived QPTIFF writer prototype, not an Akoya-certified exporter.

The vendor ValidationCode is removed because it cannot authenticate modified images. Vendor image-integrity checks should not be expected to pass.
Original XML is otherwise retained, including acquisition history. A JSON manifest records crop provenance.
Original slide label/macro images are not copied; the overview shows the crop instead.
Pyramids are rebuilt using nearest-neighbour sampling rather than copied from the source.

The exporter supports flat QPI TIFF IFD layouts with one grayscale page per channel and resolution. It rejects missing metadata, RGB baseline data, inconsistent pyramid channels, and layouts it cannot identify.

Actual PhenoCycler files may require adaptation. Compatibility with QuPath and Phenochart must be checked on the actual input files before batch processing. Pixel equality alone does not establish compatibility with vendor software.

Detection method

Detection uses one named low-resolution channel (DAPI by default), Gaussian smoothing, Otsu thresholding, morphological closing, hole filling, and connected-component analysis.

Solidity and circularity describe the processed detection mask, not untouched tissue anatomy.

Known limitations include:

Fragmented cores may be detected as multiple candidates.
Touching cores may merge into one candidate.
Overlapping crop boxes are flagged and should be corrected so a crop cannot include a neighbouring core.
Missing tissue cannot be reconstructed automatically.
Merged samples are not automatically separated.
Parameter tuning

Tune the following options using the detection preview:

--min-area — minimum candidate area.
--max-area — maximum candidate area.
--close-radius — morphological closing radius.
--sigma — Gaussian smoothing.
--threshold — threshold adjustment.
--margin — crop margin in full-resolution pixels.

Area and morphology parameters are expressed in preview-pixel units where applicable; --margin is in full-resolution pixels.

The shape thresholds --min-circularity .65 and --min-solidity .85 are provisional flags, not validated biological cutoffs.

One channel crop is read at a time, and TIFF/Zarr reads intersecting tiles. Small pyramid arrays are held across channels.

Reference

Akoya Phenochart 2.2.0 User Guide, Appendix B:

https://akoyabio.boldfocus.com/wp-content/uploads/PhenochartUserManual_2_2_0_rev0.pdf
