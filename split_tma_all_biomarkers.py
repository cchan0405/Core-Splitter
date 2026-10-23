#!/usr/bin/env python3
"""Detect TMA cores and export ALL biomarker channels as named OME-TIFF planes.
Install: pip install numpy tifffile imagecodecs zarr scipy scikit-image pillow
Run: python split_tma_all_biomarkers.py slide.qptiff core_output --channel DAPI --export
For reviewed boxes: add --review-csv core_output/cores.csv --export
DAPI selects the detection image ONLY; every full-resolution channel is exported.
core_preview.png is a single-channel detection preview, not the biomarker image.
Open core_XXX.ome.tif in an OME-aware viewer such as QuPath or Fiji/Bio-Formats.
Source pixels, intensity values and bit depth are preserved without masking.
Output is OME-TIFF, not vendor QPTIFF; vendor-specific metadata is not transferred.
Input layout restrictions from the original splitter still apply.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import uuid
import xml.etree.ElementTree as ET
import numpy as np
import tifffile as tf
import zarr
from scipy import ndimage as ndi
from skimage import filters, measure, morphology
from PIL import Image, ImageDraw


def xml(page):
    try:
        root = ET.fromstring(page.description)
    except (ET.ParseError, TypeError):
        return None
    return root if root.tag == 'PerkinElmer-QPI-ImageDescription' else None


def inventory(tif):
    groups = {}
    for p in tif.pages:
        root = xml(p)
        if root is not None and root.findtext('ImageType') in ('FullResolution', 'ReducedResolution'):
            if len(p.shape) != 2 or p.samplesperpixel != 1:
                raise ValueError('Only separate grayscale channel IFDs are supported.')
            groups.setdefault(tuple(p.shape), []).append(p)
    if not groups:
        raise ValueError('No QPI grayscale channels found; refusing generic TIFF export.')
    sizes = sorted(groups, key=lambda s: s[0]*s[1], reverse=True)
    full = groups[sizes[0]]
    names = [xml(p).findtext('Name') for p in full]
    if any(n is None for n in names):
        raise ValueError('Missing channel Name metadata; inspect this file layout first.')
    for shape in sizes:
        if [xml(p).findtext('Name') for p in groups[shape]] != names:
            raise ValueError(f'Channel order/count differs at resolution {shape}.')
    return groups, sizes, names


def region(page, box=None):
    # Zarr decodes only intersecting TIFF tiles, not the whole full-resolution slide.
    store = page.aszarr()
    try:
        arr = zarr.open(store, mode='r')
        return np.asarray(arr[:] if box is None else arr[box[1]:box[3], box[0]:box[2]])
    finally:
        store.close()


def display(a):
    a = np.nan_to_num(a.astype(np.float32))
    lo, hi = np.percentile(a, [1, 99.8])
    return np.uint8(np.clip((a-lo)/max(float(hi-lo), 1e-12), 0, 1)*255)


def detect(image, fullshape, args):
    signal = ndi.gaussian_filter(image.astype(float), args.sigma)
    threshold = args.threshold if args.threshold is not None else filters.threshold_otsu(signal)
    mask = signal > threshold
    if args.close_radius:
        mask = morphology.closing(mask, morphology.disk(args.close_radius))
    mask = ndi.binary_fill_holes(mask)
    labels = measure.label(mask)
    sy, sx = fullshape[0]/image.shape[0], fullshape[1]/image.shape[1]
    rows = []
    for p in measure.regionprops(labels):
        if p.area < args.min_area:
            continue
        y0, x0, y1, x1 = p.bbox
        circularity = min(1., 4*math.pi*p.area/max(p.perimeter_crofton**2, 1.))
        reasons = []
        if circularity < args.min_circularity: reasons.append('low_circularity')
        if p.solidity < args.min_solidity: reasons.append('low_solidity')
        if args.max_area and p.area > args.max_area: reasons.append('large_area_possible_merge')
        if x0 == 0 or y0 == 0 or x1 == image.shape[1] or y1 == image.shape[0]: reasons.append('slide_edge')
        box = [max(0, math.floor(x0*sx)-args.margin), max(0, math.floor(y0*sy)-args.margin),
               min(fullshape[1], math.ceil(x1*sx)+args.margin), min(fullshape[0], math.ceil(y1*sy)+args.margin)]
        rows.append(dict(core_id='', x0=box[0], y0=box[1], x1=box[2], y1=box[3],
                         area_preview=int(p.area), circularity=round(circularity,4), solidity=round(p.solidity,4),
                         flags=';'.join(reasons), export='true', mask_label=p.label))
    rows.sort(key=lambda r: (r['y0'], r['x0']))
    for i,r in enumerate(rows,1): r['core_id'] = f'core_{i:03d}'
    for i,r in enumerate(rows):
        for j,s in enumerate(rows):
            if i != j and max(r['x0'],s['x0']) < min(r['x1'],s['x1']) and max(r['y0'],s['y0']) < min(r['y1'],s['y1']):
                if 'overlapping_crop' not in r['flags']: r['flags'] += ';overlapping_crop'
    return rows, labels


def description(page, kind, identifier, baseline=False):
    root = xml(page)
    # Preserve acquisition/channel metadata, never reuse vendor integrity values.
    for tag in ('ValidationCode',):
        for node in root.findall(tag): root.remove(node)
    if not baseline:
        for node in root.findall('ScanProfile'): root.remove(node)
    for tag,value in [('ImageType',kind),('Identifier',identifier)]:
        node = root.find(tag)
        if node is None: node = ET.SubElement(root,tag)
        node.text = value
    if kind in ('Thumbnail','Overview'):
        for tag in ('Name','Color','Responsivity','ExposureTimeArray'):
            for node in root.findall(tag): root.remove(node)
    return ET.tostring(root,encoding='unicode')


def rational(tag):
    val = tag.value
    return val[0]/val[1] if isinstance(val,tuple) else float(val)


def write_page(writer, a, source, kind, identifier, scale=1, baseline=False, box=None):
    rx = rational(source.tags['XResolution'])
    ry = rational(source.tags['YResolution'])
    unit = int(source.tags['ResolutionUnit'].value)
    extra = []
    if box is not None:
        for code,key,offset,res in [(286,'XPosition',box[0],rx),(287,'YPosition',box[1],ry)]:
            pos = rational(source.tags[key]) if key in source.tags else 0.
            if kind == 'Overview': offset += (a.shape[1] if code == 286 else a.shape[0])*scale/2
            extra.append((code,'2I',1,(int(round((pos+offset/res)*1000000)),1000000),False))
    # TIFF sample bounds use the actual sample type.
    dtypecode = {'u':'B' if a.dtype.itemsize == 1 else 'H','i':'h','f':'f'}[a.dtype.kind]
    extra += [(340,dtypecode,1,a.min().item(),False),(341,dtypecode,1,a.max().item(),False)]
    writer.write(a, photometric='rgb' if a.ndim == 3 else 'minisblack',
                 tile=(512,512) if max(a.shape[:2]) > 2048 and kind not in ('Thumbnail','Overview') else None,
                 compression='lzw', metadata=None, description=description(source,kind,identifier,baseline),
                 software='PerkinElmer-QPI (derived crop; split_tma_qptiff prototype)',
                 resolution=(rx/scale,ry/scale), resolutionunit=unit,
                 subfiletype=0 if kind == 'FullResolution' else 1,
                 datetime=source.tags['DateTime'].value if 'DateTime' in source.tags else None,
                 extratags=extra)


def export_core(tif, pages, names, row, dest):
    """Write all source channel crops into one explicitly named OME channel axis."""
    box = tuple(int(row[k]) for k in ('x0', 'y0', 'x1', 'y1'))
    h, w = pages[0].shape
    if not (0 <= box[0] < box[2] <= w and 0 <= box[1] < box[3] <= h):
        raise ValueError(f'Invalid bounds for {row["core_id"]}')
    dtype = pages[0].dtype
    if any(page.dtype != dtype for page in pages):
        raise ValueError('Mixed channel dtypes cannot be preserved in one OME series.')
    target = dest / (row['core_id'] + '.ome.tif')
    tmp = target.with_suffix('.tif.partial')
    if target.exists() or tmp.exists():
        raise FileExistsError(target if target.exists() else tmp)
    metadata = {'axes': 'CYX', 'Channel': {'Name': names},
                'Description': 'Unmodified channel pixels cropped from ' + str(Path(tif.filename).name)
                               + '; source XY bounds: ' + str(box)}
    # QPI Color is an RGB comma-separated triplet. OME Color is signed RGBA.
    colors = []
    for page in pages:
        try:
            rgb = [int(x.strip()) for x in xml(page).findtext('Color').split(',')]
            if len(rgb) != 3 or any(x < 0 or x > 255 for x in rgb):
                raise ValueError('Invalid RGB')
            rgba = (rgb[0] << 24) | (rgb[1] << 16) | (rgb[2] << 8) | 255
            colors.append(rgba if rgba < 2**31 else rgba - 2**32)
        except (AttributeError, ValueError):
            colors = []
            break
    if colors:
        metadata['Channel']['Color'] = colors
    source = pages[0]
    if all(key in source.tags for key in ('XResolution', 'YResolution', 'ResolutionUnit')):
        unit = int(source.tags['ResolutionUnit'].value)
        microns = {2: 25400., 3: 10000.}.get(unit)
        if microns:
            for axis in ('X', 'Y'):
                resolution = rational(source.tags[axis + 'Resolution'])
                if resolution > 0:
                    metadata['PhysicalSize' + axis] = microns / resolution
                    metadata['PhysicalSize' + axis + 'Unit'] = 'µm'
    try:
        with tf.TiffWriter(tmp, bigtiff=True, ome=True) as writer:
            writer.write(data=(region(page, box) for page in pages),
                         shape=(len(pages), box[3]-box[1], box[2]-box[0]),
                         dtype=dtype, photometric='minisblack', compression='lzw',
                         metadata=metadata)
        with tf.TiffFile(tmp) as out:
            root = ET.fromstring(out.ome_metadata)
            channels = root.findall('.//{*}Channel')
            if [ch.get('Name') for ch in channels] != names:
                raise ValueError('OME channel name/count verification failed')
            if len(out.pages) != len(pages):
                raise ValueError('Exported channel page count differs from source')
            for src, exported in zip(pages, out.pages):
                expected, actual = region(src, box), exported.asarray()
                if expected.dtype != actual.dtype or not np.array_equal(expected, actual, equal_nan=True):
                    raise ValueError('Full-resolution pixel verification failed')
        tmp.rename(target)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    return dict(file=target.name, box=list(box), channels=names,
                pixel_verification='all full-resolution pixels equal', format='OME-TIFF')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path); p.add_argument('output',type=Path)
    p.add_argument('--export',action='store_true',help='Export after detection, or from --review-csv')
    p.add_argument('--review-csv',type=Path,help='Edited CSV with full-resolution boxes and export=true/false')
    p.add_argument('--channel',default='DAPI',help='Exact channel name for core detection only')
    p.add_argument('--preview-max',type=int,default=4000)
    p.add_argument('--threshold',type=float,help='Manual threshold in original preview intensity units')
    p.add_argument('--sigma',type=float,default=2)
    p.add_argument('--close-radius',type=int,default=5,help='Preview pixels; may merge nearby fragments')
    p.add_argument('--min-area',type=int,default=200,help='Preview pixels, debris cutoff')
    p.add_argument('--max-area',type=int,default=0,help='Preview pixels; zero disables upper-area flag')
    p.add_argument('--min-circularity',type=float,default=.65)
    p.add_argument('--min-solidity',type=float,default=.85)
    p.add_argument('--margin',type=int,default=64,help='Full-resolution pixels')
    args=p.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    with tf.TiffFile(args.input) as tif:
        groups,sizes,names=inventory(tif); pages=groups[sizes[0]]
        print('Channels:',names,flush=True)
        if args.review_csv:
            with args.review_csv.open() as f: rows=list(csv.DictReader(f))
        else:
            if args.channel not in names: raise ValueError(f'Choose --channel from {names}')
            shape=next((s for s in sizes if max(s)<=args.preview_max),sizes[-1])
            if max(shape)>args.preview_max:
                raise ValueError('No sufficiently small pyramid level; increase --preview-max explicitly.')
            image=region(groups[shape][names.index(args.channel)])
            rows,labels=detect(image,pages[0].shape,args)
            csvpath=args.output/'cores.csv'
            if csvpath.exists(): raise FileExistsError('Use a new output folder or --review-csv; existing review CSV protected.')
            fields=list(rows[0]) if rows else ['core_id','x0','y0','x1','y1','area_preview','circularity','solidity','flags','export','mask_label']
            with csvpath.open('w',newline='') as f:
                wr=csv.DictWriter(f,fieldnames=fields);wr.writeheader();wr.writerows(rows)
            tf.imwrite(args.output/'core_labels.tif',labels.astype(np.uint32),metadata=None)
            preview=Image.fromarray(display(image)).convert('RGB');draw=ImageDraw.Draw(preview)
            sy,sx=pages[0].shape[0]/shape[0],pages[0].shape[1]/shape[1]
            for r in rows:
                b=[r['x0']/sx,r['y0']/sy,r['x1']/sx,r['y1']/sy]
                color='orange' if r['flags'] else 'lime'
                draw.rectangle(b,outline=color,width=2);draw.text((b[0],b[1]),r['core_id'],fill=color)
            preview.save(args.output/'core_preview.png')
            print(f'{len(rows)} candidates. Review cores.csv and core_preview.png.',flush=True)
        if args.export:
            manifest=[]
            for row in rows:
                if row['export'].strip().lower() not in ('true','1','yes'): continue
                core_id=row['core_id']
                if not core_id or Path(core_id).name != core_id or '/' in core_id or '\\' in core_id:
                    raise ValueError('core_id must be a filename-safe ID')
                print('Exporting',core_id,row.get('flags',''),flush=True)
                record=export_core(tif,pages,names,row,args.output)
                record['source']=str(args.input.resolve())
                record['qc_flags']=row.get('flags','')
                manifest.append(record)
                (args.output/'export_manifest.json').write_text(json.dumps(manifest,indent=2))

if __name__=='__main__': main()
