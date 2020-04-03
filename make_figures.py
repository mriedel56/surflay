#!/usr/bin/env python3
"""
Based on
https://github.com/BIDS-Apps/example/blob/aa0d4808974d79c9fbe54d56d3b47bb2cf4e0a0d/run.py
"""
import os
import os.path as op
import shutil
import subprocess
import argparse
from glob import glob
from nilearn import surface
from nilearn import plotting
from nilearn import datasets
from custom_colors import color_dict
import nibabel as nib
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


def change_image_width(workdir, prefix, plane, i):
    img = Image.open(op.join(workdir, '{prefix}_{plane}_{cut}.png'.format(prefix=prefix, plane=plane, cut=i)))
    img_bk = Image.open(op.join(workdir, '{prefix}_{plane}_{cut}_bk.png'.format(prefix=prefix, plane=plane, cut=i)))
    imageComponents = img_bk.split()
    rgbImage = Image.new("RGB", img_bk.size, (0,0,0))
    rgbImage.paste(img_bk, mask=imageComponents[3])
    imageBox = rgbImage.getbbox()
    cropped=img.crop(imageBox)
    cropped.save(op.join(workdir, '{prefix}_{plane}_{cut}.png'.format(prefix=prefix, plane=plane, cut=i)))

    img = Image.open(op.join(workdir, '{prefix}_{plane}_{cut}.png'.format(prefix=prefix, plane=plane, cut=i)))

    wpercent = (480/float(img.size[0]))
    hsize = int((float(img.size[1])*float(wpercent)))
    img = img.resize((480,hsize), Image.ANTIALIAS)
    img.save(op.join(workdir, '{prefix}_vol_{plane}_{cut}.png'.format(prefix=prefix, plane=plane, cut=i)))


def make_label_image(workdir, plane, dir, i):
    # create a text file with the slice annotation
    img = Image.new('RGB', (480, 40), color = (255, 255, 255))
    fnt = ImageFont.truetype('/Library/Fonts/Arial.ttf', 25)
    d = ImageDraw.Draw(img)
    d.text((220,10), '{0} = {1}'.format(dir, i), font=fnt, fill=(0, 0, 0))
    img.save(op.join(workdir, '{plane}_label_{cut}.png'.format(plane=plane, cut=i)))


def slicer(plane, slices, img, workdir, prefix, thresh, cmap, dim):
    if plane == "sag":
        dir="x"
    elif plane == "cor":
        dir="y"
    elif plane == "ax":
        dir="z"

    for i in slices:
        plotting.plot_stat_map(
            img,
            output_file=op.join(workdir, '{prefix}_{plane}_{cut}.png'.format(prefix=prefix, plane=plane, cut=i)),
            threshold=thresh,
            draw_cross=False,
            annotate=False,
            colorbar=False,
            cmap=cmap,
            display_mode=dir,
            cut_coords=[int(i)],
            dim=dim,
            black_bg=False)

        plotting.plot_stat_map(
            img,
            output_file=op.join(workdir, '{prefix}_{plane}_{cut}_bk.png'.format(prefix=prefix, plane=plane, cut=i)),
            threshold=thresh,
            draw_cross=False,
            annotate=False,
            colorbar=False,
            cmap=cmap,
            display_mode=dir,
            cut_coords=[int(i)],
            dim=dim,
            black_bg=True)

        # cropping, resizing image
        change_image_width(workdir, prefix, plane, i)

        make_label_image(workdir, plane, dir, i)


def run(command, env={}):
    merged_env = os.environ
    merged_env.update(env)
    process = subprocess.Popen(command, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, shell=True,
                               env=merged_env)
    while True:
        line = process.stdout.readline()
        #line = str(line).encode('utf-8')[:-1]
        line=str(line, 'utf-8')[:-1]
        print(line)
        if line == '' and process.poll() is not None:
            break

    if process.returncode != 0:
        raise Exception("Non zero return code: {0}\n"
                        "{1}\n\n{2}".format(process.returncode, command,
                                            process.stdout.read()))


def get_parser():
    parser = argparse.ArgumentParser(description='This script will generate axials, surface medial and surface lateral view images with the specified overlay.')
    parser.add_argument('-f', '--filename', required=True, dest='filename',
                        help=('Full path to input NIFTI file.'))
    parser.add_argument('-w', '--workdir', required=False, dest='workdir',
                        help='Path to a working directory.')
    parser.add_argument('-x', required=False, dest='sag', nargs='*',
                        help='enter x-coordinates for axial slices; should be MNI152-space x-locations')
    parser.add_argument('-y', required=False, dest='cor', nargs='*',
                        help='enter y-coordinates for axial slices; should be MNI152-space y-locations')
    parser.add_argument('-z', required=False, dest='ax', nargs='*',
                        help='enter z-coordinates for axial slices; should be MNI152-space z-locations')
    parser.add_argument('-a', required=False, dest='autocoord', default=None,
                        help='Will automatically select slices based on (1) peak coordinates, (2) sub-peaks, with axial slices receiving preferential treatment.')
    parser.add_argument('--prefix', required=False, dest='prefix',
                        help='prefix name.')
    parser.add_argument('-o', '--outdir', required=False, dest='outdir',
                        help='Path to output directory.')
    parser.add_argument('-c', '--colormap', required=False, dest='cmap', default='coolwarm',
                        help='Colormap options: Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, black_blue, black_blue_r, black_green, black_green_r, black_pink, black_pink_r, black_purple, black_purple_r, black_red, black_red_r, blue_orange, blue_orange_r, blue_red, blue_red_r, blue_transparent, blue_transparent_full_alpha_range, bone, bone_r, brg, brg_r, brown_blue, brown_blue_r, brown_cyan, brown_cyan_r, bwr, bwr_r, cividis, cividis_r, cold_hot, cold_hot_r, cold_white_hot, cold_white_hot_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, cyan_copper, cyan_copper_r, cyan_orange, cyan_orange_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, green_transparent, green_transparent_full_alpha_range, hot, hot_black_bone, hot_black_bone_r, hot_r, hot_white_bone, hot_white_bone_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_hot, ocean_hot_r, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, purple_blue, purple_blue_r, purple_green, purple_green_r, rainbow, rainbow_r, red_transparent, red_transparent_full_alpha_range, roy_big_bl, seismic, seismic_r, spring, spring_r, summer, summer_r, tab10, tab10_r, tab20, tab20_r, tab20b, tab20b_r, tab20c, tab20c_r, terrain, terrain_r, twilight, twilight_r, twilight_shifted, twilight_shifted_r, videen_style, viridis, viridis_r, winter, winter_r.')
    parser.add_argument('--opacity', required=False, dest='opacity', default=1, type=float,
                        help='opacity (0-1, default 1).')
    parser.add_argument('--dim', required=False, dest='dim', type=float, default=-0.3,
                        help='dim (-2-2, default -0.3).')
    parser.add_argument('--threshold', required=False, dest='thresh', default=None, type=float,
                        help='threshold (default is to automatically identify smallest absolute value number).')
    return parser


def main(argv=None):
    args = get_parser().parse_args(argv)

    if not op.exists(args.filename):
        raise ValueError('Argument "filename" must be an existing file.')

    if args.workdir is None:
        currdir = op.dirname(op.realpath(__file__))
        workdir = op.join(currdir, 'tmp')
    else:
        workdir=args.workdir

    if op.isdir(workdir):
        shutil.rmtree(workdir)
    os.makedirs(workdir)
    if args.prefix is None:
        prefix=op.basename(args.filename).split('.')[0]
    else:
        prefix = args.prefix

    if args.outdir is None:
        outdir = op.dirname(args.filename)
    else:
        if not op.isdir(args.outdir):
            raise ValueError('Argument "outdir" must be an existing directory.')
        else:
            outdir = args.outdir

    if args.dim is not None:
        if args.dim < -2 or args.dim > 2:
            raise ValueError('Argument "dim" must be between -2 and 2.')
        else:
            dim=args.dim


    if args.opacity is not None:
        if args.opacity > 1:
            raise ValueError('Argument "opacity" must be less than 1.')
        else:
            opacity=args.opacity


    colormaps=['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap', 'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r', 'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2', 'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r', 'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn', 'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral', 'Spectral_r', 'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r', 'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'black_blue', 'black_blue_r', 'black_green', 'black_green_r', 'black_pink', 'black_pink_r', 'black_purple', 'black_purple_r', 'black_red', 'black_red_r', 'blue_orange', 'blue_orange_r', 'blue_red', 'blue_red_r', 'blue_transparent', 'blue_transparent_full_alpha_range', 'bone', 'bone_r', 'brg', 'brg_r', 'brown_blue', 'brown_blue_r', 'brown_cyan', 'brown_cyan_r', 'bwr', 'bwr_r', 'cividis', 'cividis_r', 'cold_hot', 'cold_hot_r', 'cold_white_hot', 'cold_white_hot_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r', 'cyan_copper', 'cyan_copper_r', 'cyan_orange', 'cyan_orange_r', 'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar', 'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg', 'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'green_transparent', 'green_transparent_full_alpha_range', 'hot', 'hot_black_bone', 'hot_black_bone_r', 'hot_r', 'hot_white_bone', 'hot_white_bone_r', 'hsv', 'hsv_r', 'inferno', 'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_hot', 'ocean_hot_r', 'ocean_r', 'pink', 'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'purple_blue', 'purple_blue_r', 'purple_green', 'purple_green_r', 'rainbow', 'rainbow_r', 'red_transparent', 'red_transparent_full_alpha_range', 'roy_big_bl', 'seismic', 'seismic_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r', 'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'twilight', 'twilight_r', 'twilight_shifted', 'twilight_shifted_r', 'videen_style', 'viridis', 'viridis_r', 'winter', 'winter_r']


    if args.cmap is not None:
        if args.cmap in colormaps:
            cmap = args.cmap
        elif args.cmap in [x for x in color_dict.keys()]:
            cmap = color_dict[args.cmap]
        else:
            raise ValueError('Argument "colormap" is not a valid colormap.')


    fsaverage = datasets.fetch_surf_fsaverage()
    fs_pial_lh=fsaverage['pial_left']
    fs_pial_rh=fsaverage['pial_right']
    fs_sulc_lh=fsaverage['sulc_left']
    fs_sulc_rh=fsaverage['sulc_right']

    if args.thresh is None:
        img = nib.load(args.filename)
        tmpimg = np.abs(img.get_fdata())
        thresh = np.amin(tmpimg[tmpimg>0])-np.finfo(float).eps
    else:
        thresh = args.thresh - np.finfo(float).eps

    img = nib.load(args.filename)
    imgdata = img.get_fdata()

    maxval = np.amax(np.abs(imgdata))

    if len(imgdata[imgdata > 0]) != 0:
        imgdata[imgdata > 0]  = imgdata[imgdata > 0]/np.amax(imgdata[imgdata > 0])
    if len(imgdata[imgdata < 0]) != 0:
        imgdata[imgdata < 0]  = (imgdata[imgdata < 0]/np.amin(imgdata[imgdata < 0]))*-1
    imgdata = imgdata*maxval
    newimg = nib.Nifti1Image(imgdata, img.affine)

    if args.sag is not None:
        sag=[args.sag][0]
        slicer('sag', sag, newimg, workdir, prefix, thresh, cmap, dim)

    if args.cor is not None:
        cor=[args.cor][0]
        slicer('cor', cor, newimg, workdir, prefix, thresh, cmap, dim)

    if args.ax is not None:
        ax=[args.ax][0]
        slicer('ax', ax, newimg, workdir, prefix, thresh, cmap, dim)

    vol_figs = sorted(glob(op.join(workdir, '{prefix}_vol*.png'.format(prefix=prefix))), reverse=True)
    if vol_figs:
        newimg_w = 0
        newimg_h = []
        for vol_fig in vol_figs:
            img = Image.open(vol_fig)
            newimg_w = newimg_w + img.width
            newimg_h.append(img.height)

        dst = Image.new('RGB', (newimg_w, max(newimg_h)), 'white')
        newimg_w = 0
        for vol_fig in vol_figs:
            img = Image.open(vol_fig)
            dst.paste(img, (newimg_w, (dst.height - img.height) // 2))
            newimg_w = newimg_w + img.width

        dst.save(op.join(workdir, '{prefix}_volfigs.png'.format(prefix=prefix)))

    label_figs = sorted(glob(op.join(workdir, '*label*.png'.format(prefix=prefix))), reverse = True)
    if label_figs:
        newimg_w = 0
        newimg_h = []
        for label_fig in label_figs:
            img = Image.open(label_fig)
            newimg_w = newimg_w + img.width
            newimg_h.append(img.height)

        dst = Image.new('RGB', (newimg_w, max(newimg_h)), 'white')
        newimg_w = 0
        for label_fig in label_figs:
            img = Image.open(label_fig)
            dst.paste(img, (newimg_w, (dst.height - img.height) // 2))
            newimg_w = newimg_w + img.width

        dst.save(op.join(workdir, '{prefix}_labelfigs.png'.format(prefix=prefix)))

    if vol_figs:
        volimg = Image.open(op.join(workdir, '{prefix}_volfigs.png'.format(prefix=prefix)))
        labelimg = Image.open(op.join(workdir, '{prefix}_labelfigs.png'.format(prefix=prefix)))

        dst = Image.new('RGB', (volimg.width, volimg.height + labelimg.height), 'white')
        dst.paste(volimg, (0,0))
        dst.paste(labelimg, (0, volimg.height))
        dst.save(op.join(workdir, '{prefix}_slices.png'.format(prefix=prefix)))


    lh_surf = surface.vol_to_surf(newimg, fs_pial_lh, interpolation='nearest')
    rh_surf = surface.vol_to_surf(newimg, fs_pial_rh, interpolation='nearest')

    lh_surf[(lh_surf<thresh) & (lh_surf > 0)] = 0
    lh_surf[(lh_surf>thresh*-1) & (lh_surf < 0)] = 0
    rh_surf[(rh_surf<thresh) & (rh_surf > 0)] = 0
    rh_surf[(rh_surf>thresh*-1) & (rh_surf < 0)] = 0

    for view in ['medial', 'lateral']:
        plotting.plot_surf_stat_map(fs_pial_lh, lh_surf, hemi='left', cmap=cmap, bg_map=fs_sulc_lh, bg_on_data = True, alpha= opacity, colorbar=False, threshold = thresh, view=view, output_file=op.join(workdir, '{prefix}_surf_left_{view}.png'.format(prefix=prefix, view=view)))

        img = Image.open(op.join(workdir, '{prefix}_surf_left_{view}.png'.format(prefix=prefix, view=view)))
        if view == "lateral":
            cropped = img.crop([135, 70, 545, 400])
        else:
            cropped = img.crop([105, 70, 530, 400])
        cropped.save(op.join(workdir, '{prefix}_surf_left_{view}.png'.format(prefix=prefix, view=view)))

        plotting.plot_surf_stat_map(fs_pial_rh, rh_surf, hemi='right', cmap=cmap, bg_map=fs_sulc_rh, bg_on_data = True, alpha=opacity, colorbar=False, threshold = thresh, view=view, output_file=op.join(workdir, '{prefix}_surf_right_{view}.png'.format(prefix=prefix, view=view)))

        img = Image.open(op.join(workdir, '{prefix}_surf_right_{view}.png'.format(prefix=prefix, view=view)))
        if view == 'lateral':
            cropped = img.crop([105, 70, 530, 400])
        else:
            cropped = img.crop([135, 70, 545, 400])
        cropped.save(op.join(workdir, '{prefix}_surf_right_{view}.png'.format(prefix=prefix, view=view)))


    plotting.plot_surf_stat_map(fs_pial_rh, rh_surf, hemi='right', cmap=cmap, bg_map=fs_sulc_rh, bg_on_data = True, alpha=opacity, colorbar=True, threshold = thresh, view=view, output_file=op.join(workdir, '{prefix}_colorbar.png'.format(prefix=prefix)))

    img = Image.open(op.join(workdir, '{prefix}_colorbar.png'.format(prefix=prefix)))
    cropped = img.crop([565, 70, img.width, 400])
    cropped.save(op.join(workdir, '{prefix}_colorbar.png'.format(prefix=prefix)))

    if vol_figs:
        img = Image.open(op.join(workdir, '{prefix}_slices.png'.format(prefix=prefix)))
        hpercent = (360/float(img.size[1]))
        wsize = int((float(img.size[0])*float(hpercent)))
        img = img.resize((wsize, 360), Image.ANTIALIAS)
        img.save(op.join(workdir, '{prefix}_slices.png'.format(prefix=prefix)))

    left_lat = Image.open(op.join(workdir, '{prefix}_surf_left_lateral.png'.format(prefix=prefix)))
    right_lat = Image.open(op.join(workdir, '{prefix}_surf_right_lateral.png'.format(prefix=prefix)))
    left_med = Image.open(op.join(workdir, '{prefix}_surf_left_medial.png'.format(prefix=prefix)))
    right_med = Image.open(op.join(workdir, '{prefix}_surf_right_medial.png'.format(prefix=prefix)))
    colorbar = Image.open(op.join(workdir, '{prefix}_colorbar.png'.format(prefix=prefix)))
    newwidth = left_lat.width + right_lat.width + left_med.width + right_med.width + colorbar.width
    newheight = 330
    if vol_figs:
        sliceimg = Image.open(op.join(workdir, '{prefix}_slices.png'.format(prefix=prefix)))
        newheight = sliceimg.height
        newwidth = newwidth + sliceimg.width

    finalimg = Image.new('RGB', (newwidth, newheight), 'white')
    finalimg.paste(left_lat, (0,0))
    finalimg.paste(left_med, (left_lat.width, 0))
    if vol_figs:
        finalimg.paste(sliceimg, (left_lat.width + left_med.width, 0))
        finalimg.paste(right_med, (left_lat.width + left_med.width + sliceimg.width, 0))
        finalimg.paste(right_lat, (left_lat.width + left_med.width + sliceimg.width + right_med.width, 0))
        finalimg.paste(colorbar, (left_lat.width + left_med.width + sliceimg.width + right_med.width + right_lat.width, 0))
    else:
        finalimg.paste(right_med, (left_lat.width + left_med.width, 0))
        finalimg.paste(right_lat, (left_lat.width + left_med.width + right_med.width, 0))
        finalimg.paste(colorbar, (left_lat.width + left_med.width + right_med.width + right_lat.width, 0))

    finalimg.save(op.join(outdir, '{prefix}.png'.format(prefix=prefix)))

    shutil.rmtree(workdir)

if __name__ == '__main__':
    main()
