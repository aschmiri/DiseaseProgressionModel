#! /usr/bin/env python
# print __doc__
import argparse
import os.path
from subprocess import call
from src.common import adni_tools as adni

EXEC_RVIEW = 'rview'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('study', type=str, help='the study, should be ADNI1, ADNI2, or ADNIGO')
    parser.add_argument('viscode', type=str, help='the visit code, e.g. bl, m12, m24, ...')
    parser.add_argument('trans', type=str, help='the transformation model, e.g. ffd, svffd, sym, or ic')
    parser.add_argument('-s', '--spacing', dest='sx', type=str, default='10')
    parser.add_argument('-r', '--rid', type=str, default=None)
    a = parser.parse_args()

    lut = '/vol/medic02/users/cl6311/Neuro_Atlas/config/lut.csv'

    data_folder = os.path.join(adni.data_folder, a.study)
    image_folder = os.path.join(data_folder, 'native/images_unstripped')
    image_files = adni.get_images(image_folder, a.study, a.viscode)

    seg_folder = os.path.join('/vol/medic01/users/aschmidt/projects/Data/ADNI/data', a.study, 'native')
    seg_folder_bl = os.path.join(seg_folder, 'seg_138regions_baseline')
    seg_folder_fu = adni.make_dir(seg_folder, 'seg_138regions_followup_' + a.trans + '_' + a.sx + 'mm')

    print 'Found ' + str(len(image_files)) + ' images:'
    for i in range(len(image_files)):
        image = image_files[i]
        image_base = os.path.basename(image)

        seg = os.path.join(seg_folder_bl, image_base)
        if not os.path.isfile(seg):
            seg = os.path.join(seg_folder_fu, image_base)

        if a.rid == None or image.find('_S_' + a.rid) > 0:
            if not os.path.isfile(seg):
                print 'No segmentation found for', image
            else:
                print '--------------------'
                print 'Image:      ', image
                print 'Segmenation:', seg

                call([EXEC_RVIEW, image, '-seg', seg, '-lut', lut])


if __name__ == '__main__':
    main()
