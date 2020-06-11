import os
import sys
import csv
from multiprocessing import Pool

"""Usage:
   python3 lims_file_parse.py <WOID>
"""


def query(path, work_order):

    """__author__ = 'Andrew Emory'
       generate LIMS files
    """

    path = path + '/'
    print("Querying work-order: %s." % work_order)

    illumina_info = "illumina_info -report library_index_summary --format csv -woid %s --incomplete --output-file-name " \
                    "%sillumina_info.csv" % (work_order, path)
    wo_info = "wo_info --report bam_path -wo %s --format csv --output-file-name %swo_info.csv" % (work_order, path)
    limfo = "limfo e sequence-allocation-path --report result_set_search --setup-wo-id %s  --format csv " \
        "--output-file-name %slimfo.csv" % (work_order, path)
    wo_sample = "wo_info --report sample -wo %s --format csv --output-file-name %swo_sample.csv" % (work_order, path)
    wo_billing = "wo_info --report billing -wo %s --format csv --output-file-name %swo_billing.csv" % (work_order, path)

    processes = (illumina_info, wo_info, limfo, wo_sample, wo_billing)
    pool = Pool(processes=5)
    pool.map(os.system, processes)

    check = os.listdir(path)

    if not "illumina_info.csv" in check:
        print("Something went wrong with the illumina_info query.")
        sys.exit()
    elif not "wo_info.csv" in check:
        print("Something went wrong with the wo_info query.")
        sys.exit()
    elif not "limfo.csv" in check:
        print("Something went wrong with the limfo query.")
        sys.exit()
    elif not "wo_sample.csv" in check:
        print("Something went wrong with the wo_sample query.")
        sys.exit()

    m = {"illumina_info.csv": path + "illumina_info.csv", "wo_info.csv": path + "wo_info.csv",
         "limfo.csv": path + "limfo.csv", "wo_sample.csv": path + "wo_sample.csv",
         "wo_billing.csv": path + "wo_billing.csv"}

    return m


def illumina_limfo_file_merge(illumina_file, limfo_file, woid):

    """__author__ = 'Lee Trani'
       Combine illumina_info and limfo query result csv files
    """

    with open(limfo_file, 'r') as limfo, open(illumina_file, 'r') as illumina, open(
            '{}.illumina.limfo.csv'.format(woid), 'w') as outfile:

        limfo_reader = csv.reader(limfo)
        limfo_data = {}
        limfo_header_found = False

        for l in limfo_reader:

            if len(l) == 0 or 'Entity' in l[0] or '-' in l[0]:
                continue

            if not limfo_header_found and '#' in l[0]:
                limfo_header = l
                limfo_header_found = True
                continue

            if limfo_header_found:
                l_dict = dict(zip(limfo_header, l))
                fastq_file = l_dict['Absolute Path String'].split('/')[-1]
                if len(fastq_file) != 0:
                    if l_dict['Full Name'] not in limfo_data:
                        limfo_data[l_dict['Full Name']] = {fastq_file: l_dict}
                    else:
                        limfo_data[l_dict['Full Name']][fastq_file] = l_dict

        illumina_reader = csv.reader(illumina)
        illumina_data = {}
        illumina_header_found = False

        for l in illumina_reader:

            if len(l) == 0 or 'Index' in l[0] or '-' in l[0]:
                continue

            if not illumina_header_found and '#' in l[0]:
                illumina_header = l
                illumina_header_found = True
                continue

            if illumina_header_found:
                l_dict = dict(zip(illumina_header, l))
                illumina_data[l_dict['Library']] = l_dict

        outfile_writer = csv.DictWriter(outfile, fieldnames=illumina_header[1:] + limfo_header[1:],
                                        extrasaction='Ignore')
        outfile_writer.writeheader()

        for library in illumina_data.keys():
            sample_name = library.split('-lib')[0]
            if sample_name in limfo_data:
                for fastq, line_dict in limfo_data[sample_name].items():
                    if illumina_data[library]['Index Sequence'] in fastq:
                        outfile_writer.writerow({**illumina_data[library], **line_dict})
                    if library in fastq:
                        outfile_writer.writerow({**illumina_data[library], **line_dict})

    return '{}.illumina.limfo.csv'.format(woid)


def wo_sample_merge(illumina_limfo_file, wo_sample_file, woid):

    """__author__ = 'Lee Trani'
       Combine illumina_info_limfo merge file and wo_sample.csv
       Header fields have _wo_sample appended
    """

    with open(wo_sample_file, 'r') as wosf, open(illumina_limfo_file, 'r') as ilf, open(
            '{}.illumina.limfo.sample.csv'.format(woid), 'w') as outfile:

        wosf_reader = csv.reader(wosf)
        wosf_data = {}
        wosf_header_found = False

        for l in wosf_reader:

            if len(l) == 0 or 'Sample Overview' in l[0] or '-' in l[0]:
                continue

            if '#' in l[0]:
                wosf_header = ['{}_wo_sample'.format(x) for x in l]
                wosf_header_found = True
                continue

            if wosf_header_found:
                l_dict = dict(zip(wosf_header, l))
                wosf_data[l_dict['Sample Full Name_wo_sample']] = l_dict

        ilf_reader = csv.DictReader(ilf)

        outfile_writer = csv.DictWriter(outfile, fieldnames=ilf_reader.fieldnames + wosf_header[1:],
                                        extrasaction='Ignore')
        outfile_writer.writeheader()

        for l in ilf_reader:
            if l['Full Name'] in wosf_data:
                outfile_writer.writerow({**l, **wosf_data[l['Full Name']]})

    os.remove(illumina_limfo_file)
    return '{}.illumina.limfo.sample.csv'.format(woid)


def wo_info_merge(illumina_limfo_sample_file, wo_info_file, woid):

    """__author__ = 'Lee Trani'
       Combine illumina_info_limfo_sample merge file and wo_info.csv
       Header fields have _wo_info appended
    """

    with open(wo_info_file, 'r') as wif, open(illumina_limfo_sample_file, 'r') as ilsf, open(
            '{}.illumina.limfo.sample.info.csv'.format(woid), 'w') as outfile:

        wif_reader = csv.reader(wif)
        wif_data = {}
        wif_header_found = False

        for l in wif_reader:

            if len(l) == 0 or 'Illumina Sequencing Bam Path' in l[0] or '-' in l[0]:
                continue

            if not wif_header_found and '#' in l[0]:
                wif_header = ['{}_wo_info'.format(x) for x in l]
                wif_header_found = True
                continue

            if wif_header_found:
                l_dict = dict(zip(wif_header, l))
                if l_dict['Library_wo_info'] not in wif_data:
                    wif_data[l_dict['Library_wo_info']] = l_dict
                else:
                    wif_data[l_dict['Library_wo_info']]['FlowCell_wo_info'] = ','.join(
                        [wif_data[l_dict['Library_wo_info']]['FlowCell_wo_info'], l_dict['FlowCell_wo_info']])
                    wif_data[l_dict['Library_wo_info']]['Lane_wo_info'] = ','.join(
                        [wif_data[l_dict['Library_wo_info']]['Lane_wo_info'], l_dict['Lane_wo_info']])

        ilsf_reader = csv.DictReader(ilsf)

        outfile_writer = csv.DictWriter(outfile, fieldnames=ilsf_reader.fieldnames + wif_header[1:],
                                        extrasaction='Ignore')
        outfile_writer.writeheader()

        for l in ilsf_reader:
            if l['Library'] in wif_data:
                outfile_writer.writerow({**l, **wif_data[l['Library']]})

    os.remove(illumina_limfo_sample_file)
    return '{}.illumina.limfo.sample.info.csv'.format(woid)


def main():

    woid = sys.argv[1]
    # generate LIMS files
    file_dict = query(os.getcwd(), woid)
    # combine illumina_info and limfo files
    illumina_limfo_merged_file = illumina_limfo_file_merge(file_dict['illumina_info.csv'], file_dict['limfo.csv'], woid)
    # combine illumina_limfo file and wo_sample files
    illum_limfo_merged_sample_file = wo_sample_merge(illumina_limfo_merged_file, file_dict['wo_sample.csv'], woid)
    # combine illumina_limfo_sample and wo_info files
    merged_file_final = wo_info_merge(illum_limfo_merged_sample_file, file_dict['wo_info.csv'], woid)
    print('File parsing complete:\n{}'.format(merged_file_final))


if __name__ == '__main__':
    main()
