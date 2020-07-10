import os
import sys
import csv
from multiprocessing import Pool
import smartsheet
import requests
import glob
import yaml

__author__ = 'Lee Trani'

"""Usage:
   python3 lims_file_parse.py <WOID>
"""


def query(path, work_order):

    """
    __author__ = 'Andrew Emory'
    generate LIMS files
    :param path: directory path
    :param work_order: input woid
    :return: file path dict
    """

    path = path + '/'
    print("Querying work-order: %s." % work_order)

    illumina_info = "illumina_info -report library_index_summary --format csv -woid %s --incomplete " \
                    "--output-file-name %sillumina_info.csv > /dev/null 2>&1 " % (work_order, path)
    wo_info = "wo_info --report bam_path -wo %s --format csv --output-file-name %swo_info.csv > /dev/null 2>&1" \
              % (work_order, path)
    limfo = "limfo e sequence-allocation-path --report result_set_search --setup-wo-id %s  --format csv " \
            "--output-file-name %slimfo.csv > /dev/null 2>&1" % (work_order, path)
    wo_sample = "wo_info --report sample -wo %s --format csv --output-file-name %swo_sample.csv > /dev/null 2>&1" \
                % (work_order, path)
    wo_billing = "wo_info --report billing -wo %s --format csv --output-file-name %swo_billing.csv> /dev/null 2>&1" \
                 % (work_order, path)
    wo_library = "limfo e library-summary --setup-wo-id %s --format csv --output-file-name %swo_library > " \
                 "/dev/null 2>&1" % (work_order, path)

    processes = (illumina_info, wo_info, limfo, wo_sample, wo_billing, wo_library)
    pool = Pool(processes=6)
    pool.map(os.system, processes)

    m = {"illumina_info.csv": path + "illumina_info.csv", "wo_info.csv": path + "wo_info.csv",
         "limfo.csv": path + "limfo.csv", "wo_sample.csv": path + "wo_sample.csv",
         "wo_billing.csv": path + "wo_billing.csv", "wo_library.csv": path + "wo_library.csv"}

    # edited by ltrani
    # sys.exit will only occur if illumina_info or limfo file query fails
    for f, p in m.items():
        if not os.path.isfile(p):
            print("Something went wrong with the {} query.".format(f))
            if f in ['illumina_info.csv', 'limfo.csv']:
                sys.exit('{} file query failed'.format(f))

    return m


def illumina_limfo_file_merge(illumina_file, limfo_file, woid):

    """
    Combine illumina_info and limfo query result csv files
    :param illumina_file: illumina_info file
    :param limfo_file: limfo file
    :param woid: input woid
    """

    with open(limfo_file, 'r') as limfo, open(illumina_file, 'r') as illumina, open(
            '{}.illumina.limfo.merged.csv'.format(woid), 'w') as outfile:

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
                        limfo_data[l_dict['Full Name']] = {l_dict['Absolute Path String']: l_dict}
                    else:
                        limfo_data[l_dict['Full Name']][l_dict['Absolute Path String']] = l_dict

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
            for sample_name in limfo_data.keys():
                for fastq, line_dict in limfo_data[sample_name].items():
                    if illumina_data[library]['Index Sequence'] in fastq:
                        outfile_writer.writerow({**illumina_data[library], **line_dict})
                    if library in fastq:
                        outfile_writer.writerow({**illumina_data[library], **line_dict})

        # file check
        if sum(1 for l in open('{}.illumina.limfo.merged.csv'.format(woid))) == 1:
            sys.exit('illumina_info/limfo merge failed, no matches found.')

    return '{}.illumina.limfo.merged.csv'.format(woid)


def wo_sample_merge(illumina_limfo_file, wo_sample_file):

    """
    Combine illumina_info_limfo merge file and wo_sample.csv
    Header fields have _wo_sample appended
    :param illumina_limfo_file: merged illumina/limfo file
    :param wo_sample_file: wo_sample.csv
    """

    with open(wo_sample_file, 'r') as wosf, open(illumina_limfo_file, 'r') as ilf, open('temp.csv', 'w') as outfile:

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

    os.rename('temp.csv', illumina_limfo_file)


def wo_info_merge(illumina_limfo_sample_file, wo_info_file):

    """
    Combine illumina_info_limfo merge file and wo_info.csv
    Header fields have _wo_info appended
    :param illumina_limfo_sample_file: mmerged illumina/limfo file
    :param wo_info_file: wo_info.csv
    """

    with open(wo_info_file, 'r') as wif, open(illumina_limfo_sample_file, 'r') as ilsf, \
            open('temp.csv', 'w') as outfile:

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

    os.rename('temp.csv', illumina_limfo_sample_file)


def wo_library_merge(illumina_limfo_file, wo_library_file):

    """
    Combine illumina_info_limfo merge file and wo_library.csv
    Header fields have _wo_library appended
    :param illumina_limfo_file: merge illumina/limfo file
    :param wo_library_file: wo_library.csv
    """

    with open(illumina_limfo_file, 'r') as ilf, open(wo_library_file, 'r') as wlf, \
            open('temp.csv', 'w') as outfile:

        wlf_reader = csv.reader(wlf)
        wlf_data = {}
        wlf_header_found = False

        for l in wlf_reader:

            if len(l) == 0 or 'Entity Data for Class' in l[0] or '-' in l[0]:
                continue

            if not wlf_header_found and '#' in l[0]:
                wlf_header = ['{}_wo_library'.format(x.replace('"', '')) for x in l]
                wlf_header_found = True
                continue

            if wlf_header_found:
                l_dict = dict(zip(wlf_header, l))
                wlf_data[l_dict['Full Name_wo_library']] = l_dict

        ilf_reader = csv.DictReader(ilf)

        outfile_writer = csv.DictWriter(outfile, fieldnames=ilf_reader.fieldnames + wlf_header[1:],
                                        extrasaction='Ignore')
        outfile_writer.writeheader()

        for l in ilf_reader:
            if l['Library'] in wlf_data:
                outfile_writer.writerow({**l, **wlf_data[l['Library']]})

    os.rename('temp.csv', illumina_limfo_file)


def get_smartsheet_dt(illumina_limfo_file, woid):

    """
    Get corresponding work order information from Smartsheet
    :param illumina_limfo_file: merged file
    :param woid: woid to query
    :return: exit after matching work order row
    """
    ss_columns = ['If Requires Transfer to GTAC Check Box', 'Work Order ID', 'Method of Transfer',
                  'Manual Demux', 'Analysis/Transfer Instructions', 'Data Transfer Files', 'Data Recipients']

    data = dict.fromkeys(ss_columns, 'NA')

    ss = smartsheet.Smartsheet(os.environ.get('SMRT_API'))
    ss.errors_as_exceptions()

    col_dict = {}
    for col in ss.Sheets.get_columns(5216932677871492).data:
        if col.title in ss_columns:
            col_dict[col.title] = col.id

    for row in ss.Sheets.get_sheet(sheet_id=5216932677871492, column_ids=list(col_dict.values())).rows:

        for cell in row.cells:
            for title, id_ in col_dict.items():
                if cell.column_id == id_:
                    data[title] = cell.value

        if data['Work Order ID'] == woid:

            with open(illumina_limfo_file, 'r') as limfo, open('temp.file.csv', 'w') as outfile:

                limfo_reader = csv.DictReader(limfo)
                outfile_writer = csv.DictWriter(outfile, fieldnames=limfo_reader.fieldnames + list(data.keys()))
                outfile_writer.writeheader()

                for l in limfo_reader:
                    outfile_writer.writerow({**l, **data})

            os.rename('temp.file.csv', illumina_limfo_file)

            attachments = []
            for r in ss.Attachments.list_row_attachments(5216932677871492, row.id).data:
                for c in [' ', '(', ')']:
                    r.name = r.name.replace(c, '')
                open(r.name, 'wb').write(requests.get(ss.Attachments.get_attachment(5216932677871492, r.id).url)
                                         .content)
                attachments.append(r.name)

            if len(attachments) > 0:
                print('Smartsheet attachments downloaded:\n{}'.format('\n'.join(attachments)))

            return


def flowcell_lane_update(illumina_limfo_file):

    with open(illumina_limfo_file, 'r') as limfo, open('temp.file.csv', 'w') as outfile:

        limfo_reader = csv.DictReader(limfo)
        outfile_writer = csv.DictWriter(outfile, fieldnames=limfo_reader.fieldnames + ['flowcell_yaml', 'lane_yaml'])
        outfile_writer.writeheader()

        for l in limfo_reader:
            l['flowcell_yaml'] = 'NA'
            l['lane_yaml'] = 'NA'
            meta_yaml_file = glob.glob('{}/metadata.*.yaml'.format(os.path.dirname(l['Absolute Path String'])))
            if len(meta_yaml_file) == 1 and os.path.isfile(meta_yaml_file[0]):
                with open(meta_yaml_file[0], 'r') as f:
                    yaml_data = yaml.load(f, Loader=yaml.FullLoader)
                    l['flowcell_yaml'] = yaml_data['index_illumina']['flow_cell_id']
                    l['lane_yaml'] = yaml_data['index_illumina']['lane']
            outfile_writer.writerow(l)

        os.rename('temp.file.csv', illumina_limfo_file)


def main():

    # clean up directory
    process_files = ['illumina_info.csv', 'limfo.csv', 'wo_sample.csv', 'wo_info.csv', 'wo_library.csv']
    for file in process_files:
        if os.path.isfile(file):
            os.remove(file)

    woid = sys.argv[1]

    # generate LIMS files
    file_dict = query(os.getcwd(), woid)

    # combine illumina_info and limfo files
    illumina_limfo_merged_file = illumina_limfo_file_merge(file_dict['illumina_info.csv'], file_dict['limfo.csv'], woid)

    # combine illumina_limfo file and wo_sample files
    if os.path.isfile(file_dict['wo_sample.csv']):
        wo_sample_merge(illumina_limfo_merged_file, file_dict['wo_sample.csv'])

    # combine illumina_limfo file and wo_info files
    if os.path.isfile(file_dict['wo_info.csv']):
        wo_info_merge(illumina_limfo_merged_file, file_dict['wo_info.csv'])

    # combine illumina_limfo file and wo_library files
    if os.path.isfile(file_dict['wo_library.csv']):
        wo_library_merge(illumina_limfo_merged_file, file_dict['wo_library.csv'])

    # update flow cell/lane using meta file
    flowcell_lane_update(illumina_limfo_merged_file)

    # get smartsheet data
    get_smartsheet_dt(illumina_limfo_merged_file, woid)

    print('Merged file:\n{}\nFile parsing complete'.format(illumina_limfo_merged_file))


if __name__ == '__main__':
    main()
