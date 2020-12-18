import os
import shutil
import logging
import atexit
import re
import argparse

from mpi4py import MPI
from mpi4py.futures import MPIPoolExecutor


def get_files():
    """

    :return: all the files to be concatenated
    """
    return ordered_files(['data/{}'.format(f) for f in os.listdir('data')])

def ordered_files(vcfs):
    """
    Make sure the vcf list is in the right order
    :param vcfs:
    :return:
    """
    vcfs.sort(key=lambda f: int(re.sub('\D', '', f)))
    return vcfs

def header_size(vcf_file):
    """

    :param vcf_file: A vcf file
    :return: the size of its header in byte
    """
    the_size = 0
    with open(vcf_file, 'rb') as fp:
        for line in fp:
            if not line.startswith(b'##'):
                break
            the_size = fp.tell()
    logger.debug('AAAAAAAAAAA {} '.format(the_size))
    return the_size


def send(offset, vcf):
    """
    send to worker
    """

class MpiHandle(object):
    FH = None
    output_file = None

    def __init__(self, output_file):
        self.open(output_file)

    def Write_at(self, offset, data, write_all=True):
        # if write_all:
            self.FH.Write_at_all(offset, data)
        # else:
        #     self.FH.Write_at(offset, data)

    @classmethod
    def open(cls, output_file, write_all=True):
        cls.output_file = output_file
        if not cls.FH:
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            amode = MPI.MODE_WRONLY
            logger.info('opening {}'.format(output_file))
            cls.FH = MPI.File.Open(comm, cls.output_file, amode)
            logger.debug('done for {} , rank = {}'.format(output_file, rank))
            atexit.register(exit_handler)

    @classmethod
    def close(cls):
        if cls.FH is not None:
            logger.debug('closing {}'.format(cls.output_file))
            cls.FH.Close()
            cls.FH = None
        else:
            logger.debug('{} already closed'.format(cls.output_file))


def exit_handler():
    MpiHandle.close()


def serf_write(input_file, output_file, header, offset, buffer_size=None, write_all=True):
    """
    Skip the header of input file and write only its data to the output file
    the input file can be large
    :param input_file: The input file with a header's path
    :param output_file: The output file path
    :param header: The size of the header in byte
    :param offset: The position in byte at witch to start writing in the output file
    :param buffer_size: The size of the buffer to read the input file (Default = 4 GB)
    :return: None
    """
    rank = comm.Get_rank()
    h = MpiHandle(output_file)

    if buffer_size is None:
        buffer_size = 4*(1024**3)  # 4 GB

    logger.debug('opening {} writing {}, offset={}'.format(output_file,  input_file, offset))
    with open(input_file, 'rb') as fp:
        fp.seek(header)
        while True:
            data = fp.read(buffer_size)
            if not data:
                break
            h.Write_at(offset, data, write_all=write_all)
            offset = offset + buffer_size

    logger.debug('write {}'.format(offset, data, rank))
    # h.close()
    return input_file


def seigneur(output_file, vcfs, comm=None, n_workers=1):
    '''
    Concatenated vcf that have all the same header, only works with uncompress vcf for now
    :vcfs: a ordered (in chr and position) list of vcf files
    :return:
    '''
    logger.info('starting main')
    # Seigneur

    file_num = len(vcfs)
    header = header_size(vcfs[0])
    total_size = sum([os.path.getsize(filename=f) for f in vcfs]) - (file_num-1) * header
    # the master get to write the chunk one with the header

    logger.debug('copy file header and first chunk')
    shutil.copyfile(vcfs[0], output_file)
    logger.info('start serf pool')

    with MPIPoolExecutor(max_workers=n_workers) as exe:
        offset = os.path.getsize(vcfs[0])
        offsets = [offset]
        for vcf in vcfs[1:-1]:
            offset = offset + os.path.getsize(vcf) - header
            offsets.append(offset)

        logger.debug('map call')
        chunks = int((file_num-1)/(n_workers))
        logger.debug(chunks)
        logger.debug(file_num-1)
        # chunks = 1
        logger.debug('Number of chunks {}'.format(chunks))
        for i in range(chunks):
            first = 1 + i*(n_workers)
            last = first + (n_workers)
            # last = len(offsets)
            logger.debug(vcfs[first:last])
            results = exe.map(serf_write, vcfs[first:last], [output_file for i in offsets[first:last]],
                              [header for i in offsets[first:last]],
                              offsets[first:last])
            for r in results:
                logger.info('results {}'.format(r))
            logger.info('Chunk {} written'.format(i+1))
        # and what is left
        logger.debug('last vcf {} '.format(last))
        results = exe.map(serf_write, vcfs[last:], [output_file for i in offsets[last:]],
                          [header for i in offsets[last:]],
                          offsets[last:], [False for i in offsets[last:]])
        for r in results:
            logger.info('results {}'.format(r))
        logger.info('Last Chunk done')


        return results


def main(comm=None):

    parser = argparse.ArgumentParser(description="Concatenate vcf files that have the same header".format(
        file=os.path.basename(
        __file__)))

    parser.add_argument('--file', '-f', nargs='+')
    parser.add_argument('--output', '-o', default='concateated.vcf', required=False)
    parser.add_argument('--n_workers', '-n', default=1, type=int, required=False, help='Number of worker for the ' \
                                                                                      'process')

    parsed = parser.parse_args()
    output = parsed.output
    n_workers = parsed.n_workers

    in_vcfs_files = parsed.file

    vcfs = []
    for in_file in in_vcfs_files:
        vcfs += [line.strip() for line in open(in_file, 'r')]

    seigneur(output_file=output, vcfs=vcfs, comm=comm, n_workers=n_workers)


if __name__ == '__worker__':
    comm = MPI.COMM_WORLD
    logger = logging.getLogger("Serf_{1}_{0}".format(__file__, comm.rank))
    logging.basicConfig(level=logging.DEBUG)
    logger.info('Start'.format(comm.rank))


if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    logger = logging.getLogger("Seigneur {}".format(__file__))
    logging.basicConfig(level=logging.DEBUG)
    main(comm)
