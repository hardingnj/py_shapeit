__author__ = 'Nicholas Harding'

import re
import os
import sh
import uuid
import tempfile
from utils import parse_command, create_sh_script, md5_for_file
import subprocess
import yaml
import sys


class ShapeIt():

    @staticmethod
    def _get_version(executable):

        p = subprocess.Popen(args=[executable, '--version'],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, _ = p.communicate()

        p = re.compile("Version : (.+)\n")  # parentheses for capture groups
        m = p.search(out)
        if m:
            return m.group(1)
        else:
            print(out)
            raise Exception('Version not parsed.')

    def __init__(self, outdir, executable=os.getenv('SHAPEIT', 'shapeit'),
                 ligate_exec=os.getenv('LIGATEHAPLOTYPES', 'ligateHAPLOTYPES'),
                 duohmm_exec=os.getenv('DUOHMM', 'duohmm'),
                 version=None):

        self.executable = executable
        self.ligate_bin = ligate_exec
        self.duohmm_bin = duohmm_exec
        if version is None:
            version = self._get_version(self.executable)
        self.version = version
        self.name = 'ShapeIt'
        self.checksum_file = None

        self.run_id = self.name + '_' + str(uuid.uuid4().get_hex().upper()[0:8])

        self.outdir = os.path.join(outdir, self.name, self.version, self.run_id)
        self.basedir = outdir
        self.si_job_list = list()
        self.ligate_script = None
        self.duohmm_script = None
        self.h5_script = None
        self.si_script = None
        self.graph_convert = "{shapeit} -convert --input-graph {input} " + \
            "--output-max {phased} {sample} --aligned"
        self.dirs = {d: os.path.join(self.outdir, d) for d in ('log', 'script')}
        [sh.mkdir(d, '-p') for d in self.dirs.values() if not os.path.isdir(d)]

        self.haplotypes_f = os.path.join(self.outdir, self.run_id +
                                         '_final.haps.gz')
        self.phased_f = os.path.join(self.outdir, self.run_id +
                                     '_final.samples.gz')

        self.h5out = os.path.join(self.outdir, self.run_id + '.h5')

        self.duohmm_haps = None
        self.duohmm_sample = None
        self.duohmm_geno_e = None
        self.duohmm_recomb = None

        self.param_f = os.path.join(self.outdir, self.run_id + '.yaml')
        self.settings = {'run_id': self.run_id, 'version': self.version,
                         'executable': self.executable, 'outdir': self.outdir,
                         'basedir': outdir, 'haps': self.haplotypes_f,
                         'samples': self.phased_f}

    def setup_region_jobs(self, parameters, regions=None, contig=None,
                          vcf_file=None, pirs=None, graphs=False,
                          duohmm=False, sample_file=None):

        """
        Setup the directory structure and scripts for a shapeIt job that
          automatically splits the job into regions and optionally uses PIRs
          and duoHMM.
        :param parameters: List of parameters to be passed to shapeIt. See
         SHAPEIT docs for full details.
        :param regions: A list of tuples of type (start, end) that represents
         the regions the job is split into. Should have at least 5kb overlap,
         so ligateHAPLOTYPES can work
        :param contig: An ID for the contig being phased
        :param vcf_file:
          VCF format 4.0 file. In gzipped or plaintext. IMPORTANT: shape it does
          not support bzgipped format VCF files.
        :param pirs:
          A list of PIR files that contain the phase informative reads for each
          region in the region list.
        :param graphs:
          boolean, whether to produce haplotype graphs for each region.
        :param duohmm:
          boolean, whether to use the duoHMM tool to correct haplotypes using
          pedigree information.
        :param sample_file:
          if duohmm, the sample file containing the pedigree information.
        :return:
        """

        # basically, set up a shapeIT job for each region. This may be several
        # lines long.
        region_dir = os.path.join(self.outdir, 'region_outputs')
        os.mkdir(region_dir)

        parameters = [str(x) for x in parameters] + ['--input-vcf', vcf_file]
        self.checksum_file = vcf_file

        if pirs is not None:
            parameters.insert(0, '-assemble')
            if graphs:
                raise ValueError("PIRs not compatible with graphs")

        hap_files = list()
        for i, region in enumerate(regions):
            start, stop = [str(x) for x in region]
            region_stem = os.path.join(self.outdir, region_dir,
                                       str(i) + '_' + self.run_id)
            haps = region_stem + '.haps.gz'
            hap_files.append(haps)
            samples = region_stem + '.sample.gz'

            pir_string = ''
            if pirs is not None:
                pir_f = pirs.format(start, stop)
                assert os.path.isfile(pir_f)
                pir_string = " ".join(['--input-pir', pir_f])

            files = ['--input-from', start, '--input-to', stop, pir_string]

            if graphs:
                graph_fn = region_stem + '.graph'
                output = ['--output-graph', graph_fn]

                convert_cmd = self.graph_convert.format(shapeit=self.executable,
                                                        input=graph_fn,
                                                        phased=haps,
                                                        sample=samples)
            else:
                output = ['--output-max', haps, samples]
                convert_cmd = ''

            cmd_shape_it = " ".join([self.executable] + parameters +
                                    files + output)

            script_name = os.path.join(self.dirs['script'], str(i) +
                                       '_shapeIt.sh')
            exit_check = "\n".join(["{",
                                    "if [ -f %s.ok ]; then",
                                    "\techo \"Already completed ok!\"",
                                    "\texit 0",
                                    "fi",
                                    "}"])

            create_sh_script(
                filename=script_name,
                commands=['cd ' + region_dir, exit_check % haps,
                          cmd_shape_it, convert_cmd],
                outfile=haps)

            # name, script, mem, dependency
            self.si_job_list.append(script_name)

        # set up ligateHaplotypes
        tmp = tempfile.NamedTemporaryFile(delete=False)
        gunzip = "gunzip -c {0} > {1}"

        self.ligate_script = os.path.join(self.dirs['script'], 'ligatehaps.sh')
        cmd_ligate = [self.ligate_bin, '--vcf', tmp.name, '--chunks',
                      " ".join(hap_files), '--output', self.haplotypes_f,
                      self.phased_f]

        create_sh_script(self.ligate_script,
                         [gunzip.format(vcf_file, tmp.name),
                          'cd ' + self.outdir,
                          " ".join(cmd_ligate),
                          'rm ' + tmp.name],
                         self.haplotypes_f)

        # set up duohmm if specified
        if duohmm:
            tmp_duo = tempfile.NamedTemporaryFile(delete=False)
            self.duohmm_script = os.path.join(self.dirs['script'],
                                              'duohmm.sh')

            duohmm_root = os.path.join(self.outdir, self.run_id + '_duohmm')
            self.duohmm_haps = duohmm_root + '.haps'
            self.duohmm_sample = duohmm_root + '.sample'
            self.duohmm_geno_e = duohmm_root + '.GE.txt'
            self.duohmm_recomb = duohmm_root + '.RC.txt'

            cmd_duohmm = [self.duohmm_bin, '-H', tmp_duo.name,
                          '-O', duohmm_root,
                          '-G', self.duohmm_geno_e,
                          '-R', self.duohmm_recomb]

            create_sh_script(self.duohmm_script,
                             [gunzip.format(self.haplotypes_f,
                                            tmp_duo.name + '.haps'),
                              "cp {0} {1}".format(sample_file,
                                                  tmp_duo.name+'.sample'),
                              'cd ' + self.outdir,
                              " ".join(cmd_duohmm),
                              " ".join(["gzip",
                                        self.duohmm_haps,
                                        self.duohmm_sample,
                                        self.duohmm_geno_e,
                                        self.duohmm_recomb])],
                             outfile=self.duohmm_haps + '.gz')

        # finally run shapeIt to hdf5 script
        self.h5_script = os.path.join(self.dirs['script'], 'convert_h5.sh')
        h5_root = "{0} {1}/bin/shapeit_2_hdf5.py {2} {3} --chr {4} --out {5}"
        if duohmm:
            h5_cmd = h5_root.format(sys.executable, sys.exec_prefix,
                                    self.duohmm_haps + '.gz',
                                    self.duohmm_sample + '.gz',
                                    contig, self.h5out)
        else:
            h5_cmd = h5_root.format(sys.executable, sys.exec_prefix,
                                    self.haplotypes_f, self.phased_f,
                                    contig, self.h5out)

        create_sh_script(self.h5_script, [h5_cmd], self.h5out)

        self.settings['params'] = parse_command(parameters)

    def run_region_qsub(self, si_args, li_args, dm_args=None):
        """
        Run the previously created jobs via qsub. This is designed for use with
        SGE, perhaps a later version will have more general dependency support.
        :param si_args: Qsub arguments passed to shapeIt jobs.
        :param li_args: Qsub arguments passed to the ligateHAPLOTYPE job.
        :param dm_args: Qsub arguments passed to the duohmm job.
        :return: None
        """
        qsub_parameters = ['-S', '/bin/bash',
                           '-j', 'y',
                           '-o', self.dirs['log']]
        # get checksum
        self.settings['checksum'] = md5_for_file(self.checksum_file)
        yaml.dump(self.settings,
                  stream=open(self.param_f, 'w'))

        for job in self.si_job_list:
            print sh.qsub('-N', self.run_id, qsub_parameters, si_args, job)

        print sh.qsub('-hold_jid', self.run_id, '-N', 'ligate' + self.run_id,
                      qsub_parameters + li_args, self.ligate_script)

        if self.duohmm_script is not None:
            print sh.qsub('-hold_jid', 'ligate' + self.run_id,
                          '-N', 'duohmm' + self.run_id,
                          qsub_parameters + dm_args, self.duohmm_script)
            print sh.qsub('-hold_jid', 'duohmm' + self.run_id,
                          '-N', 'si2hdf5' + self.run_id, '-l', 'h_vmem=4G',
                          qsub_parameters, self.h5_script)
        else:
            print sh.qsub('-hold_jid', 'ligate' + self.run_id,
                          '-N', 'si2hdf5' + self.run_id, '-l', 'h_vmem=4G',
                          qsub_parameters, self.h5_script)

    def setup_single_job(self, parameters, vcf_file, duohmm, sample_file=None):
        # basically, set up a shapeIT job running the whole file as one.
        parameters = [str(x) for x in parameters] + ['--input-vcf', vcf_file]
        self.checksum_file = vcf_file

        cmd_shape_it = " ".join([self.executable] + parameters +
                                ['--output-max', self.haplotypes_f,
                                 self.phased_f])

        self.si_script = os.path.join(self.dirs['script'], 'shapeIt.sh')
        create_sh_script(filename=self.si_script,
                         commands=['cd ' + self.outdir, cmd_shape_it],
                         outfile=self.haplotypes_f)

        if duohmm:
            gunzip = "gunzip -c {0} > {1}"
            tmp_duo = tempfile.NamedTemporaryFile(delete=False)
            self.duohmm_script = os.path.join(self.dirs['script'], 'duohmm.sh')

            duohmm_root = os.path.join(self.outdir, self.run_id + '_duohmm')
            self.duohmm_haps = duohmm_root + '.haps'
            self.duohmm_sample = duohmm_root + '.sample'
            self.duohmm_geno_e = duohmm_root + '.GE.txt'
            self.duohmm_recomb = duohmm_root + '.RC.txt'

            cmd_duohmm = [self.duohmm_bin, '-H', tmp_duo.name,
                          '-O', duohmm_root,
                          '-G', self.duohmm_geno_e,
                          '-R', self.duohmm_recomb]

            create_sh_script(self.duohmm_script,
                             [gunzip.format(self.haplotypes_f,
                                            tmp_duo.name+'.haps'),
                              "cp {0} {1}".format(sample_file,
                                                  tmp_duo.name+'.sample'),
                              'cd ' + self.outdir, " ".join(cmd_duohmm),
                              " ".join(["gzip",
                                        self.duohmm_haps,
                                        self.duohmm_sample,
                                        self.duohmm_geno_e,
                                        self.duohmm_recomb])],
                             outfile=self.duohmm_haps + '.gz')
        self.settings['params'] = parse_command(parameters)

    def run_single(self, si_args, dm_args):

        qsub_parameters = ['-S', '/bin/bash',
                           '-j', 'y',
                           '-o', self.dirs['log']]
        # get checksum
        self.settings['checksum'] = md5_for_file(self.checksum_file)
        yaml.dump(self.settings,
                  stream=open(self.param_f, 'w'))

        print sh.qsub('-N', self.run_id, qsub_parameters,
                      si_args, self.si_script)
        if self.duohmm_script is not None:
            print sh.qsub('-hold_jid', self.run_id,
                          '-N', 'duohmm' + self.run_id,
                          qsub_parameters + dm_args, self.duohmm_script)